##
# MoA - VIKING on RNA-Seq Data
# ----------------------------
# system.time({source("sources/viking.R")})

# BiocManager::install( c("EmpiricalBrownsMethod","doSNOW","qvalue","Rfast") )
# install.packages("../DeMAND2/",type = "source",repos = NULL)
library(DeMAND2)

require(tidyverse)
require(ggplot2)
require(ggrepel)
require(viper)
library(doParallel)
registerDoParallel(cores=4)

source("sources/utils.R")

preppi_str_threshold <- 5
preppi_walk_n_neighbors <- 5

create_workspace(run_dir = "viking-results-for-cell-chemical-biology-paper")

## Loading Data ----
nbl_data <- list()
nbl_data$expmat <- readRDS("data/trueseq-cpm.rds")
nbl_data$metadata <- readRDS("data/trueseq-metadata.rds")
str(nbl_data,1)
dim(nbl_data$expmat)

## Loading Kinases Tables from: http://kinase.com/human/kinome/ ----
kinome_df <- readxl::read_excel( "data/Kincat_Hsap.08.02.xls" , sheet = 1 )
kinases <- sort( kinome_df$Entrez_Symbol[ kinome_df$`Pseudogene?` == "N" ] )
kinases <- kinases[ kinases != "" ]

preppi.demand.network_kinases <- generateDemandNetworkFromPrePPI( regulators = kinases , likelihood.threshold = 0.9 , structural_score.threshold = preppi_str_threshold , n_neighbors = preppi_walk_n_neighbors )
preppi.demand.network_TF <- generateDemandNetworkFromPrePPI( regulators = "TF" , likelihood.threshold = 0.9 , structural_score.threshold = preppi_str_threshold , n_neighbors = preppi_walk_n_neighbors )
preppi.demand.network_coTF <- generateDemandNetworkFromPrePPI( regulators = "coTF" , likelihood.threshold = 0.9 , structural_score.threshold = preppi_str_threshold , n_neighbors = preppi_walk_n_neighbors )
preppi.demand.network_SIG <- generateDemandNetworkFromPrePPI( regulators = "SIG" , likelihood.threshold = 0.9 , structural_score.threshold = preppi_str_threshold , n_neighbors = preppi_walk_n_neighbors ) # summary(preppi.table$str_score)$median is 5.5
preppi.demand.network <- rbind(preppi.demand.network_TF,preppi.demand.network_coTF,preppi.demand.network_SIG,preppi.demand.network_kinases)
nrow(preppi.demand.network)
head(preppi.demand.network)

print_msg_info(">>> Using PrePPI network")
{
  tmp <- as.data.frame( preppi.demand.network )
  tmp$stringForDuplicates <- make.names( paste0(tmp$Gene1,"-",tmp$Gene2) )
  a <- nrow(tmp)
  tmp <- tmp[ !(duplicated(tmp$stringForDuplicates)) , ]
  b <- nrow(tmp)
  message("- Found " , (a-b) , " duplicated edges - removing ...")
  tmp$stringForDuplicates <- NULL # removing temorpary variables used to delete duplicate edges
  
  nbl.demand.network <- as.matrix(tmp)
  rm(tmp)
  nrow(nbl.demand.network) # the merged network  
}

## Loading Gene Regulatory Network (TARGET NBL) ----
print_msg_info(">>> Loading Gene Regulatory Network")
interactome.target_db <- readRDS( "data/target-nbl-network.rds")
print_msg_info(">>> Running VIPER Analysis")
nbl_data$vpmat <- viper( nbl_data$expmat , regulon = pruneRegulon(interactome.target_db,100) , method = "rank" )
print(dim(nbl_data$vpmat))

controls_tags <- c("Null")
print(controls_tags)
drug_names <- unique(nbl_data$metadata$treatment)
drug_names <- setdiff(drug_names,controls_tags)
print(drug_names)
drug_names <- drug_names[ !grepl( "MOCK" , drug_names ) ]
print(drug_names)

nrow(nbl.demand.network)
head(nbl.demand.network)

# -------

print_msg_info(">>> Preparing for DeMAND Analysis ...")
demandAnalysis <- list()
demandAnalysis$demandObject <- list()
demandAnalysis$demandResults <- list()
demandAnalysis$computationTime <- list()

## Run VIKING on several drugs ----
for( aDrug in drug_names )
{
  # aDrug <- "Isop"
  
  ## Running DeMAND ----
  print_msg_info(">>> >> Running DeMAND on drug - " , aDrug )
  {
    demand_treatment_sample_names <- nbl_data$metadata$sample_id[ nbl_data$metadata$treatment %in% aDrug ]
    demand_control_sample_names <- nbl_data$metadata$sample_id[ nbl_data$metadata$treatment %in% controls_tags ]
    
    demand_treatment_samples <- nbl_data$vpmat[,demand_treatment_sample_names]
    demand_control_samples <- nbl_data$vpmat[,demand_control_sample_names]
    
    demand_exp <- cbind( demand_treatment_samples , demand_control_samples )
    print(colnames(demand_exp))
    
    treatment_index <- 1:ncol(demand_treatment_samples)
    control_index <- (ncol(demand_treatment_samples)+1):(ncol(demand_treatment_samples)+ncol(demand_control_samples))
    
    print(treatment_index)
    print(control_index)
    
    print(dim(demand_exp))
    demand_exp <- demand_exp[ rowVars(demand_exp) > 0 , ]
    print(dim(demand_exp))      
    
    demand_annotation <- cbind( rownames(demand_exp) , rownames(demand_exp) )
    
    demandAnalysis$demandObject[[aDrug]] <- demandClass( exp = demand_exp , 
                                                         anno = demand_annotation , 
                                                         network = nbl.demand.network )
    demandAnalysis$computationTime[[aDrug]] <- print( system.time( demandAnalysis$demandResults[[aDrug]] <- runDeMAND( demandAnalysis$demandObject[[aDrug]] , 
                                                                                                                       fgIndex = treatment_index , bgIndex = control_index
    ) ) )
    # }
    
    filename <- file.path(processed_data.dir,paste0("viking-analyisis-on-",aDrug,"-vs-null-analog-from-truseq-data.rds") )
    saveRDS(object = demandAnalysis,file = filename )
    
    df <- demandAnalysis$demandResults[[aDrug]]@moa %>% as_tibble()
    df <- df %>% arrange(adjustedPvalue)
    
    write_csv( df , file.path(processed_data.dir,paste0("demand-on-viper-analyisis-on-",aDrug,"-vs-null-analog-from-truseq.csv")) )
  }
  
}
print_msg_warn("*** [VIKING Analysis Completed] ***")

viking_table <- tibble()
## Plotting VIKING results from several drugs ----
for( aDrug in drug_names )
{
  filename <- file.path(processed_data.dir,paste0("viking-analyisis-on-",aDrug,"-vs-null-analog-from-truseq-data.rds") )
  ## Load VIKING data ----
  print_msg_info(">>> Loading VIKING data and plotting VIKING Scatter Plot for " , aDrug )
  {
    demandAnalysis <- readRDS( file = filename )
    demand_table <- as_tibble(demandAnalysis$demandResults[[aDrug]]@moa)
    demand_table$bonferroni <- p.adjust(demand_table$Pvalue,method = "bonferroni")
    print_msg_warn("*** Using Benjamini Hockberg FDR Correction ***")
    demand_table$FDR <- p.adjust(demand_table$Pvalue,method = "BH")
    index <- grepl(aDrug,colnames(nbl_data$vpmat))
    print_msg_info(">>> >> Stouffer on: " , colnames(nbl_data$vpmat)[index] )
    viper_activity <- nbl_data$vpmat[,index]
    viper_activity[ grepl("CSNK2|BLK|LYN|PIK|TOR", rownames(viper_activity) ) , ]
    y <- doStouffer(viper_activity)
    viper_stouffer_zscore <- y
    print(viper_stouffer_zscore[c("CSNK2A1","MYCN","TEAD4","MYC","FYN","LYN","BLK","CDK2","CDK6")])
    
    demand_table <- cbind( demand_table , viper_score = viper_stouffer_zscore[ match( demand_table$moaGene , names(viper_stouffer_zscore) ) ] )
    demand_table <- demand_table %>% dplyr::filter( !is.na( demand_table$viper_score ) )
    demand_table <- as_tibble(demand_table)
    demand_table
    print( demand_table[ demand_table$moaGene %in% c("CSNK2A1","CSNK2A2","BLK","LYN","FYN") , ] )
    
    ## VIKING Plot ----
    print_msg_info(">>> VIKING Plot for " , aDrug )
    {
      viper_threshold <- -1
      demand_bonferroni_threshold <- 0.05
      demand_table_filtered <- demand_table %>% dplyr::filter( moaGene %in% kinases )
      
      print_msg_warn( ">>> >> VIPER SCORE min: " , min(demand_table$viper_score) )
      print_msg_warn( ">>> >> DeMAND SCORE min: " , min(demand_table$bonferroni) )
      
      demand_min_limit <- min(log10(demand_table$bonferroni))
      viper_min_limit <- min(demand_table$viper_score)
      
      demand_table_filtered %>% 
        filter( bonferroni < demand_bonferroni_threshold ) %>%
        filter( viper_score < viper_threshold )
      
      demand_plot <- ggplot( data = demand_table_filtered , aes( x = log10(bonferroni) , y = viper_score ) ) +
        
        annotate("rect", xmin = -Inf , xmax =  0 , ymin = viper_threshold , ymax = Inf, fill = "gray" , alpha = 0.35 ) +
        annotate("rect", xmin = log10(demand_bonferroni_threshold) , xmax =  0 , ymin = viper_threshold , ymax = -Inf , fill = "gray" , alpha = 0.35 ) +
        
        geom_point( data = demand_table_filtered %>% filter( bonferroni >= demand_bonferroni_threshold | viper_score >= viper_threshold ) , alpha = 0.25 , stroke = 1 , size = 2 , color = "black" , shape = 20 ) +
        
        geom_text_repel( data = demand_table_filtered %>% filter( bonferroni < demand_bonferroni_threshold & viper_score < viper_threshold ) , 
                         aes( label = moaGene ) , 
                         size = 6 , max.overlaps = 20 , force = TRUE , force_pull = FALSE , segment.color = "gray85" ,
                         seed = 666 , min.segment.length = unit(0, 'lines') , nudge_x = -3 , nudge_y = -1.5 ) +
        geom_point( data = demand_table_filtered %>% filter( bonferroni < demand_bonferroni_threshold & viper_score < viper_threshold ) , alpha = 1 , stroke = 1 , size = 2 , fill = NA , shape = 20 ) +        
        
        geom_vline( xintercept = log10(demand_bonferroni_threshold) , color = "red" , alpha = 0.75 , linetype = "longdash") + 
        geom_hline( yintercept = viper_threshold , color = "red" , alpha = 0.75 , linetype = "longdash") + 
        xlim(c(demand_min_limit-5,0)) +
        ylim(c(viper_min_limit-5,0)) +
        xlab( paste0( "DeMAND Bonferroni adjusted pvalue (log10)") ) +
        ylab( paste0( "VIPER Score") ) +
        ggtitle( paste0( "VIKING Plot for " , aDrug ) ) +
        theme_minimal() +
        theme( 
          plot.title = element_text(size=15,face="bold") , # family="Open Sans") ,
          axis.title.x = element_text(size=12,face="bold") , # family="Open Sans") ,
          axis.title.y = element_text(size=12,face="bold") ) # family="Open Sans") )
      
      print(demand_plot)
      
      pdf(file.path(reports.dir,paste0("viking-plot-for-",aDrug,"-vs-null-analog-on-rnaseq.pdf")) )
      print(demand_plot)
      dev.off()
    }
    
  }
  
  demand_table_filtered$drug <- aDrug
  viking_table <- rbind(viking_table,demand_table_filtered)  
  filename <- file.path(processed_data.dir,"viking-table-vs-null-analog-from-rnaseq.rds" )
  saveRDS(viking_table,filename)
}

print_msg_info("***** Analysis Completed *****")

