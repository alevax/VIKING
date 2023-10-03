
require(crayon)
create_workspace <- function( run_dir = "my-run" ,
                              isWorkspaceToClean = FALSE ,
                              experiments_dir = "experiments" )
{
    if (isWorkspaceToClean)
    {
    	vars_to_clean <- ls()
    	vars_to_clean <- vars_to_clean[ !(vars_to_clean %in% c("experiments_dir","run_dir")) ]
    	rm(list = vars_to_clean )
    }

    experiments.dir <- file.path( getwd() , experiments_dir ) 
    dir.create( experiments.dir , showWarnings = FALSE )
    # experiment.name <- paste( Sys.Date() , format( Sys.time(), "%H%M" ) , experiment.tag , sep = "-" )
    experiment.name <- paste( Sys.Date() , run_dir , sep = "-" )
    experiment.dir <- file.path( experiments.dir , experiment.name )
    reports.dir <<- file.path( experiment.dir , "reports" )
    processed_data.dir <<- file.path( experiment.dir , "processed_data" )
    dir.create( experiment.dir , showWarnings = FALSE )
    dir.create( reports.dir , showWarnings = FALSE )
    dir.create( processed_data.dir , showWarnings = FALSE )
    
    message( black$bgYellow$bold( ">>> Created directory: " , experiment.dir) )
    message( black$bgCyan( ">>> >> Created subdirectory: " , reports.dir) )
    message( black$bgCyan( ">>> >> Created subdirectory: " , processed_data.dir) )
}

	print_msg_warn <- function(...) {
		info_msg <- black $ bgYellow $ bold
		text <- list(...)
		# text <- str_c(text,sep = "")
		cat(info_msg(text),"\n",sep = "")
	}	
	print_msg_info <- function(...) {
		info_msg <- yellow $ bold
		text <- list(...)
		# text <- str_c(text,sep = "")
		cat(info_msg(text),"\n",sep = "")
	}		

  require(data.table)
  preppi_table <- fread( "data/preppi_final600.txt" , header = TRUE ) #, skipNul = TRUE , sep = "\t")
	generateDemandNetworkFromPrePPI <- function( likelihood.threshold = 0.9 , regulators = "TF" , structural_score.threshold = 0 , n_neighbors = 0 )
  {
    library(data.table)
	  if ( is.vector(regulators) && length(regulators) > 1 ){ 
	    .tfs <- regulators
	  } else {
	    .tfs <- getRegulators( kind_of_id = "symbols" , kind_of_regulators = regulators )  
	  }
    
    .dt <- preppi_table[ symbol_prot1 %in% .tfs | symbol_prot2 %in% .tfs
      ][ likelihood >= likelihood.threshold
				][ str_score >= structural_score.threshold
				  ][ !is.na(symbol_prot1) & !is.na(symbol_prot2)
				    ][ symbol_prot1 != "" & symbol_prot2 != ""
          	  ][ , list("Gene1"=symbol_prot1,"Gene2"=symbol_prot2) ]
    nrow(.dt)
    .dt <- .dt %>% setkey(Gene1,Gene2)
    .dt <- unique(.dt)
    nrow(.dt)
    
    if (is.numeric(n_neighbors) && n_neighbors > 0 )
    {
      .dt_inner_inner <- data.table()
      .other_prots_to_expand <- unique(c(.dt$Gene1,.dt$Gene2))
      for ( i in 1:n_neighbors )
      {
        .dt_inner <- preppi_table[ symbol_prot1 %in% .other_prots_to_expand | symbol_prot2 %in% .other_prots_to_expand
          ][ likelihood >= likelihood.threshold
            ][ str_score >= structural_score.threshold
              ][ !is.na(symbol_prot1) & !is.na(symbol_prot2)
                ][ symbol_prot1 != "" & symbol_prot2 != ""
                  ][ , list("Gene1"=symbol_prot1,"Gene2"=symbol_prot2) ]
        nrow(.dt_inner)
        .dt_inner <- .dt_inner %>% setkey(Gene1,Gene2)
        .dt_inner <- unique(.dt_inner)
        nrow(.dt_inner) 
        .other_prots_to_expand <- unique(c(.dt_inner$Gene1,.dt_inner$Gene2))
      }
      
      .dt <- do.call(rbind,list(.dt,.dt_inner))
      .dt <- .dt %>% setkey(Gene1,Gene2)
      .dt <- unique(.dt)
    }
    
    .df <- as.data.frame(.dt)
    .df$ppi <- 1
    .demand.network <- as.matrix(.df)
		print( message("- Network with " , nrow(.demand.network) , " interactions" ) )
    return(.demand.network)
  }

	expandPrePPIClosestNeighbor <- function( symbol , likelihood.threshold = 0.9 , structural_score.threshold = 0 , with.final.score=FALSE)
  {
    require(data.table)

	  if (with.final.score){
	    .dt_1 <- preppi_table[ symbol_prot1 %in% symbol
	      ][ likelihood >= likelihood.threshold
	        ][ str_score >= structural_score.threshold
	          ][ !is.na(symbol_prot1) & !is.na(symbol_prot2)
	            ][ , list("query_symbol"=symbol_prot1,"closest_neighbor"=symbol_prot2,"score"=final_score,"likelihood"=likelihood,"str_score"=str_score) ]
	    .dt_2 <- preppi_table[ symbol_prot2 %in% symbol
	      ][ likelihood >= likelihood.threshold
	        ][ str_score >= structural_score.threshold
	          ][ !is.na(symbol_prot1) & !is.na(symbol_prot2)
	            ][ , list("query_symbol"=symbol_prot2,"closest_neighbor"=symbol_prot1,"score"=final_score,"likelihood"=likelihood,"str_score"=str_score) ]	    
	    .dt <- rbind(.dt_1,.dt_2)
	    
	  } else {
	    .dt_1 <- preppi_table[ symbol_prot1 %in% symbol
	      ][ likelihood >= likelihood.threshold
	        ][ str_score >= structural_score.threshold
	          ][ !is.na(symbol_prot1) & !is.na(symbol_prot2)
	            ][ , list("query_symbol"=symbol_prot1,"closest_neighbor"=symbol_prot2) ]
	    .dt_2 <- preppi_table[ symbol_prot2 %in% symbol
	      ][ likelihood >= likelihood.threshold
	        ][ str_score >= structural_score.threshold
	          ][ !is.na(symbol_prot1) & !is.na(symbol_prot2)
	            ][ , list("query_symbol"=symbol_prot2,"closest_neighbor"=symbol_prot1) ]	    
	    .dt <- rbind(.dt_1,.dt_2)
	  }

    .dt <- unique(.dt)
    .df <- as.data.frame(.dt)
		print( message("- Network with " , nrow(.df) , " interactions" ) )
    return(.df)
  }	
	
	expandPrePPIClosestNeighbor_multiple <- function( symbol.list , likelihood.threshold = 0.9 , structural_score.threshold = 100 , with.final.score=FALSE ) {
	  res <- as_tibble( do.call( rbind , lapply( unique(symbol.list) , expandPrePPIClosestNeighbor , likelihood.threshold = likelihood.threshold , structural_score.threshold = structural_score.threshold , with.final.score=with.final.score ) ) )
	  return(res)
	}	
	
	doStouffer <- function( aMatrix , weights = NULL ) {
		if ( is.null(weights) )
		{
			res <- rowSums(aMatrix) / sqrt( ncol(aMatrix) )
		}
		else {
			
			res <- (aMatrix %*% weights) / sqrt( sum(weights**2) )
		}
		return(res)
	}

	 getRegulators <- function( kind_of_id = "entrez" , kind_of_regulators = "all" , filename = "data/regulators.csv" )
  {
  	regulators.df <- read.csv2( file = filename , stringsAsFactors =  FALSE )

  	.id <- match.arg(arg = kind_of_id,c("entrez","symbols"))
  	if ( .id == "entrez" ) .id <- "entrez.id" else .id <- "gene.symbols"
  	.regulators <- match.arg(arg = kind_of_regulators,c("all","TF","coTF","SIG"))

  	require(dplyr)
  	if(.regulators=="all")
  	{
			ret <- regulators.df %>% dplyr::select(.id)
  	} else {
  		ret <- regulators.df %>% dplyr::select(.id) %>% dplyr::filter( regulators.df$type.of.regulators == .regulators )
  	}
  	return( (ret[,1]) )
  }
	