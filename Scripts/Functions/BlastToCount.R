BlastToCount <- function(name, sep="\t", colnames, filter, tax.sep, tax.colnames, tax.filter, summarise.by){
  
  #' @title BLASTn file to count data converter
  #'
  #' @description This function allows to filter raw BLASTn files and converting them to count data. The function returns a list containing first, the count data, second, the number of taxa excluded at each filter step, and third, the mean BLASTn values of the table. 
  #' 
  #' @section Required packages:
  #' 'dplyr'
  #' 'tidyr' 
  #'
  #' @param name (\emph{character}). A complete filename \strong{including the path (!)} to a raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}). 
  #' @param sep (\emph{character}). A field separator character that separates the columns of the raw BLASTn file. By default, sep = '\t'.
  #' @param colnames (\emph{character vector}). A vector with names to specify the content of the columns of the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}). Possible names are: 'qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', and 'qcovhsp'. 
  #' @param filter (\emph{named nested list}). A nested list that is used to filter the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}). Each column of the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}) can be used to exclude data. To this end, please specify per filter (i.e. per column) a list that is named after one column of the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}) and contains 2 elements: 1.) one of the characters 'e', 'ne', 'l', 'g', 'le' or 'ge' (e=equal, ne=not equal, l=less than, g=greater than, le=less than or equal to, ge=greater than or equal to), 2.) a variable that defines a threshold. All 'filter lists' are stored in another list and passed to the argument 'filter'. Example to exclude pident<80 and qlen<150: filter=list(pident=list('l',80), qlen=list('l',150)).
  #' @param tax.sep (\emph{character}). A field separator character that separates the column 'sacc' containing the taxonomic information.
  #' @param tax.colnames (\emph{character vector}). A vector with names to specify the contents of the taxonomy columns of the separated 'sacc' column of the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}).
  #' @param tax.filter (\emph{named nested list}). A nested list that is used to filter the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}) by the taxonomy columns. Each taxonomy column of the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}) can be used to exclude data. To this end, please specify per taxonomy filter (i.e. per taxonomy column) a list that is named after one taxonomy column of the raw BLASTn file (\href{https://www.metagenomics.wiki/tools/blast/blastn-output-format-6}{output format 6}) and contains 2 elements: 1.) one of the characters 'e' or 'ne' (e=equal, ne=not equal), 2.) a variable that defines a taxonomic name to filter by. All 'taxonomic filter lists' are stored in another list and passed to the argument 'tax.filter'. Example to exclude Bacteria: tax.filter=list(Kingdom=list('e','Bacteria')).
  #' @param summarise.by (\emph{character vector}). A vector with names to specify the columns across which the BLASTn table is summarised and the frequencies are counted.
  #'
  #' @return
  #' @export
  #'
  #' @examples
  #' 
  #############################################################################################
  #################################### Check prerequisites ####################################
  #------------------------- Check if needed packages are installed --------------------------#
  
  if(!"dplyr" %in% rownames(installed.packages())){
    stop("Please install the package 'dplyr'.")
  }
  
  if(!"tidyr" %in% rownames(installed.packages())){
    stop("Please install the package 'tidyr'.")
  }
  
  #--------------------------- Check if needed arguments are valid ---------------------------#
  
  # Check if argument'name' is provided and valid
  if(missing(name)){
    stop("The argument 'name' is missing.\nPlease provide a complete file name (including the path!) to a raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).")
  } else if(class(name) != "character"){
    stop("The argument 'name' is not a character.\nPlease provide a complete file name (including the path!) to a raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).")
  } else if(length(name) != 1){
    stop("The argument 'name' has more tha one element.\nPlease provide a complete file name (including the path!) to a raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).")
  }
  
  # Check if argument 'sep' is valid
  if(class(sep) != "character"){
    stop("The argument 'sep' is not a character.\nPlease specify a field separator character that separates the columns (e.g. 'qseqid', 'evalue', 'pident', etc.) of the raw BLASTn file. By default, sep = '\\t'.")
  } else if(length(sep) != 1){
    stop("The argument 'sep' has more tha one element.\nPlease specify a field separator character that separates the columns (e.g. 'qseqid', 'evalue', 'pident', etc.) of the raw BLASTn file. By default, sep = '\\t'.")
  }
  
  # Check if argument 'colnames' is valid
  if(!missing(colnames)){
    if(!(is.vector(colnames) & is.atomic(colnames))){
      stop("The argument 'colnames' is not a vector.\nPlease provide a vector with names to specify the content of the columns of the raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).\nPossible names are: 'qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', and 'qcovhsp'.\nBy default, colnames = c('qseqid', 'qlen', 'sacc', 'bitscore', 'evalue','length', 'nident', 'pident').")
    } else if(class(colnames) != "character"){
      stop("The argument 'colnames' is not a character vector.\nPlease provide a vector with names to specify the content of the columns of the raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).\nPossible names are: 'qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', and 'qcovhsp'.\nBy default, colnames = c('qseqid', 'qlen', 'sacc', 'bitscore', 'evalue','length', 'nident', 'pident').")
    } else if(!all(colnames %in% c("qseqid", "qgi", "qacc", "qaccver", "qlen", "sseqid", "sallseqid", "sgi", "sallgi", "sacc", "saccver", "sallacc", "slen", "qstart", "qend", "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", "positive", "gapopen", "gaps", "ppos", "frames", "qframe", "sframe", "btop", "staxids", "sscinames", "scomnames", "sblastnames", "sskingdoms", "stitle", "salltitles", "sstrand", "qcovs", "qcovhsp"))){
      stop("Not all names defined in the argument 'colnames' are valid.\nPlease provide a vector with names to specify the content of the columns of the raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).\nPossible names are: 'qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', and 'qcovhsp'.\nBy default, colnames = c('qseqid', 'qlen', 'sacc', 'bitscore', 'evalue','length', 'nident', 'pident').")
    }
    
    #----------------- Load data -----------------#
    # Load file, rename columns
    BlastFile <- setNames(data.frame(readLines(name)), "V1")
    suppressWarnings(BlastFile <- data.frame(tidyr::separate(BlastFile, V1, paste0("V", seq(200)), sep = sep)))
    BlastFile <- BlastFile[,colSums(is.na(BlastFile)) != nrow(BlastFile)]
    
    # Assign column names
    if(ncol(BlastFile) == length(colnames)){
      colnames(BlastFile) <- colnames
    } else {
      stop(paste("The number of columns after the separation of the blast file by the defined field separator is not equal to the number of specified new column names. Please specify", ncol(BlastFile) ,"column names in the argument 'colnames'." ))
    }
  } else {
    stop("The argument 'colnames' is missing.\nPlease provide a vector with names to specify the content of the columns of the raw BLASTn file (output format 6, see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).\nPossible names are: 'qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', and 'qcovhsp'.\nBy default, colnames = c('qseqid', 'qlen', 'sacc', 'bitscore', 'evalue','length', 'nident', 'pident').")
  }
  
  # Check if argument 'filter' is valid
  if(!missing(filter)){
    if(class(filter) != "list"){
      stop("The argument 'filter' is not a list.\nPlease provide a nested list that is used to filter the raw BLASTn file.\nEach column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per filter (i.e. per column) a list that is named after one column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e', 'ne', 'l', 'g', 'le' or 'ge' (e=equal, ne=not equal, l=less than, g=greater than, le=less than or equal to, ge=greater than or equal to),\n2.) a variable that defines a threshold.\nAll 'filter lists' are stored in another list and passed to the argument 'filter'.\nExample to exclude pident<80 and qlen<150: filter=list(pident=list('l',80), qlen=list('l',150)).")
    }else if(!all(names(filter) %in% colnames)){
      stop("The names of the nested lists in the 'filter list' cannot all be found in the column names of the BLASTn table.\nPlease provide a nested list that is used to filter the raw BLASTn file.\nEach column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per filter (i.e. per column) a list that is named after one column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e', 'ne', 'l', 'g', 'le' or 'ge' (e=equal, ne=not equal, l=less than, g=greater than, le=less than or equal to, ge=greater than or equal to),\n2.) a variable that defines a threshold.\nAll 'filter lists' are stored in another list and passed to the argument 'filter'.\nExample to exclude pident<80 and qlen<150: filter=list(pident=list('l',80), qlen=list('l',150)).")
    }else if(any(sapply(filter, class) != "list")){
      stop("The argument 'filter' is not a named nested list (a list that contains named lists).\nPlease provide a nested list that is used to filter the raw BLASTn file.\nEach column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per filter (i.e. per column) a list that is named after one column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e', 'ne', 'l', 'g', 'le' or 'ge' (e=equal, ne=not equal, l=less than, g=greater than, le=less than or equal to, ge=greater than or equal to),\n2.) a variable that defines a threshold.\nAll 'filter lists' are stored in another list and passed to the argument 'filter'.\nExample to exclude pident<80 and qlen<150: filter=list(pident=list('l',80), qlen=list('l',150)).")
    }else if(length(unlist(filter))/length(filter) != 2){
      stop("The nested lists within the 'filter llist' do not all contain 2 variables respectively.\nPlease provide a nested list that is used to filter the raw BLASTn file.\nEach column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per filter (i.e. per column) a list that is named after one column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e', 'ne', 'l', 'g', 'le' or 'ge' (e=equal, ne=not equal, l=less than, g=greater than, le=less than or equal to, ge=greater than or equal to),\n2.) a variable that defines a threshold.\nAll 'filter lists' are stored in another list and passed to the argument 'filter'.\nExample to exclude pident<80 and qlen<150: filter=list(pident=list('l',80), qlen=list('l',150)).")
    }
  } else {
    message("Data was not filtered by any column of the blast table.")
  }
  
  # Check if argument 'tax.sep' is valid
  if(!missing(tax.sep)){
    if(class(tax.sep) != "character"){
      stop("The argument 'tax.sep' is not a character.\nPlease specify a field separator character that separates the column 'sacc' containing the taxonomic information.")
    } else if(length(tax.sep) != 1){
      stop("The argument 'tax.sep' has more tha one element.\nPlease specify a field separator character that separates the column 'sacc' containing the taxonomic information.")
    }
    
    # Check if argument 'tax.colnames' is valid
    if(!missing(tax.colnames)){
      if(!(is.vector(tax.colnames) & is.atomic(tax.colnames))){
        stop("The argument 'tax.colnames' is not a character vector.\nPlease provide a vector with names to specify the contents of the taxonomy columns of the separated 'sacc' column of the raw BLASTn file.")
      } else if(class(tax.colnames) != "character"){
        stop("The argument 'tax.colnames' is not a character vector.\nPlease provide a vector with names to specify the contents of the taxonomy columns of the separated 'sacc' column of the raw BLASTn file.")
      }
    } else {
      stop("The argument 'tax.colnames' is missing.\nHowever, the argument 'tax.colnames' is required if the argument 'tax.sep' is defined.\nThus, please provide a vector with names to specify the contents of the taxonomy columns of the separated 'sacc' column of the raw BLASTn file.")
    }
    
    # Check if column 'sacc' is given
    if(!"sacc" %in% colnames){
      stop("The BLASTn table can not be separated or filtered by the taxonomy as the column 'sacc' (this column must contain all taxonomic information) is missing.")
    }
    
    # Check if argument 'tax.filter' is valid
    if(!missing(tax.filter)){
      if(class(tax.filter) != "list"){
        stop("The argument 'tax.filter' is not a list.\nPlease provide a nested list that is used to filter the raw BLASTn file by the taxonomy columns.\nEach taxonomy column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per taxonomy filter (i.e. per taxonomy column) a list that is named after one taxonomy column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e' or 'ne' (e=equal, ne=not equal),\n2.) a variable that defines a taxonomic name to filter by.\nAll 'taxonomic filter lists' are stored in another list and passed to the argument 'tax.filter'.\nExample to exclude Bacteria: tax.filter=list(Kingdom=list('e','Bacteria')).")
      }else if(!all(names(tax.filter) %in% tax.colnames)){
        stop("The names of the nested lists in the list 'tax.filter' cannot all be found in the column names of the taxonomy column.\nPlease provide a nested list that is used to filter the raw BLASTn file by the taxonomy columns.\nEach taxonomy column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per taxonomy filter (i.e. per taxonomy column) a list that is named after one taxonomy column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e' or 'ne' (e=equal, ne=not equal),\n2.) a variable that defines a taxonomic name to filter by.\nAll 'taxonomic filter lists' are stored in another list and passed to the argument 'tax.filter'.\nExample to exclude Bacteria: tax.filter=list(Kingdom=list('e','Bacteria')).")
      }else if(any(sapply(tax.filter, class) != "list")){
        stop("The argument 'tax.filter' is not a named nested list (a list that contains named lists).\nPlease provide a nested list that is used to filter the raw BLASTn file by the taxonomy columns.\nEach taxonomy column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per taxonomy filter (i.e. per taxonomy column) a list that is named after one taxonomy column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e' or 'ne' (e=equal, ne=not equal),\n2.) a variable that defines a taxonomic name to filter by.\nAll 'taxonomic filter lists' are stored in another list and passed to the argument 'tax.filter'.\nExample to exclude Bacteria: tax.filter=list(Kingdom=list('e','Bacteria')).")
      }else if(length(unlist(tax.filter))/length(tax.filter) != 2){
        stop("The named nested lists within the list 'tax.filter' do not all contain 2 variables respectively.\nPlease provide a nested list that is used to filter the raw BLASTn file by the taxonomy columns.\nEach taxonomy column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per taxonomy filter (i.e. per taxonomy column) a list that is named after one taxonomy column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e' or 'ne' (e=equal, ne=not equal),\n2.) a variable that defines a taxonomic name to filter by.\nAll 'taxonomic filter lists' are stored in another list and passed to the argument 'tax.filter'.\nExample to exclude Bacteria: tax.filter=list(Kingdom=list('e','Bacteria')).")
      }else if(any(sapply(unlist(tax.filter), class) != "character")){
        stop("The named nested lists within the list 'tax.filter' do not all contain 2 character (!) variables respectively.\nPlease provide a nested list that is used to filter the raw BLASTn file by the taxonomy columns.\nEach taxonomy column of the raw BLASTn file can be used to exclude data.\nTo this end, please specify per taxonomy filter (i.e. per taxonomy column) a list that is named after one taxonomy column of the raw BLASTn file and contains 2 elements:\n1.) one of the characters 'e' or 'ne' (e=equal, ne=not equal),\n2.) a variable that defines a taxonomic name to filter by.\nAll 'taxonomic filter lists' are stored in another list and passed to the argument 'tax.filter'.\nExample to exclude Bacteria: tax.filter=list(Kingdom=list('e','Bacteria')).")
      }
    } else {
      message("Data was not filtered by the taxonomy columns.")
    }
    
  } else {
    message("Data was not separated or filtered by the taxonomy columns.")
  }
  
  # Check if argument 'summarise.by' is valid
  if(!missing(summarise.by)){
    if(!(is.vector(summarise.by) & is.atomic(summarise.by))){
      stop("The argument 'summarise.by' is not a vector.\nPlease provide a vector with names to specify the columns across which the BLASTn table is summarised and the frequencies are counted.")
    } else if(class(summarise.by) != "character"){
      stop("The argument 'summarise.by' is not a character vector.\nPlease provide a vector with names to specify the columns across which the BLASTn table is summarised and the frequencies are counted.")
    }
  } else {
    message("Blast table was not summarised.")
  }
  
  #############################################################################################
  ####################### Separate, filter and summaries the blast file #######################
  
  # Convert class, format numeric columns
  for (x in colnames(BlastFile)) {
    if(x %in% c("qgi", "qlen", "sgi", "sallgi", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", "positive", "gapopen", "gaps", "ppos", "staxids", "qcovs", "qcovhsp")){
      BlastFile[,x] <- as.numeric(BlastFile[,x])
    }
  }
  
  # Get total number of assigned taxa, save them in list
  NrTaxa <- list(Total=nrow(BlastFile))
  
  # Filter -------------------------------------------------------------------------------------
  # Check if 'filter' are provided and valid
  if(!missing(filter)){
    
    # Loop over the filter
    for (i in seq(filter)) {
      
      # Check if filter is valid
      if(!filter[[i]][[1]] %in% c("e", "ne", "l", "g", "le", "ge")){
        stop(paste0("Filter '", filter[[i]][[1]], "' not valid, please choose between 'e', 'ne', 'l', 'g', 'le', or 'ge' (e=equals, ne=not equal to, l=less than, g=greater than, le=less than or equal to, ge=greater than or equal to)."))
      }
      
      # Filter
      if(filter[[i]][[1]] == "e"){
        # Equals
        BlastFile <- BlastFile[!(BlastFile[,names(filter)[i]] == filter[[i]][[2]]),]
        
        # Get total number of assigned taxa
        NrTaxa[[i+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
        names(NrTaxa)[i+1] <- paste(names(filter)[i],"==", filter[[i]][[2]])
      } else if(filter[[i]][[1]] == "ne"){
        # Not equal to
        BlastFile <- BlastFile[!(BlastFile[,names(filter)[i]] != filter[[i]][[2]]),]
        
        # Get total number of assigned taxa
        NrTaxa[[i+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
        names(NrTaxa)[i+1] <- paste(names(filter)[i],"!=", filter[[i]][[2]])
      } else if(filter[[i]][[1]] == "l"){
        # Less than
        if(class(filter[[i]][[2]])!="numeric"){
          stop(paste("Filter 'less than' can not be used for the as character (!) defined filter", names(filter)[i], "<", filter[[i]][[2]]))
        }
        
        BlastFile <- BlastFile[!(BlastFile[,names(filter)[i]] < filter[[i]][[2]]),]
        
        # Get total number of assigned taxa
        NrTaxa[[i+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
        names(NrTaxa)[i+1] <- paste(names(filter)[i],"<", filter[[i]][[2]])
      } else if(filter[[i]][[1]] == "g"){
        # Greater than
        if(class(filter[[i]][[2]])!="numeric"){
          stop(paste("Filter 'less than' can not be used for the as character (!) defined filter", names(filter)[i], "<", filter[[i]][[2]]))
        }
        
        BlastFile <- BlastFile[!(BlastFile[,names(filter)[i]] > filter[[i]][[2]]),]
        
        # Get total number of assigned taxa
        NrTaxa[[i+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
        names(NrTaxa)[i+1] <- paste(names(filter)[i],">", filter[[i]][[2]])
      } else if(filter[[i]][[1]] == "le"){
        # Less than or equal to
        if(class(filter[[i]][[2]])!="numeric"){
          stop(paste("Filter 'less than' can not be used for the as character (!) defined filter", names(filter)[i], "<", filter[[i]][[2]]))
        }
        
        BlastFile <- BlastFile[!(BlastFile[,names(filter)[i]] <= filter[[i]][[2]]),]
        
        # Get total number of assigned taxa
        NrTaxa[[i+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
        names(NrTaxa)[i+1] <- paste(names(filter)[i],"<=", filter[[i]][[2]])
      } else {
        # Greater than or equal to
        if(class(filter[[i]][[2]])!="numeric"){
          stop(paste("Filter 'less than' can not be used for the as character (!) defined filter", names(filter)[i], "<", filter[[i]][[2]]))
        }
        
        BlastFile <- BlastFile[!(BlastFile[,names(filter)[i]] >= filter[[i]][[2]]),]
        
        # Get total number of assigned taxa
        NrTaxa[[i+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
        names(NrTaxa)[i+1] <- paste(names(filter)[i],">=", filter[[i]][[2]])
      }
    }
  }
  
  # Separate taxonomy --------------------------------------------------------------------------
  if(!missing(tax.sep)){
    
    # Split taxonomy column 'SubjectAccession' 
    suppressWarnings(BlastFile <- data.frame(tidyr::separate(BlastFile, sacc, paste0("V", seq(200)), sep = tax.sep)))
    BlastFile <- BlastFile[,colSums(is.na(BlastFile)) != nrow(BlastFile)]
    
    # Assign new column names
    if(ncol(BlastFile)-(length(colnames)-1) > length(tax.colnames)){
      colnames(BlastFile)[which(colnames(BlastFile) == "V1"):(length(tax.colnames)-1+which(colnames(BlastFile) == "V1"))] <- tax.colnames 
    } else {
      colnames(BlastFile)[which(colnames(BlastFile) == "V1"):(length(tax.colnames)-1+which(colnames(BlastFile) == "V1"))] <- 
        tax.colnames[seq(ncol(BlastFile)-(length(colnames)-1))]
    }

    # Filter -------------------------------------------------------------------------------------
    # Check if 'filter' are provided and valid
    if(!missing(tax.filter)){
      
      # Loop over the tax.filter
      for (i in seq(tax.filter)) {
   
        # Check if filter is valid
        if(!tax.filter[[i]][[1]] %in% c("e", "ne")){
          stop(paste0("Filter '", tax.filter[[i]][[1]], "' not valid, please choose between 'e', 'ne' (e=equals or ne=not equal to)."))
        }
        
        # Change class of the filter column
        BlastFile[,names(tax.filter)[i]] <- as.character(BlastFile[,names(tax.filter)[i]])
        
        # Filter
        if(tax.filter[[i]][[1]] == "e"){
          # Equals
          
          BlastFile <- BlastFile[!(BlastFile[,names(tax.filter)[i]] == tax.filter[[i]][[2]]),]
          
          # Get total number of assigned taxa
          NrTaxa[[length(NrTaxa)+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
          names(NrTaxa)[length(NrTaxa)] <- paste(names(tax.filter)[i],"==", tax.filter[[i]][[2]])
        } else {
          # Not equal to
          
          BlastFile <- BlastFile[!(BlastFile[,names(tax.filter)[i]] != tax.filter[[i]][[2]]),]
          
          # Get total number of assigned taxa
          NrTaxa[[length(NrTaxa)+1]] <- (NrTaxa$Total - (sum(unlist(NrTaxa)) - NrTaxa$Total))-nrow(BlastFile)
          names(NrTaxa)[length(NrTaxa)] <- paste(names(tax.filter)[i],"!=", tax.filter[[i]][[2]])
        } 
      }
    }
  }
  
  # Get means of the numeric columns  ----------------------------------------------------------
  stat <- BlastFile[,sapply(BlastFile,class) %in% c("integer","numeric")]

  # Unlist NrTaxa
  NrTaxa <- unlist(NrTaxa)
  
  # Create count table -------------------------------------------------------------------------
  if(!missing(summarise.by)){
    if(all(summarise.by %in% colnames(BlastFile))){
      
      # Summarise
      CountData <- data.frame(dplyr::count(dplyr::group_by_all(BlastFile[,summarise.by])))

    } else {
      warning("not all column names specified in th argument 'summarise.by' can be found in the blast table, thus the table is not summarised.")
      CountData <- BlastFile
    }
  } else {
    CountData <- BlastFile
  }
  
  # Export results -----------------------------------------------------------------------------
  return(list(CountData, NrTaxa, stat))
}
