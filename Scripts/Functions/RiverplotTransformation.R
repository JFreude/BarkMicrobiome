riverplot.transformation<-function(counts, taxonomy, key="Order", nr.taxa=10, taxa.abundance="mean", total,
                                   order.taxonomic.ranks=T, color.palette="viridis", color.others="darkgray", 
                                   not.shown, labels="names.rel.abundance", decimal.places=1, show.others=T){
  #' @title Sankey diagram
  #'
  #' @description 
  #' 
  #' The function creates a Sankey diagram showing the relative proportion of the most abundant taxa across their taxonomic ranks.
  #' 
  #' One can choose \strong{how many taxa} are displayed (argument \strong{\emph{nr.taxa}}) and also the lowest taxonomic rank (argument \strong{\emph{key}}). For example, if we show the 10 most abundant orders of a community, we can see to which higher taxonomic ranks (clade, phyla, etc.) they belong, and also, based on the thickness of the branches, what proportion of the total community they represent. Taxa that are not among the most abundant are grouped into 'Others'. 
  #' 
  #' To \strong{select the most abundant taxa}, one can either use their mean abundance or their total (summed) abundance across all samples (argument \strong{\emph{taxa.abundance}}). For the calculation of the \strong{relative abundance of taxa}, one can choose to calculate the mean relative abundance across all samples (mean(abundance of one taxon)/mean(abundance of all taxa)) or the summed relative abundance across all samples (sum(abundance of one taxon)/sum(abundance of all taxa)). In addition, one can define the total value (argument \strong{\emph{total}}) for the calculation of the relative abundance (abundance of taxon/total). For example, if we want to calculate 2 Sankey diagrams, one for the prokaryotes and one for the eukaryotes, but still show the relative mean abundance of the taxa in the entire community (prokaryotes + eukaryotes), we can pass the mean total abundance of prokaryotes + eukaryotes to the argument \strong{\emph{total}}.  
  #' 
  #' The \strong{colors} that are used for the Sankey diagram can either be selected from a 'Viridis Color Palette' or customized. If you want to customize the colors (argument \strong{\emph{color.palette}} = 'customize'), you are asked for colors for the highest taxonomic rank (e.g. Phylum) and the lowest taxonomic rank (defined by the argument \strong{\emph{key}}). The color for the taxa that are not among the most abundant can be specified with the argument \strong{\emph{color.others}}. Color gradients are created between the colors of the highest and lowest taxonomic ranks.
  #' 
  #' The most abundant taxa are shown in \strong{alphabetical order}. For this, the taxa are sorted alphabetically from the lowest taxonomic rank (e.g. Species) to the highest taxonomic rank (e.g. Phylum). The argument \strong{\emph{not.shown}} gives the option to sort by a taxonomic rank or group, but not to show it in the diagram. For example, if we want to create a Sankey diagram for the eukaryotes where the taxa of fungi, protists, and metazoa are grouped together, we can add a column with this assignment to the taxonomy and pass the name of this column to the argument \strong{\emph{not.shown}}, thus not showing this column itself in the diagram.   
  #' 
  #' At last, we can define the labels for the Sankey diagram with the argument \strong{\emph{labels}}. We can choose to display no labels, the taxa names, the total abundance, the relative abundance, or a combination of names and total/relative abundance. If the mean total/relative abundance is calculated and displayed, also the standard deviation is shown. The number of decimal places for the total/relative abundance can be defined with the argument \strong{\emph{decimal.places}}. 
  #'
  #' @section Required packages:
  #' 'dplyr', riverplot, viridis 
  #' 
  #' @param counts (\emph{data.frame/matrix}). A data.frame/matrix with count data (numerics). Each row contains all information to one taxon for each sample. Column names must be sample-IDs and row names taxa-IDs. Taxa-IDs must match the taxa-IDs of the taxonomy data.
  #' @param taxonomy (\emph{data.frame/matrix}). A data.frame/matrix with taxonomic ranks, formatted as characters/factors. Each column is one taxonomic rank, providing the classification of each taxon. Each row contains all taxonomic information for one taxon. Column names must be taxonomic ranks, and row names taxa-IDs. Taxa-IDs must match the taxa-IDs of the taxonomy data.
  #' @param key (\emph{character}). A character that specifies the one column name of the \strong{taxonomy}. The defined taxonomic rank is the lowest taxonomic rank shown in the Sankey diagram.
  #' @param nr.taxa (\emph{integer}). An integer specifying the number of most abundant taxa shown in the Sankey diagram. Taxa that are not among the most abundant are binned into 'Others'. 
  #' @param taxa.abundance (\emph{character}). A character that specifies how the relative abundance are calculated, one of "mean", "sum". To select the most abundant taxa, one can either use their mean abundance ("mean") or their total abundance ("sum") across all samples. For the calculation of the relative abundance of taxa, either the mean relative abundance across all samples ("mean") are calculated or the summed relative abundance across all samples ("sum"). 
  #' @param total (\emph{numeric}). An optional numeric that define the total value for the calculation of the relative abundance (abundance of taxon/total).
  #' @param order.taxonomic.ranks (\emph{logical}). A logical, if TRUE the columns of the taxonomy table are sorted according to the number of unique entries, columns with more entries than the key column are excluded. If FALSE, the columns are not sorted according to the number of unique entries but the order is retained, all columns after the key column are excluded.
  #' @param color.palette (\emph{character}). A character that specifies the colors (colors from Viridis Color Palette or customized colors) that are used for the Sankey diagram, one of "magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo" or "customize". If you want to customize the colors (color.palette = "customize"), you are asked for colors for the highest taxonomic rank (e.g. Phylum) and the lowest taxonomic rank.
  #' @param color.others (\emph{character}). A character that specifies the color for the taxa that are not among the most abundant taxa.
  #' @param not.shown (\emph{character}). An optional character vector that contains one ore more column name of the \strong{taxonomy}. The most abundant taxa are sorted alphabetically from the lowest taxonomic rank (e.g. Species) to the highest taxonomic rank (e.g. Phylum). This argument defines column(s) of the taxonomy that is/are used to sort the taxa alphabetically, but is/are not shown in the Sankey diagram.
  #' @param labels (\emph{character}}). A character that specifies the labels for the Sankey diagram, one of "none", "names", "abundance", "rel.abundance", "names.abundance" or "names.rel.abundance". We can choose to display no labels ("none"), the taxa names ("names"), the total abundance ("abundance"), the relative abundance ("rel.abundance"), or a combination of names and total/relative abundance ("names.abundance", "names.rel.abundance"). If the mean total/relative abundance is calculated and displayed, also the standard standard deviation is shown. 
  #' @param decimal.places (\emph{integer}). An integer specifying the number of decimal places used for the labels for the relative abundance. 
  #' @param show.others (\emph{logical}). A logical, if TRUE taxa that are not among the most abundant are shown as 'Others' in the Sankey diagram.

  #############################################################################################
  #################################### Check prerequisites ####################################
  #-------------------------- Check if needed packages are installed -------------------------#
  if(!"dplyr" %in% rownames(installed.packages())){
    stop("Please install the package 'dplyr'.")
  }
  
  if(!"riverplot" %in% rownames(installed.packages())){
    stop("Please install the package 'riverplot'.")
  }
  
  if(!"viridis" %in% rownames(installed.packages())){
    stop("Please install the package 'viridis'.")
  }
  
  #--------------------------- Check if needed arguments are valid ---------------------------#
  # Check if 'counts' are provided
  if(missing(counts)){
    stop("argument counts is missing.\nPlease provide a data.frame/matrix with count data (numerics).\nColumn names must be sample-IDs and row names taxa-IDs.\nTaxa-IDs must match the taxa-IDs of the taxonomy data.")
  } else if(!class(counts) %in% c("data.frame","matrix")){
    stop("argument counts is not a data.frame or matrix.\nPlease provide a data.frame/matrix with count data (numerics).\nColumn names must be sample-IDs and row names taxa-IDs.\nTaxa-IDs must match the taxa-IDs of the taxonomy data.")
  } else if(!all(sapply(counts, is.numeric))){
    stop("argument counts, columns are not numeric.\nPlease provide a data.frame/matrix with count data (numerics).\nColumn names must be sample-IDs and row names taxa-IDs.\nTaxa-IDs must match the taxa-IDs of the taxonomy data.")
  }
  counts <- data.frame(counts)
  
  # Check if 'taxonomy' are provided
  if(missing(taxonomy)){
    stop("argument taxonomy is missing.\nPlease provide a data.frame/matrix with taxonomic ranks, formatted as characters/factors.\nColumn names must be taxonomic ranks, and row names taxa-IDs.\nTaxa-IDs must match the taxa-IDs of the taxonomy data.")
  } else if(!class(taxonomy) %in% c("data.frame","matrix")){
    stop("argument taxonomy is not a data.frame or matrix.\nPlease provide a data.frame/matrix with taxonomic ranks, formatted as characters/factors.\nColumn names must be taxonomic ranks, and row names taxa-IDs.\nTaxa-IDs must match the taxa-IDs of the taxonomy data.")
  }
  taxonomy <- data.frame(dplyr::mutate_if(taxonomy, is.factor, as.character))
  
  # Check if all row names of taxonomy and count data match
  if(!(all(rownames(taxonomy) %in% rownames(counts)) & all(rownames(counts) %in% rownames(taxonomy)))){
    if(!all(rownames(taxonomy) %in% rownames(counts))){
      stop("Not all taxa-IDs of the taxonomy data can be found in the count data.")
    } else {
      stop("Not all taxa-IDs of the count data can be found in the taxonomy data.")
    }
  }
  taxonomy <- taxonomy[match(rownames(counts), rownames(taxonomy)),]

  # Check if the value of the argument 'key' is valid
  if(class(key) != "character"){
    stop("The argument 'key' is not a character. Please provide a character that specifies one column name of the taxonomy table.\nThe defined taxonomic rank is the lowest taxonomic rank shown in the Sankey diagram.")
  } else if(!(key %in% colnames(taxonomy))){
    stop("The argument 'key' is not valid. The specified name is not found in the column names of the taxonomy table.\nPlease provide a character that specifies one column name of the taxonomy table.\nThe defined taxonomic rank is the lowest taxonomic rank shown in the Sankey diagram.")
  }
  
  # Check if the value of the argument 'nr.taxa' is valid
  if(!is.numeric(nr.taxa) | length(nr.taxa) != 1){
    stop("The argument 'nr.taxa' is not valid. Please provid an integer specifying the number of most abundant taxa (in terms of reads or OTUs) shown in the Sankey diagram.")
  } else if(nr.taxa%%1 != 0){
    stop("The argument 'nr.taxa' is not an integer. Please provid an integer specifying the number of most abundant taxa (in terms of reads or OTUs) shown in the Sankey diagram.")
  }
  
  # Check if the value of the argument 'taxa.abundance' is valid
  if(class(taxa.abundance) != "character"){
    stop("The argument 'taxa.abundance' is not a character. Please provide a character that specifies how the relative abundance are calculated, one of 'mean' or 'sum'.")
  } else if(!taxa.abundance %in% c("mean", "sum")){
    stop("The argument 'taxa.abundance' is not valid. Please provide a character that specifies how the relative abundance are calculated, one of 'mean' or 'sum'.")
  }

  # Check if the value of the argument 'color.palette' is valid
  if(class(color.palette) != "character"){
    stop("The argument 'color.palette' is not a character. Please provide a character that specifies the colors (colors from Viridis Color Palette or customized colors) that are used for the Sankey diagram, one of 'magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo' or 'customize'.")
  } else if(!color.palette %in% c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo", "customize")){
    stop("The argument 'color.palette' is not valid. . Please provide a character that specifies the colors (colors from Viridis Color Palette or customized colors) that are used for the Sankey diagram, one of 'magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo' or 'customize'.")
  }

  # Check if the value of the argument 'color.others' is valid
  if(class(color.others) != "character"){
    stop("The argument 'color.others' is not a character. Please provide a character that specifies the color for the taxa that are not among the most abundant taxa.")
  } else if(!tryCatch(is.matrix(col2rgb(color.others)), error = function(e) FALSE)){
    stop("The argument 'color.others' is not a color. Please provide a character that specifies the color for the taxa that are not among the most abundant taxa.")
  }

  # Check if the value of the argument 'labels' is valid
  if(class(labels) != "character"){
    stop("The argument 'labels' is not a character. Please provide a character that specifies the labels for the Sankey diagram, one of 'none', 'names', 'abundance', 'rel.abundance', 'names.abundance' or 'names.rel.abundance'.")
  } else if(!labels %in% c("none", "names", "abundance", "rel.abundance", "names.abundance", "names.rel.abundance")){
    stop("The argument 'labels' is not valid. Please provide a character that specifies the labels for the Sankey diagram, one of 'none', 'names', 'abundance', 'rel.abundance', 'names.abundance' or 'names.rel.abundance'.")
  }
  
  # Check if the value of the argument 'decimal.places' is valid
  if(!is.numeric(decimal.places) | length(decimal.places) != 1){
    stop("The argument 'decimal.places' is not valid. Please provid an integer specifying the number of decimal places used for the labels for the relative abundance.")
  } else if(decimal.places%%1 != 0){
    stop("The argument 'decimal.places' is not an integer. Please provid an integer specifying the number of decimal places used for the labels for the relative abundance.")
  }

  # Check if the value of the argument 'show.others' is valid
  if(class(show.others) != "logical"){
    stop("The argument 'show.others' is not valid. Please provid a logical (TRUE or FALSE) specifying if taxa that are not among the most abundant are shown as 'Others' in the Sankey diagram.")
  } 


  #-------------------------- Check if optional arguments are valid --------------------------#
  # Check if the value of the argument 'total' is valid
  if(!missing(total)){
    if(!is.numeric(total) | length(total) != 1){
      stop("The argument 'total' is not valid. Optionally, one can provid a numeric that define the total value for the calculation of the relative abundance (abundance of taxon/total).")
    } 
  }

  # Check if the value of the argument 'not.shown' is valid
  if(!missing(not.shown)){
    if(class(not.shown) != "character"){
      stop("The argument 'not.shown' is not a character vector. Optionally, one can provid a character vector that contains one ore more column names of the taxonomy that are used to sort the taxa alphabetically, but are not shown in the Sankey diagram.")
    } else if(!(all(not.shown %in% colnames(taxonomy)))){
      stop("The argument 'not.shown' is not valid as not all names are found in the taxonomy table. Optionally, one can provid a character vector that contains one ore more column names of the taxonomy that are used to sort the taxa alphabetically, but are not shown in the Sankey diagram.")
    }
  }

  #############################################################################################
  ###################################### Calculate Sankey #####################################
  
  #---------------------------- Summarise counts to most abundant ----------------------------#
  # Exclude all columns with higher taxonomic ranks than the new key column, exclude numeric columns
  taxonomy <- taxonomy[,sapply(taxonomy, is.character)]
  col.names <- sort(unlist(lapply(taxonomy,function(x){length(unique(x))})))
  if(order.taxonomic.ranks){
    taxonomy <- taxonomy[,names(col.names)[col.names <= length(unique(taxonomy[,key]))]]
  } else {
    taxonomy <- taxonomy[,seq(which(colnames(taxonomy) == key))]
  }
  
  # Summaries counts across all taxonomic level
  counts <- cbind(taxonomy, counts) 
  counts <- data.frame(dplyr::summarise_if(dplyr::group_by_if(counts, is.character), is.numeric, sum))
  
  # Calculate total reads or summed reads and sort in decreasing order
  if(taxa.abundance == "mean"){
    counts <- cbind(counts, TOTAL_COUNTS=apply(counts[,sapply(counts, is.numeric)],1, mean))
  } else {
    counts <- cbind(counts, TOTAL_COUNTS=apply(counts[,sapply(counts, is.numeric)],1, sum))
  }
  counts <- counts[order(counts$TOTAL_COUNTS, decreasing = T),]
  counts$TOTAL_COUNTS <- NULL
  
  # Sum all taxa below nr.taxa to Others
  if(length(unique(taxonomy[,key])) > nr.taxa+1){
    
    # Get others
    Others <- counts[(nr.taxa+1):nrow(counts),]
    if(missing(not.shown)){
      Others <- data.frame(dplyr::summarise_all(
        dplyr::group_by_if(Others[,c(1, which(sapply(Others, is.numeric)))], is.character), sum))
    } else {
      Others <- data.frame(dplyr::summarise_all(
        dplyr::group_by_if(Others[,c(which(colnames(Others) == 
                                             colnames(Others[,!(sapply(Others, is.numeric) | 
                                                                  colnames(Others) %in% not.shown)])[1]),
                                     which(sapply(Others, is.numeric)))], is.character), sum))
    }

    # Get counts except others
    counts <- counts[1:nr.taxa,] 
    
    # Add 'Others' only to branches with multiple taxa 
    Others[Others[,1] %in% unique(counts[,1]),1] <- 
      paste("ZZZZZZZZZZZZZZZZZZZZZZZZZZOther", Others[Others[,1] %in% unique(counts[,1]),1])

    # Merge both data sets
    counts <- merge(counts, Others, all = T)
    
    # Reorder columns
    counts <- counts[,c(colnames(taxonomy), colnames(counts)[!colnames(counts) %in% colnames(taxonomy)])]

    # Replace NAs in last columns, therefore iterate over rows
    for (row in seq(nrow(counts))) {
      if(any(is.na(counts[row,]))){
        column.nr.name <-  which(sapply(counts, is.character) & (!is.na(counts[row,])))
        counts[row, is.na(counts[row,])] <- counts[row, column.nr.name]
        counts[row, column.nr.name] <- sub("ZZZZZZZZZZZZZZZZZZZZZZZZZZOther ", "", counts[row, column.nr.name])
      }
    }
  }
  
  #------------------------------- Order columns alphabetically ------------------------------#
  # Order columns alphabetically
  for(column in rev(seq(ncol(counts)))){
    counts <- counts[order(counts[,column]),]
  }
  
  # Remove unwanted column
  if(!missing(not.shown)){
    counts <- counts[,!colnames(counts) %in% not.shown]
  }
 
  # Replace pattern ZZZZZZZZZZZZZZZZZZZZZZZZZZOther by Others
  counts <- data.frame(sapply(counts[,sapply(counts,is.character)], 
                              function(x) gsub("ZZZZZZZZZZZZZZZZZZZZZZZZZZOther", "Other", x)),
                       counts[,sapply(counts,is.numeric)])
  
  # Duplicate last taxonomic column to create a blank column
  counts <- data.frame(counts[,sapply(counts, is.character)], 
                       LastColumn=counts[,rev(which(sapply(counts, is.character)))[1]], 
                       counts[,sapply(counts, is.numeric)])
  
  #------------------------------ Assign unique letters to taxa ------------------------------#
  # Create a vector with letters 
  letters <- sort(c(sapply(LETTERS, sub, pattern = " ", 
                               x = sprintf("%2s", sapply(LETTERS, sub, pattern = " ", 
                                                         x = sprintf("%3s", LETTERS))))))
  
  # Get column names
  col.names <- colnames(counts)[sapply(counts, is.character)]
  
  # Assign a ID column
  counts$ID_COLUMN <- letters[seq(nrow(counts))]
  
  # Iterate over the taxonomic columns
  for(column in seq(col.names)){

    # Create copy with only taxonomic ranks
    counts.copy <- setNames(data.frame(counts[,seq(column)]), col.names[seq(column)])

    # Delete duplicated rows 
    counts.copy <- unique(counts.copy)
    
    # Add new column with IDs
    counts.copy[,paste("NEWWWWW", colnames(counts.copy)[column])] <- 
      paste(counts.copy[,column], letters[seq(nrow(counts.copy))])
    
    # Merge count table and the copy, order count table
    counts <- merge(counts, counts.copy, all = T)
    counts <- counts[order(counts$ID_COLUMN),]
    
    # Delete used IDs
    letters <- letters[-seq(nrow(counts.copy))]
  }
  
  # Delete duplicated column and ID column
  counts <- counts[,!colnames(counts) %in% c(col.names, "ID_COLUMN")]
  colnames(counts) <- sub("^NEWWWWW ", "", colnames(counts))
  counts <- data.frame(counts[,sapply(counts, is.character)], 
                       counts[,sapply(counts, is.numeric)])

  #------------------------------------- Get frequencies -------------------------------------#
  # Get frequencies for all taxonomic level, therefore iterate over character columns
  frequencies <- setNames(data.frame(matrix(ncol = 4, nrow = 0)) , c("node", "sum", "mean", "sd"))
  for (column in  which(sapply(counts, is.character))) {
    freq <- data.frame(dplyr::summarise_all(
      dplyr::group_by_if(counts[,c(column,which(sapply(counts, is.numeric)))], is.character), sum))
    
    # Assign ID column (letters)
    colnames(freq)[1] <- "ID"
    freq$ID <- sub(".* ", "", freq$ID)
    
    # Calculate sum, mean and sd across numeric columns
    freq <- data.frame(node=freq[,1], 
                       sum=apply(freq[,sapply(freq, is.numeric)],1, sum),
                       mean=apply(freq[,sapply(freq, is.numeric)],1, mean),
                       sd=apply(freq[,sapply(freq, is.numeric)],1, sd)) 

    frequencies <- rbind(frequencies, freq)
  }
  
  # Delete unused frequency column
  if(taxa.abundance == "mean"){
    frequencies$sum <- NULL
  } else {
    frequencies$mean <- NULL
    frequencies$sd <- NULL
  }
  
  # Assign percentage, if total is defined divide by total
  if(taxa.abundance == "mean"){
    if(missing(total)){
      frequencies$pct <- frequencies$mean/mean(colSums(counts[,sapply(counts, is.numeric)]))*100
      frequencies$pct_sd <- frequencies$sd/mean(colSums(counts[,sapply(counts, is.numeric)]))*100
    } else {
      frequencies$pct <- frequencies$mean/total*100
      frequencies$pct_sd <- frequencies$sd/total*100
    }
  } else {
    if(missing(total)){
      frequencies$pct <- frequencies$sum/sum(counts[,sapply(counts, is.numeric)])*100
    } else {
      frequencies$pct <- frequencies$sum/total*100
    }
  }
  
  #------------------------------------- Make node table -------------------------------------#
  # Create node table
  counts <- counts[,sapply(counts, is.character)]
  nodes <- unique(cbind(x=sort(rep(seq(ncol(counts)), nrow(counts))), stack(counts)))

  # Assign letters as new ID column to nodes table and delete letters from value
  nodes$ID <- sub(".* ", "", nodes$values)
  nodes$values <- sub(" ...$", "", nodes$values)

  #-------------------------------------- Assign colors --------------------------------------#
  # Create a copy for the count table
  counts.copy <- counts
  counts.copy$LastColumn <- NULL
  counts.copy <- counts.copy[!grepl("Other", counts.copy[,ncol(counts.copy)]),]

  #------------- Customized colors
  if(color.palette == "customize"){
    
    message(paste0("The sankey diagram has ", length(unique(counts.copy[,1])), " branches that split into the ", 
                   nr.taxa, " most abundant taxa of the taxonomic rank '", key, "'. "), 
            ifelse(any(table(counts[,1]) > 1), 
                   paste0(paste0(table(counts[,1])[table(counts[,1]) != 1]-1, " separate taxon/taxa belong to the branch '", 
                                 sub(" ...$", "", names(table(counts[,1]))[table(counts[,1]) != 1]), "'", collapse = " & "), ". "), ""), 
            ifelse(any(table(counts[,1]) == 1), 
                   paste0("No separate taxa are shown for the branche(s) ", 
                          paste(sub(" ...$", "", names(table(counts[,1]))[table(counts[,1]) == 1]), collapse = " & "), ". "), ""),
            paste0("\nYou are now asked to specify one color for each of the ", length(unique(counts.copy[,1])), 
                   " branches as well as one color for each of the ", nr.taxa, " most abundant taxa. "),
            "For the sankey diagram, color gradients are created between the colors of the individual taxa and the corresponding branch colors. ",
            "You can either use hexadecimal colors (e.g. #FF0000) or named colors (e.g. red).")
    
    #--------- Assign colors to the branches
    all.colors <- c()
    for (branch in unique(counts.copy[,1])) {
      
      # Get color for the branch
      while(TRUE){
        color <- readline(paste0("Please specify the color for the branch ", sub(" ...$", "", branch),":") )
        color <- gsub("\"", "", color)
        if(tryCatch(is.matrix(col2rgb(color)), error = function(e) FALSE)){
          break
        } else {
          message("The specified color is not valid, please use either hexadecimal colors (e.g. #FF0000) or named colors (e.g. red).")
        }
      }
      
      # Convert color to rgb color
      color <- rgb(col2rgb(color)[,1][1], col2rgb(color)[,1][2], col2rgb(color)[,1][3], maxColorValue=255)
      all.colors <- c(all.colors, color)
      names(all.colors)[length(all.colors)] <- branch
    }
    
    # Create a data frame with the branch colors
    all.colors <- setNames(data.frame(all.colors), "col")
    all.colors$ID <- rownames(all.colors)
    
    #--------- Get colors for the individual taxa
    
    # Iterate over the branches
    for (branch in unique(counts.copy[,1])) {
      
      #--- Assign colors for only one taxon per branch (no others)
      if(sum(counts.copy[,1] %in% branch) == 1 | sum(counts[,1] %in% branch) == 1){
        
        # Specify the printed text
        if(sum(counts[,1] %in% branch) == 1){
          text <- paste0("Please specify the color for last node (lowest taxonomic rank) of the from the branch ", sub(" ...$", "", branch), ":")
          
        } else {
          text <- paste0("Please specify the color for the taxon '", 
                         sub(" ...$", "", counts.copy[counts.copy[,1] == branch,ncol(counts.copy)]), "' from the branch ", sub(" ...$", "", branch), ":")
        }
        
        # Get color for the taxon of the branch
        while(TRUE){
          color <- readline(text)
          color <- gsub("\"", "", color)
          if(tryCatch(is.matrix(col2rgb(color)), error = function(e) FALSE)){
            break
          } else {
            message("The specified color is not valid, please use either hexadecimal colors (e.g. #FF0000) or named colors (e.g. red).")
          }
        }
        
        # Create a color sequence for the taxonomic ranks between the branch and the taxon color 
        palette <- colorRampPalette(colors=c(all.colors[all.colors$ID == branch,"col"], color))
        color <- palette(ncol(counts.copy))
        
        # Add the colors to the overall data frame
        all.colors <- rbind(all.colors, setNames(data.frame(color, unlist(counts.copy[counts.copy[,1] == branch,])), colnames(all.colors)))
        
        #--- Assign colors for multiple taxa per branch  
      } else {
        
        # Create open data frame for the colors of the taxa of one branch
        colors.branch <- setNames(data.frame(matrix(ncol = ncol(counts.copy), nrow = 0)), colnames(counts.copy))
        
        # Iterate across each taxon that belongs to the branch
        for (taxon in counts.copy[counts.copy[,1] == branch, ncol(counts.copy)]) {
          # Get color for the taxon
          while(TRUE){
            color <- readline(paste0("Please specify the color for the taxon '", sub(" ...$", "", taxon),
                                     "' from the branch ", sub(" ...$", "", branch), ":"))
            color <- gsub("\"", "", color)
            if(tryCatch(is.matrix(col2rgb(color)), error = function(e) FALSE)){
              break
            } else {
              message("The specified color is not valid, please use either hexadecimal colors (e.g. #FF0000) or named colors (e.g. red).")
            }
          }
          
          # Create a color sequence for the taxonomic ranks between the branch and the taxon color
          palette <- colorRampPalette(colors=c(all.colors[all.colors$ID == branch,"col"], color))
          color <- palette(ncol(counts.copy))
          
          # Add the colors to the color branch  data frame
          colors.branch <- rbind(colors.branch, setNames(data.frame(t(color)),colnames(counts.copy)))
        }
        
        # Create a copy from the count table for the branch
        counts.branch <- counts.copy[counts.copy[,1] %in% branch,]
        
        # Iterate over taxonomic ranks from highest to lowest
        for (column in rev(seq(ncol(counts.branch))[-1])) {
          
          # Check if taxa are duplicated at next higher taxonomic rank than the current
          if(any(duplicated(counts.branch[,column-1]))){
            
            # Iterate over duplicated taxa names
            for(taxon in unique(counts.branch[,column-1])){
              
              # Check if this taxon is duplicated
              if(sum(counts.branch[,column-1] == taxon) > 1){
                
                # Calculate mean of colors between the taxa
                means <- rowMeans(sapply(colors.branch[counts.branch[,column-1] == taxon,column], FUN=function(x) col2rgb(x)))
                means <- rgb(red=means[1], green=means[2], blue=means[3], maxColorValue=255)
                
                # Generate a new palette  
                palette <- colorRampPalette(colors=c(colors.branch[counts.branch[,column-1] == taxon,1][1], means))
                color <- palette(column)
                
                # Assign new colors to the branch color table
                for (i in seq(column-1)) {
                  colors.branch[counts.branch[,column-1] == taxon,i] <- color[i]
                }
                
              } 
            }
          } 
          
          # Add colors to overall color table
          all.colors <- rbind(all.colors, setNames(data.frame(colors.branch[,column],counts.branch[,column]), colnames(all.colors)))
        }
      }
    }
    
    # Remove duplicated rows
    all.colors <- unique(all.colors)
    all.colors$ID <- sub(".* ", "", all.colors$ID)
    
    #------------- Color paletts
  } else {
    # Assign colors
    all.colors <- c()
    for (row in seq(nrow(counts.copy))) {
      for (column in seq(ncol(counts.copy))) {
        all.colors <- c(all.colors, paste(colnames(counts.copy)[column], counts.copy[row,column]))
      }
    }
    all.colors <- unique(all.colors)
    
    # Assign colors
    all.colors <- data.frame(ID=all.colors,col=viridis::viridis(length(all.colors), option = color.palette))
    all.colors$ID <- sub(".* ","", all.colors$ID)
  }
  
  # Merge colors with node table
  nodes <- merge(nodes, all.colors, by = 'ID', all = TRUE)

  # Assign color to other and first and last column
  nodes[nodes$ind == "LastColumn", "col"] <- "white" 
  nodes[is.na(nodes$col), "col"] <- color.others
  
  # Merge node and frequencies table
  nodes <- merge(nodes, frequencies, by.x = 'ID', by.y = 'node', all = TRUE)

  # Assign labels
  if(labels=="none"){
    nodes$labels <- ""
  } else if(labels=="names"){
    nodes$labels <- nodes$values
  } else if(labels %in% c("abundance", "names.abundance")){
    if(taxa.abundance == "mean"){
      nodes$labels <- paste0(round(nodes$mean, decimal.places), " ± ", round(nodes$sd, decimal.places))#±
    } else {
      nodes$labels <- round(nodes$sum, decimal.places)
    }
  } else if(labels %in% c("rel.abundance", "names.rel.abundance")){
    if(taxa.abundance == "mean"){
      nodes$labels <- paste0(round(nodes$pct, decimal.places), " ± ", round(nodes$pct_sd, decimal.places),  "%")#±
    } else {
      nodes$labels <- paste0(round(nodes$pct, decimal.places), "%")
    }
  }
  if(labels %in% c("names.abundance", "names.rel.abundance")){
    nodes$labels <- paste0(nodes$values, " (", nodes$labels, ")")
  }

  # Assign labels
  nodes[(nodes$x != length(unique(nodes$x))-1) & (grepl("Other ", nodes$values)), "labels"] <- ""
  nodes[nodes$ind == "LastColumn", "labels"] <- ""

  #------------------------------------- Make edge table -------------------------------------#
  # Create edge table
  edges <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("N1", "N2", "Value"))
  for(i in seq(ncol(counts))[-1]){
    e <- setNames(unique(counts[,c(i-1,i)]), c("N1", "N2"))
    e$N1 <- sub(".* ","", e$N1)
    e$N2 <- sub(".* ","", e$N2)
    e <- merge(e, nodes[nodes$x==i,c("ID",taxa.abundance)], by.x = "N2", by.y = "ID", all.x = T)
    edges <- rbind(edges, setNames(data.frame(e[,c("N1", "N2", taxa.abundance)], 
                                              stringsAsFactors = F), c("N1", "N2", "Value")))
  }
  
  # Don show 'Others'
  if(show.others==F){
    edges <- edges[!(edges$N1 %in% nodes[grepl("^Other ", nodes$values),"ID"] | 
                       edges$N2 %in% nodes[grepl("^Other ", nodes$values),"ID"]),]
    nodes <- nodes[!grepl("^Other ", nodes$values),]
  }

  # Adjust node table
  n <- setNames(data.frame(nodes[,c("ID", "x", "col", "labels", "ind")], stringsAsFactors= FALSE),
                c("ID", "x", "col", "labels", "TaxonomicLevel"))
  n$TaxonomicLevel <- as.character(n$TaxonomicLevel)
  n[n$TaxonomicLevel == "LastColumn", "TaxonomicLevel"] <- ""
  n <- n[order(n$ID, decreasing = T),]
  n$ID <- factor(n$ID)
  
  # Create riverplot
  r <- riverplot::makeRiver(n, edges )
  
  # Clear node table
  nodes <- nodes[nodes$labels!="",]
  
  return(list(r=r, nodes=n, edges=edges))
}
