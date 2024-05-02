association.frequencies<-function(edge.table, taxonomy, old.key="RowNames", new.key){
  
  #' @title Frequencies of undirected associations per taxonomic rank
  #'
  #' @description The function calculates the frequency of undirected negative and positive associations (not taking into account the edge weight) at the desired taxonomic level (specified by the argument \strong{new.key}). 
  #' This requires a network inference table and the corresponding taxonomy data. If numeric columns are given in the taxonomy table, the values are summed up at desired taxonomic level.
  #' \strong{An example:} The frequency of associations per Order is calculated, i.e. \strong{new.key} = "Order". Thus, the number of respectively positive and negative associations is counted for each Order pair. Further, if the taxonomy table contains a column with the number of reads per taxon (numerics), the reads are summed up at Order level.
  #' The function returns a list with the edge and taxonomy tables.
  #'
  #' @section Required packages:
  #' dplyr
  #'
  #' @param edge.table (\emph{data.frame/matrix}). A network table with 2 node columns (characters/factors) and 1 edge weight column (numeric, only used to differentiate between negative and positive associations). 
  #' @param taxonomy (\emph{data.frame/matrix}). Taxonomy table with taxonomic ranks as columns and taxa as rows. For character/factor columns the frequency of negative and positive associations is counted per \strong{new.key}, numeric columns are summed per \strong{new.key}.
  #' @param old.key (\emph{character/"RowNames"}). A character that specifies the column name of the \strong{old.key} of the taxonomy table. The \strong{old.key} contains the taxa names (IDs) that were also used for network inference (names of Node 1 & 2). In case the taxa names are used as row names of the \strong{taxonomy}, please set the argument to \emph{"RowNames"}. By default \emph{old.key = "RowNames"}.
  #' @param new.key (\emph{character}). A character that specifies the column name of the \strong{new.key} of the \strong{taxonomy}. The \strong{new.key} is the column name for which the frequency of the associations is calculated.
  
  
  #############################################################################################
  #################################### Check prerequisites ####################################
  
  #--------------------------- Check if needed arguments are valid ---------------------------#
  
  # Check if the 'edge.table' is provided
  if(missing(edge.table) | !(class(edge.table) %in% c("data.frame", "matrix"))){
    stop("Please provide a data.frame/matrix with 2 node columns (characters/factors) and 1 edge weight column (numeric).")
  } else if(ncol(edge.table)!=3){
    stop("Please provide a data.frame/matrix with only 3 columns: 2 node columns (characters/factors) and 1 edge weight column (numeric).")
  } else {
    # Format data
    edge.table <- data.frame(edge.table)
    edge.table <- dplyr::mutate_if(edge.table, is.factor, as.character)
    
    # Check if 2 character and one numeric column are persent
    if(ncol(dplyr::select_if(edge.table, is.character)) != 2 | ncol(dplyr::select_if(edge.table, is.numeric)) != 1){
      stop("Please provide a data.frame/matrix with 2 character/factor columns (2 node colums) and 1 numeric column (edge weights).")
    }
  }
  
  # Check if the 'taxonomy' table is provided
  if(missing(taxonomy) | !(class(taxonomy) %in% c("data.frame", "matrix"))){
    stop("Please provide a taxonomy table formatted as data.frame/matrix.")
  }
  taxonomy <- data.frame(dplyr::mutate_if(taxonomy, is.factor, as.character))
  
  # Check if the value of the argument 'old.key' is valid
  if(class(old.key) != "character"){
    stop("The argument 'old.key' is not valid. Please provide a character that specifies one column name of the taxonomy table or set the argument to 'RowNames'.")
  } else if(!(old.key %in% colnames(taxonomy) | old.key == "RowNames")){
    stop("The argument 'old.key' is not valid. The specified name is not found in the column names of the taxonomy table. Please provide one column name of the taxonomy table or set the argument to 'RowNames'.")
  } else {
    ifelse(old.key == "RowNames", old.key.column <- rownames(taxonomy), old.key.column <- taxonomy[,old.key])
    
    # Check if all entries of the Node 1 & 2 can be found in the old key column
    if(!all(unlist(dplyr::select_if(edge.table, is.character)) %in% old.key.column)){
      stop("Not all entries of the node columns can be found in the taxonomy table column spcified by the argument old.key.")
    }
  }
  
  # Check if the value of the argument 'new.key' is valid
  if(missing(new.key)){
    stop("The argument new.key is missing. Please specify one column name of the taxonomy table for which the frequency of the associations is calculated.")
  } else if(!new.key %in% colnames(taxonomy)){
    stop("The argument 'new.key' is not valid. Please specify one column name of the taxonomy table for which the frequency of the associations is calculated.")
  }

  #############################################################################################
  ################################### Calculate frequencies ###################################
  
  # Add column specifying old.key 
  if(old.key == "RowNames"){
    taxonomy$Old_key_column <- rownames(taxonomy)
  }else{
    taxonomy$Old_key_column <- taxonomy[,old.key]
  }
  
  # Reduce taxonomy, keep only taxa also present in the edge table
  taxonomy <- taxonomy[taxonomy$Old_key_column %in% unlist(edge.table[,colnames(dplyr::select_if(edge.table, is.character))]),]
  
  # Exclude all columns with higher taxonomic ranks than the new key column, keep numerics
  colnames.taxonomy <- colnames(taxonomy)[sapply(taxonomy,function(x){ length(unique(x))}) <= length(unique(taxonomy[,new.key]))]
  taxonomy <- taxonomy[,colnames(taxonomy) %in% c(colnames.taxonomy, colnames(dplyr::select_if(taxonomy, is.numeric)), "Old_key_column")]

  #------------------------ Create taxonomy table with new key column -----------------------#
  # Duplicate rows of the taxonomy shall be removed
  # but only characters/factors are taken into account to detect duplicate rows, not numerics. Numerics are summed up later.
  taxonomy.new <- dplyr::select_if(taxonomy, is.character)
  taxonomy.new$Old_key_column <- NULL
  
  # Removing duplicated rows of taxonomy table
  taxonomy.new <- taxonomy.new[!duplicated(taxonomy.new),]

  # Order table by new key column
  if(class(taxonomy.new) != "data.frame"){
    taxonomy.new <- setNames(data.frame(sort(taxonomy.new)), new.key)
  } else {
    taxonomy.new <- taxonomy.new[order(taxonomy.new[,new.key]),]
  }
  
  # Add number to the the new key column
  # This is necessary because different organisms may have the identical names on one taxonomic level
  taxonomy.new$New_key_column <- paste(seq(nrow(taxonomy.new)), taxonomy.new[,new.key], sep = "_")
  
  #----------------------------- Merge old and new taxonomy table ----------------------------#
  
  # Merge both taxonomy tables
  taxonomy <- dplyr::full_join(taxonomy, taxonomy.new)
  taxonomy <- data.frame(dplyr::mutate_if(taxonomy, is.factor, as.character))
  
  #------------------------------- Check if edges are undirected -----------------------------#
  
  # Extract names of Node columns
  nodes.column.1=colnames(dplyr::select_if(edge.table, is.character))[1]
  nodes.column.2=colnames(dplyr::select_if(edge.table, is.character))[2]
  
  # Check if the node combinations Node1_Node2 also exists as Node2_Node1
  if(any(paste(edge.table[,nodes.column.1], edge.table[,nodes.column.2], sep = "_") %in%
         paste(edge.table[,nodes.column.2], edge.table[,nodes.column.1], sep = "_"))){
    stop("Some of the node combinations exist as Node1_Node2 as well as Node2_Node1. This cannot occur in an undirected networks.")
  }
  
  #--------------------------------- Assign new Keys to Nodes --------------------------------#
  
  # Merge edge table and new taxonomy
  edge.table.new <- merge(x=edge.table, y=taxonomy, by.x = nodes.column.1, by.y ="Old_key_column")
  
  # Replace old keys with new
  edge.table.new[,nodes.column.1] <- edge.table.new$New_key_column
  
  # Only keep edge properties
  edge.table.new <- edge.table.new[,colnames(edge.table)]
  
  # Merge edge table and new taxonomy
  edge.table.new <- merge(x=edge.table.new, y=taxonomy, by.x = nodes.column.2, by.y ="Old_key_column")
  
  # Replace old keys with new
  edge.table.new[,nodes.column.2] <- edge.table.new$New_key_column
  
  # Only keep edge properties
  edge.table.new <- edge.table.new[,colnames(edge.table)]
  
  # Sort all node combinations alphabetically to avoid having the combination NodeA_NodeB as well as NodeB_NodeA
  for(row in seq(nrow(edge.table.new))){
    nodes <- sort(c(edge.table.new[row,nodes.column.1], edge.table.new[row, nodes.column.2]))
    edge.table.new[row,nodes.column.1] <- nodes[1]
    edge.table.new[row,nodes.column.2] <- nodes[2]
  }
  
  #------------------------- Count number of Associations per new Key ------------------------#
  # Gent name of edge weight column
  colname.edge.weights=colnames(dplyr::select_if(edge.table.new, is.numeric))
  
  # Add column to edge.table.new to distinguish between positive and negative associations
  edge.table.new$Association <- "negative"
  edge.table.new[dplyr::select_if(edge.table.new, is.numeric) == 0,"Association"] <- "neutral"
  edge.table.new[dplyr::select_if(edge.table.new, is.numeric) > 0,"Association"] <- "positive"
  
  # Edge weights are replaced by 1, since there is one association
  edge.table.new[,colname.edge.weights] <- 1
  
  # Sum associations per new key
  edge.table.new <- data.frame(dplyr::summarise_if(dplyr::group_by_if(edge.table.new, is.character), is.numeric, sum))
  
  #------------------------- Fit taxonomy data to the new Edge table -------------------------#
  
  # Delete old key column
  taxonomy$Old_key_column <- NULL
  
  # Remove duplicated rows by character/factor column, sum numeric column up
  if(ncol(dplyr::select_if(taxonomy, is.numeric)) == 0){
    taxonomy <- taxonomy[!duplicated(taxonomy),]
  } else {
    taxonomy <- data.frame(dplyr::summarise_if(dplyr::group_by_if(taxonomy, is.character), is.numeric, sum))
  }
  
  # Set new keys as row names of the taxonomy table
  rownames(taxonomy) <- taxonomy$New_key_column
  
  # Set Key at first position
  taxonomy <- cbind(Key=taxonomy$New_key_column, taxonomy[,!colnames(taxonomy) %in% "New_key_column"])
  
  return(list(Edges=edge.table.new, Taxonomy=taxonomy))
}