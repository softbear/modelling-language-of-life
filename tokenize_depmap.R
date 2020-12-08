library("readr")
library("dplyr")
library("tidyverse")

setwd("~/Documents/orthogonal_dependency")
path <- getwd()
gene_expression <- read_csv(paste(path,"CCLE_expression.csv",sep="/"))
Achilles_gene_effect <- read_csv(paste(path,"Achilles_gene_effect.csv",sep="/"))
colnames(gene_expression)[1] <- "DepMap_ID"
#expression_and_essentiality <- merge(gene_expression,Achilles_gene_effect, by = "DepMap_ID")
gene_expression <- column_to_rownames(gene_expression,var = "DepMap_ID")
Achilles_gene_effect <- column_to_rownames(Achilles_gene_effect,var = "DepMap_ID")

common_cell_lines <- intersect(rownames(gene_expression),rownames(Achilles_gene_effect))

common_expression <- gene_expression[which(rownames(gene_expression) %in% common_cell_lines),]
common_essentiality <- Achilles_gene_effect[which(rownames(Achilles_gene_effect) %in% common_cell_lines),]


#expression <- matrix(data = NA, nrow = length(common_cell_lines),ncol = ncol(gene_expression)-1)
#expression <- vector("list", length(common_cell_lines))
#essentiality <- matrix(data = NA, nrow = length(common_cell_lines), ncol = ncol(Achilles_gene_effect)-1)
#essentiality <- vector("list", length(common_cell_lines))

#common_expression <- merge(gene_expression,common_cell_lines,by = intersect(unlist(gene_expression[,1]),common_cell_lines))

gene_name_to_id <- function(symbol){
  return(strsplit(symbol,"[()]")[[1]][2])
}

gene_name_to_symbol_and_id <- function(name){
  symbol <- strsplit(name," ")[[1]][1]
  id <- strsplit(name,"[()]")[[1]][2]
  return(list(symbol,id))
}



#debugging gene_name_to_id()
test <- colnames(gene_expression[,1:10])
sapply(test,gene_symbol_to_id)


colnames(common_expression) <- sapply(colnames(common_expression),gene_name_to_id)
colnames(common_essentiality) <- sapply(colnames(common_essentiality),gene_name_to_id)


gene_name_to_symbol_and_id <- function(name){
  symbol <- strsplit(name," ")[[1]][1]
  id <- strsplit(name,"[()]")[[1]][2]
  return(c(symbol,id))
}
#debugging gene_name_to_symbol_and_id()
test <- colnames(gene_expression[,1:10])
sapply(test,gene_name_to_symbol_and_id)

#symbol_id_conversion 
gene_names <- unique(c(colnames(gene_expression),colnames(Achilles_gene_effect)))
conversion <- t(sapply(gene_names,gene_name_to_symbol_and_id))
write.csv(conversion,file="gene_symbol_id_conversion.csv")





getOrderedExpression <- function(row){
  #remove index 1 which is the cell line name
  row <- row[-which(row<1)]
  decr_order <- order(row,decreasing = TRUE)
  res <- (row[decr_order])
  #print(res)
  return(names(res))
  
}
getOrderedEssentiality <- function(row){
  #remove index 1 which is the cell line name
  row <- row[-which(row>=-0.5)]
  incr_order <- order(row,decreasing = FALSE)
  res <- (row[incr_order])
  #print(res)
  return(names(res))
  
}


#debugging the getOrderedGenes() function
test <- common_expression[1:5,1:10]
apply(test,1,getOrderedExpression)

test <- common_essentiality[1:5,1:10]
apply(test,1,getOrderedEssentiality)

expression <- apply(common_expression,1,getOrderedExpression)
essentiality <- apply(common_essentiality,1,getOrderedEssentiality)

concatenate_list_into_string <- function(lst) {
  return(paste(lst,collapse = " "))
}

#Debug concatenate_list_into_string()
test <- expression[[1]][1:10]
concatenate_list_into_string(test)
test <- expression[1:3]
lapply(test,concatenate_list_into_string)

concat_expression <- unlist(lapply(expression,concatenate_list_into_string))
concat_essentiality <- unlist(lapply(essentiality,concatenate_list_into_string))


write.csv(concat_expression,"concat_expression.csv")
write.csv(concat_essentiality,"concat_essentiality.csv")
