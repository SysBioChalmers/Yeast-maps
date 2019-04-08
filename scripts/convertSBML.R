##################################################
## Project: Yeast-maps
## Script purpose: Convert yeast maps from SBML
##                 format into .NET
## Date: 2019-02-28
## Author: Angelo limeta
##################################################

# == INSTALL PACKAGES ==

if (!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
require(BiocManager)
BiocManager::install("SBMLR", version = "3.8")
require(SBMLR)

#install.packages("xml2")
require(xml2)
install.packages("readr")
require(readr)

# == LOAD DATA ==

# Load subsystem maps in SBML format
PATH = "~/Documents/PhD/Yeast-maps/"
setwd(PATH)
subsystem = "pyruvate metabolism"
PATH_xml = paste(PATH,"SBMLfiles/",subsystem,".xml",sep = "")
subsystem_xml = read_xml(PATH_xml)
subsystem_list = as_list(subsystem_xml)

# Load example .NET file
# This is the correct formatting for omix compatibility 
PATH_csv = paste(PATH,"netfiles/golgi_example.net",sep = "")
example_csv = as.data.frame(read_delim(PATH_csv, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE))
colnames(example_csv)

# Load yeast GEM
library(readxl)
# Reactions
yeastGEM_rxn = data.frame(read_excel("scripts/yeast_8.3.xlsx"))
# Metabolites
yeastGEM_met = read_excel(paste(PATH,"scripts/yeastGEM.xlsx",sep = ""), sheet = "METS")
colnames(yeastGEM_met)[1:2] = c("Metabolite.name", "Metabolite.description")

# Load geneGroups
geneGroups <- data.frame(read_csv("scripts/geneGroups.csv"))

# == CONSTRUCT PATHWAY DATA FRAME ==

pathway_df = data.frame("P", subsystem, subsystem, -100, -1, -1)
colnames(pathway_df) = colnames(example_csv)

# == CONSTRUCT REACTION DATA FRAME ==

# Setup empty data frame for storing data
nrxn = as.integer(length(subsystem_list$sbml$model$listOfReactions)) # Retrieve amount of metabolites
rxn_df = data.frame(matrix(NA, nrow = nrxn, ncol = 8)) # Create empty data frame

# "R" denotes that each row is a reaction
rxn_df[,1] = rep("R",nrxn)

# Store reaction IDs
for (i in 1:nrxn) {
  rxn_df[i,2] = attributes(subsystem_list$sbml$model$listOfReactions[[i]])$metaid
  rxn_df[i,3] = attributes(subsystem_list$sbml$model$listOfReactions[[i]])$metaid
}

# Store Pathway ID
rxn_df[,4] = rep(subsystem,nrxn)

# Store reaction coordinates
nalias = length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases)
for (i in 1:nrxn) {
  # Store base reactant/product id
  base_reac = attributes(subsystem_list$sbml$model$listOfReactions[[i]]$annotation$extension$baseReactants$baseReactant)$alias
  base_prod = attributes(subsystem_list$sbml$model$listOfReactions[[i]]$annotation$extension$baseProducts$baseProduct)$alias
  # Store base rectant/product coordinates
  for (j in 1:nalias) {
    if(attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id == base_reac){
      base_reac_x = as.numeric(attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$x)
      base_reac_y = as.numeric(attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$y)
    }
    if(attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id == base_prod){
      base_prod_x = as.numeric(attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$x)
      base_prod_y = as.numeric(attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$y)
    }
  }
  # Average the x and y coordinate values for base reactant/product and set this as reaction coordinates
  # NOTE, this can map multiple reactions to one position if they share the same base reactant/product
  rxn_df[i,5] = (base_reac_x+base_prod_x)/2
  rxn_df[i,6] = (base_reac_y+base_prod_y)/2
}

# Store reversibility information
yeastGEM_rxn_bounds = yeastGEM_rxn[,c(1,5,6)] # Keep only reaction IDs, lower (LB) and upper bounds (UB)
rownames(yeastGEM_rxn_bounds) = yeastGEM_rxn_bounds[,1]
yeastGEM_rxn_bounds = yeastGEM_rxn_bounds[,-1]
rxn_rev = rowSums(yeastGEM_rxn_bounds) # Reversible reactions have -1000 and 1000 as LB/UB, i.e sums to zero if reversible
for (i in 1:length(rxn_rev)) {
  if(rxn_rev[i] == 0){
    rxn_rev[i] = "true"
    }
  else{
    rxn_rev[i] = "false"
  }
}
for (i in 1:nrxn) {
  rxn_df[i,7] = as.character(rxn_rev[rxn_df[i,2]])
}

# Store reaction type
rxn_df[,8] = rep("normal",nrxn) 

# == CONSTRUCT GENE DATA FRAME ==

gene_ids = data.frame(ID = 1:length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases), 
                      Name = 1:length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases))
for (i in 1:length(gene_ids[,1])) {
  gene_ids[i,1] = attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[i]])$species
}
# Store unique gene names
for (i in 1:length(gene_ids[,1])) {
  for (j in 1:length(subsystem_list$sbml$model$listOfSpecies)) {
    if(gene_ids[i,1] == attributes(subsystem_list$sbml$model$listOfSpecies[[j]])$id){
      gene_ids[i,2] = attributes(subsystem_list$sbml$model$listOfSpecies[[j]])$name
    }
  }
}
# Remove metabolites/proteins
gene_ids = gene_ids[grepl("g_", gene_ids$Name),]
# Remove duplicate values
gene_ids = unique(gene_ids)

# Create list object with all unique gene IDs
gene_list = list()
for (i in 1:nrow(gene_ids)) {
  gene_list[[gene_ids[i,1]]] = list()
  attributes(gene_list[[gene_ids[i,1]]])$Name = gene_ids[i,2]
}

# Create sublists for each unique gene ID which contain all duplicate IDs
for (i in 1:length(gene_list)) {
  for (j in 1:length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases)) {
    if(names(gene_list)[i] == attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$species){
      gene_list[[i]][[attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id]] = list()
      attributes(gene_list[[i]][[attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id]])$x = attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$x
      attributes(gene_list[[i]][[attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id]])$y = attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$y
    }
  }
}
# Create sublists for each duplicate ID which contain all reactions it participates in
for (i in 1:length(gene_list)) {
  for (j in 1:length(gene_list[[i]])) {
    for (k in 1:length(subsystem_list$sbml$model$listOfReactions)) {
      # Store genes for each reaction
      gr = c()
      for (m in 1:length(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfModification)) {
        if(attributes(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfModification[[m]])$aliases %in% names(gene_list[[i]])[j]){
          gr = c(gr,attributes(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfModification[[m]])$aliases)
        }
      }
      # Attribute reactions to the base reactant
      if(names(gene_list[[i]])[j] %in% gr){
        gene_list[[i]][[j]][[gsub("e", "_", attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid)]] = list()
      }
    }
  }
}

# Count up how many genes are in the subsystem
ngene = 0
for (i in 1:length(gene_list)) {
  str_genes = geneGroups$string_agg[which(geneGroups$id %in% names(gene_list[[i]][[1]]))]
  if(length(str_genes) == 1){
    ngene = ngene + length(strsplit(str_genes,split = ", ")[[1]])
  }
}

# Create empty data frame
gene_df = data.frame(matrix(NA, nrow = ngene, ncol = 6))

# Fill data frame
row_index = 1
for (i in 1:length(gene_list)) {
  str_genes = geneGroups$string_agg[which(geneGroups$id %in% names(gene_list[[i]][[1]]))]
  if(length(str_genes) == 1){
  str_genes = strsplit(str_genes,split = ", ")[[1]]
    for (m in 1:length(str_genes)) {
      gene_df[row_index,1] = "E"
      gene_df[row_index,2] = str_genes[m]
      gene_df[row_index,3] = str_genes[m]
      gene_df[row_index,4] = as.character(as.numeric(attributes(gene_list[[i]][[1]])$x) + round(runif(1,-5,5),digits = 1))
      gene_df[row_index,5] = as.character(as.numeric(attributes(gene_list[[i]][[1]])$y) + round(runif(1,-5,5),digits = 1))
      gene_df[row_index,6] = names(gene_list[[i]][[1]])
      row_index = row_index + 1
    }
  }
}

# == CONSTRUCT METABOLITE DATA FRAME ==

# Construct a nested list object containing the necessary info about reaction/metabolite connectivity and features

# Store unique metabolite IDs
met_ids = data.frame(ID = 1:length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases), 
                      Name = 1:length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases))
for (i in 1:length(met_ids[,1])) {
  met_ids[i,1] = attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[i]])$species
}

# Store unique metabolite names
for (i in 1:length(met_ids[,1])) {
  for (j in 1:length(subsystem_list$sbml$model$listOfSpecies)) {
    if(met_ids[i,1] == attributes(subsystem_list$sbml$model$listOfSpecies[[j]])$id){
      met_ids[i,2] = attributes(subsystem_list$sbml$model$listOfSpecies[[j]])$name
    }
  }
}

# Keep only metabolite names, i.e. not genes/proteins
# Metabolites do not contain the "_" character, so we can filter using it
met_ids = met_ids[!grepl("_", met_ids$Name),]
# Remove duplicate values
met_ids = unique(met_ids)

# Create list object with all unique metabolite IDs
met_list = list()
for (i in 1:nrow(met_ids)) {
  met_list[[met_ids[i,1]]] = list()
  attributes(met_list[[met_ids[i,1]]])$Name = met_ids[i,2]
}

# Create sublists for each unique metabolite ID which contain all duplicate IDs
for (i in 1:length(met_list)) {
  for (j in 1:length(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases)) {
    if(names(met_list)[i] == attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$species){
      met_list[[i]][[attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id]] = list()
      attributes(met_list[[i]][[attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id]])$x = attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$x
      attributes(met_list[[i]][[attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]])$id]])$y = attributes(subsystem_list$sbml$model$annotation$extension$listOfSpeciesAliases[[j]]$bounds)$y
    }
  }
}

# Create sublists for each duplicate ID which contain all reactions it participates in
for (i in 1:length(met_list)) {
  for (j in 1:length(met_list[[i]])) {
    for (k in 1:length(subsystem_list$sbml$model$listOfReactions)) {
      # Store base reactants for each reaction
      br = c()
      for (m in 1:length(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$baseReactants)) {
        br = c(br,attributes(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$baseReactants[[m]])$alias)
      }
      
      # Attribute reactions to the base reactant
      if(names(met_list[[i]])[j] %in% br){
        met_list[[i]][[j]][[gsub("e", "_", attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid)]] = list()
        attributes(met_list[[i]][[j]][[attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid]])$edgetype = "outEdge" # Set direction of arrow
        attributes(met_list[[i]][[j]])$collapsed = "false" # Not a cofactor
        }
      
      # Store base products for each reaction
      bp = c()
      for (m in 1:length(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$baseProducts)) {
        bp = c(bp,attributes(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$baseProducts[[m]])$alias)
      }
      # Attribute reactions to the base product
      if(names(met_list[[i]])[j] %in% bp){
        met_list[[i]][[j]][[gsub("e", "_", attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid)]] = list()
        attributes(met_list[[i]][[j]][[attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid]])$edgetype = "inEdge" # Set direction of arrow
        attributes(met_list[[i]][[j]])$collapsed = "false" # Not a cofactor
      }
      
      # Store cofactor reactants for each reaction
      lr = c()
      for (m in 1:length(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfReactantLinks)) {
        lr = c(lr,attributes(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfReactantLinks[[m]])$alias)
      }
      # Attribute reactions to the cofactor reactant
      if(names(met_list[[i]])[j] %in% lr){
        met_list[[i]][[j]][[gsub("e", "_", attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid)]] = list()
        attributes(met_list[[i]][[j]][[attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid]])$edgetype = "outEdge" # Set direction of arrow
        attributes(met_list[[i]][[j]])$collapsed = "true" # Cofactor
      }
      
      # Store cofactor products for each reaction
      lp = c()
      for (m in 1:length(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfProductLinks)) {
        lp = c(lp,attributes(subsystem_list$sbml$model$listOfReactions[[k]]$annotation$extension$listOfProductLinks[[m]])$alias)
      }
      # Attribute reactions to the cofactor product
      if(names(met_list[[i]])[j] %in% lp){
        met_list[[i]][[j]][[gsub("e", "_", attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid)]] = list()
        attributes(met_list[[i]][[j]][[attributes(subsystem_list$sbml$model$listOfReactions[[k]])$metaid]])$edgetype = "inEdge" # Set direction of arrow
        attributes(met_list[[i]][[j]])$collapsed = "true" # Cofactor
      }
      
    }
  }
}



for (i in 1:length(met_list)) {
  for (j in 1:length(met_list[[i]])) {
    if(names(met_list[[i]])[j] %in% br){print(met_list[[i]][j])}
  }
}

# Fuse incorrectly labeled duplicates, i.e. s123 "oxaloacetate [c]" and s96 "oxaloacetate[c]"
index_duplicate = c() # Create empty vector for storing indexes of duplicates
for (i in 1:length(met_list)) { # Loop over each metabolite
  isduplicate = FALSE # This is to ensure that incorrectly labeled metabolites that do not have duplicates do not get removed 
  if(grepl(" \\[",attributes(met_list[[i]])$Name)){ # Find all metabolites with a space before the square brackets
    duplicate_met = gsub(" \\[","\\[",attributes(met_list[[i]])$Name) # Store the correct label
    print(duplicate_met)
    for (j in 1:length(met_list)) { # Loop over each metabolite again to identify the correctly labeled duplicate
      if(attributes(met_list[[j]])$Name == duplicate_met){
        isduplicate = TRUE
        index_duplicate = c(index_duplicate, i) # Save index for the incorrectly labeled duplicate
        duplicate_name = gsub(" \\[","\\[",attributes(met_list[[j]])$Name)
        met_list[[j]] = append(met_list[[j]],met_list[[i]]) # Append all sub lists of the incorrectly labeled duplicate to the correct one
        attributes(met_list[[j]])$Name = duplicate_name 
        break
      }
    }
    if(!isduplicate){ # Relabel incorrectly labeled metabolites
      attributes(met_list[[i]])$Name = duplicate_met
    }
  }
}
met_list = met_list[-index_duplicate] # Remove all incorrectly labeled duplicates

# Add yeast model IDs
for (i in 1:length(met_list)) {
  attributes(met_list[[i]])$ModelID = yeastGEM_met$`REPLACEMENT ID`[which(yeastGEM_met$Metabolite.description == attributes(met_list[[i]])$Name)]
}

# Store number of rows in data frame
df_rows = 0
for (i in 1:length(met_list)) {
  df_rows = df_rows +1
  df_rows = df_rows + length(met_list[[i]])
  for (j in 1:length(met_list[[i]])) {
    df_rows = df_rows + length(met_list[[i]][[j]])
  }
}

met_df = data.frame(matrix(NA, nrow = 0, ncol = 13)) # Create empty data frame, 9 cols + 4 more for splinepoints

row_index = 1
for (i in 1:length(met_list)) {
  print(attributes(met_list[[i]])$Name)
  for (j in 1:length(met_list[[i]])) {
    if(j == 1){
      met_df[row_index,1] = "MI"
      met_df[row_index,2] = attributes(met_list[[i]])$ModelID
      met_df[row_index,3] = attributes(met_list[[i]])$Name
      met_df[row_index,4] = attributes(met_list[[i]][[j]])$collapsed
      met_df[row_index,5] = attributes(met_list[[i]][[j]])$x
      met_df[row_index,6] = attributes(met_list[[i]][[j]])$y
      met_df[row_index,7] = "Center"
      met_df[row_index,8] = "0.0"
      met_df[row_index,9] = "0.0"
    }
    else{
      met_df[row_index,1] = "MI"
      met_df[row_index,2] = NA
      met_df[row_index,3] = NA
      met_df[row_index,4] = attributes(met_list[[i]][[j]])$collapsed
      met_df[row_index,5] = attributes(met_list[[i]][[j]])$x
      met_df[row_index,6] = attributes(met_list[[i]][[j]])$y
      met_df[row_index,7] = "Center"
      met_df[row_index,8] = "0.0"
      met_df[row_index,9] = "0.0"
    }
    row_index = row_index + 1
    for (k in 1:length(met_list[[i]][[j]])) {
      met_df[row_index,1] = "L"
      met_df[row_index,2] = attributes(met_list[[i]][[j]][[k]])$edgetype
      met_df[row_index,3] = names(met_list[[i]][[j]])[k]
      met_df[row_index,4] = attributes(met_list[[i]][[j]])$collapsed
      met_df[row_index,5] = NA
      met_df[row_index,6] = NA
      row_index = row_index + 1
    }
  }
}

# Convert data frame values to character
pathway_df[] <- lapply(pathway_df, as.character)
rxn_df[] <- lapply(rxn_df, as.character)
gene_df[] <- lapply(gene_df, as.character)
met_df[] <- lapply(met_df, as.character)


# Append all data frames and export TSV file
library(dplyr)
dfToExport = bind_rows(pathway_df,rxn_df,gene_df,met_df)

write.table(dfToExport, file = paste(PATH,"netfiles/",subsystem,".net",sep = ""), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

