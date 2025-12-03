install.packages(c("FactoMineR", "dplyr", "openxlsx", "factoextra", "vegan", 
                   "tidyverse", "cowplot", "colorspace", "ggrepel", "lme4", 
                   "DHARMa", "fitdistrplus", "ggpubr", "emmeans", "TSdist", 
                   "gtools", "stringr"))


library("FactoMineR")
library(dplyr)
library(openxlsx)
library("factoextra")
library(vegan)
library(tidyverse)
library(cowplot)
library(colorspace)
library(ggrepel)
library(lme4)
library(DHARMa)
library(fitdistrplus)
library("ggpubr")
library(emmeans)
library(TSdist)
library(gtools)
library(stringr)

set.seed(12)
list <- list()


# Loop script till figure 4 and save results in "list"
for (a in 1:2) {

# Load the dataset (txt file)
xdata <- read.delim("xdata.txt", header=TRUE, stringsAsFactors = TRUE)

IDs <- xdata$Combination
tab <- table(IDs)
weights <- 1 / (length(tab) * tab[IDs])
weights <- as.numeric(weights)

#Generate random data for all FoCs
for (i in 8:343) {
  
  xdata[,i] <- as.factor(rbinom( 735 , 1 , .5))
}


############################################ Multiple Correspondence Analysis (MCA) ############################################

xdata <- xdata[,-c(1)] # For clearer results, we remove the recording names from the analysis

res.mca<-MCA(xdata,quanti.sup=c(5,6), quali.sup=c(1:4), row.w=weights) # Run the MCA

#################################################### Extract distances between utterances ####################################################

coord <- as.data.frame(res.mca$quali.sup$coord[1:26,]) # Extract the 5-dimension coordinates of each utterance type
distances <- as.matrix(vegdist(coord, method="euclidian")) # Measure the distance between each utterance type

# Create a function to assess whether a single call is part of a combination
is_in_combination <- function(cri, combination) {
  cri_parts <- unlist(strsplit(cri, "-"))
  any(cri_parts %in% unlist(strsplit(combination, "-")))
}

# Extract, for each possible pair of utterance types of the dataset, the distance between the two utterances and their relationship (the utterances are two single calls, two combinations, or one signal call is a constituent or not of the other combination) 
newdata <- data.frame(matrix(nrow = 351))
processed_pairs <- list()
k=1
for (i in 1:nrow(distances)){
  for (j in 1:ncol(distances)){
    
    pair <- sort(c(rownames(distances)[i], colnames(distances)[j]))
    pair_string <- paste(pair, collapse = "_")
    if (!(pair_string %in% processed_pairs)) {
      
      if ((nchar(rownames(distances)[i])==nchar(colnames(distances)[j])) & (nchar(rownames(distances)[i])==2)) {
        newdata$dist[k] <- distances[i,j]
        newdata$call_1[k] <- colnames(distances)[j]
        newdata$call_2[k] <- rownames(distances)[i]
        newdata$pair[k] <- pair_string
        newdata$comb[k] <- "single"
        k <- k+1
      } else if ((nchar(rownames(distances)[i])==nchar(colnames(distances)[j])) & (nchar(rownames(distances)[i])==5)) {
        newdata$dist[k] <- distances[i,j]
        newdata$call_1[k] <- colnames(distances)[j]
        newdata$call_2[k] <- rownames(distances)[i]
        newdata$pair[k] <- pair_string
        newdata$comb[k] <- "two_comb"
        k <- k+1
      } else if (is_in_combination(rownames(distances)[i], colnames(distances)[j])==TRUE) {
        newdata$dist[k] <- distances[i,j]
        newdata$call_1[k] <- colnames(distances)[j]
        newdata$call_2[k] <- rownames(distances)[i]
        newdata$pair[k] <- pair_string
        newdata$comb[k] <- "within_comb"
        k <- k+1
      } else {
        newdata$dist[k] <- distances[i,j]
        newdata$call_1[k] <- colnames(distances)[j]
        newdata$call_2[k] <- rownames(distances)[i]
        newdata$pair[k] <- pair_string
        newdata$comb[k] <- "outside_comb"
        k <- k+1
      }
      processed_pairs <- c(processed_pairs, pair_string) # To keep track of which pairs have been processed already
    }
  }
}

############################## The meaning of AB is derived from the meaning of A and B  ###################


xdata_1 <- rbind(newdata[which(newdata$comb=="within_comb"),], newdata[which(newdata$comb=="outside_comb"),]) # Substract the dataset and keep all the pairs one combination - one single type
xdata_1$comb <- as.factor(xdata_1$comb)
xdata_1$pair <- as.factor(xdata_1$pair)

# Indicate, for each pair, the name of the combination type
xdata_1$comb2 <- apply(xdata_1, 1, function(row) {
  if(nchar(row["call_1"]) > nchar(row["call_2"])) {
    return(row["call_1"])
  } else {return(row["call_2"])}
})
xdata_1$comb2 <- as.factor(xdata_1$comb2)

# Create the full model and check model assumptions
model1 <- lm(dist~comb*comb2, data=xdata_1)
simulationOutput <- simulateResiduals(fittedModel = model1)

model_null <- lm(dist~1, data=xdata_1) # Create a null model and compare it to the full model using a likelihood ratio test 

# Create two reduced models and compare them to the full model using a likelihood ratio test
model_comb <- lm(dist~comb, data=xdata_1) 
model_comb2 <- lm(dist~comb2, data=xdata_1)

model1_emm <- emmeans(model1, pairwise~comb|comb2) # Post hoc analysis: which combination is closer to its constituting calls than to other single calls of the repertoire

# Plot the results (Figure 4)
d <- as.data.frame(model1_emm$emmeans)
e <- as.data.frame(model1_emm$contrasts)
significant_contrasts <- e[e$p.value < 0.05, ]


  list[[a]] <- significant_contrasts
}

# Calculate porportion on runs with significant results
s <- sapply(list, nrow)
(length(list)-sum(s==0))/length(list)
