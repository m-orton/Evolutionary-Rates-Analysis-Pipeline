library(ape)
library(Biostrings) 
library(plotly)
library(foreach)

# Renaming columns
colnames(dfL1L2TempArthropoda)[14] <- "L1"
colnames(dfL1L2TempArthropoda)[15] <- "L2"
colnames(dfL1L2TempArthropoda)[16] <- "Outgroup"

# 2nd codon position - using regular expressions and taking every 3rd nucleotide char, 2nd codon = 1st position in aligned sequences
dfL1L2TempArthropoda$L1_2nd <- gsub("(.)..", "\\1", dfL1L2TempArthropoda$L1) 
dfL1L2TempArthropoda$L2_2nd <- gsub("(.)..", "\\1", dfL1L2TempArthropoda$L2) 
dfL1L2TempArthropoda$Outgroup_2nd <- gsub("(.)..", "\\1", dfL1L2TempArthropoda$Outgroup) 

# 3rd codon position - same but 3rd codon = 2nd position in aligned sequences
# Chopping off first char to get to second
dfL1L2TempArthropoda$L1_3rd <- substr(dfL1L2TempArthropoda$L1, 2 , nchar(dfL1L2TempArthropoda$L1))
dfL1L2TempArthropoda$L2_3rd <- substr(dfL1L2TempArthropoda$L2, 2 , nchar(dfL1L2TempArthropoda$L2))
dfL1L2TempArthropoda$Outgroup_3rd <- substr(dfL1L2TempArthropoda$Outgroup, 2 , nchar(dfL1L2TempArthropoda$Outgroup))

# 3rd codon extraction
dfL1L2TempArthropoda$L1_3rd <- gsub("(.)..", "\\1", dfL1L2TempArthropoda$L1_3rd) 
dfL1L2TempArthropoda$L2_3rd <- gsub("(.)..", "\\1", dfL1L2TempArthropoda$L2_3rd) 
dfL1L2TempArthropoda$Outgroup_3rd <- gsub("(.)..", "\\1", dfL1L2TempArthropoda$Outgroup_3rd) 

# Conversion to dnaBIN format before pairwise distance calcs
# 2nd
L1_2ndbin <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% as.DNAbin(DNAStringSet(dfL1L2TempArthropoda$L1_2nd[i]))
L2_2ndbin <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% as.DNAbin(DNAStringSet(dfL1L2TempArthropoda$L2_2nd[i]))
Outgroup_2ndbin <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% as.DNAbin(DNAStringSet(dfL1L2TempArthropoda$Outgroup_2nd[i]))

# For ingroup distance
L1L2_2ndbin <- foreach(i=1:length(L1_2ndbin)) %do% append(L1_2ndbin[[i]], L2_2ndbin[[i]])

# For dist from L1 and L2 to outgroup
L1_outgroup2ndbin <- foreach(i=1:length(L1_2ndbin)) %do% append(L1_2ndbin[[i]], Outgroup_2ndbin[[i]])
L2_outgroup2ndbin <- foreach(i=1:length(L1_2ndbin)) %do% append(L2_2ndbin[[i]], Outgroup_2ndbin[[i]])

# 3rd
L1_3rdbin <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% as.DNAbin(DNAStringSet(dfL1L2TempArthropoda$L1_3rd[[i]]))
L2_3rdbin <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% as.DNAbin(DNAStringSet(dfL1L2TempArthropoda$L2_3rd[[i]]))
Outgroup_3rdbin <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% as.DNAbin(DNAStringSet(dfL1L2TempArthropoda$Outgroup_3rd[[i]]))

# For ingroup distance
L1L2_3rdbin <- foreach(i=1:length(L1_3rdbin)) %do% append(L1_3rdbin[[i]], L2_3rdbin[[i]])

# For dist from L1 and L2 to outgroup
L1_outgroup3rdbin <- foreach(i=1:length(L1_3rdbin)) %do% append(L1_3rdbin[[i]], Outgroup_3rdbin[[i]])
L2_outgroup3rdbin <- foreach(i=1:length(L1_3rdbin)) %do% append(L2_3rdbin[[i]], Outgroup_3rdbin[[i]])

# inGroupDistances and outgroup distance calculation

# 2nd
IngroupDist_2nd <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% dist.dna(L1L2_2ndbin[[i]], model = "TN93", pairwise.deletion = TRUE)
dfL1L2TempArthropoda$IngroupDist_2nd <- as.numeric(IngroupDist_2nd)
OugroupDistL1_2nd <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% dist.dna(L1_outgroup2ndbin[[i]], model = "TN93", pairwise.deletion = TRUE)
dfL1L2TempArthropoda$OugroupDistL1_2nd <- as.numeric(OugroupDistL1_2nd)
OugroupDistL2_2nd <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% dist.dna(L2_outgroup2ndbin[[i]], model = "TN93", pairwise.deletion = TRUE)
dfL1L2TempArthropoda$OugroupDistL2_2nd <- as.numeric(OugroupDistL2_2nd)

# 3rd
IngroupDist_3rd <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% dist.dna(L1L2_3rdbin[[i]], model = "TN93", pairwise.deletion = TRUE)
dfL1L2TempArthropoda$IngroupDist_3rd <- as.numeric(IngroupDist_3rd)
OugroupDistL1_3rd <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% dist.dna(L1_outgroup3rdbin[[i]], model = "TN93", pairwise.deletion = TRUE)
dfL1L2TempArthropoda$OugroupDistL1_3rd <- as.numeric(OugroupDistL1_3rd)
OugroupDistL2_3rd <- foreach(i=1:nrow(dfL1L2TempArthropoda)) %do% dist.dna(L2_outgroup3rdbin[[i]], model = "TN93", pairwise.deletion = TRUE)
dfL1L2TempArthropoda$OugroupDistL2_3rd <- as.numeric(OugroupDistL2_3rd)

# Filter out very small ingroup distance values
zFilter_2nd <- which(dfL1L2TempArthropoda$IngroupDist_2nd > 0) 
zFilter_3rd <- which(dfL1L2TempArthropoda$IngroupDist_3rd > 0) 

# Filter out indeterminant values - only present in 3rd codon
indFilter_3rd <- which(is.na(dfL1L2TempArthropoda$IngroupDist_3rd) | is.na(dfL1L2TempArthropoda$OugroupDistL2_3rd) | is.na(dfL1L2TempArthropoda$OugroupDistL1_3rd))

dfL1L2TempArthropoda_2nd <- dfL1L2TempArthropoda[zFilter_2nd,]
dfL1L2TempArthropoda_3rd <- dfL1L2TempArthropoda[zFilter_3rd,]
dfL1L2TempArthropoda_3rd <- dfL1L2TempArthropoda[-indFilter_3rd,]

# First creating starting branch lengths
dfL1L2TempArthropoda_2nd$startingBranchLength_2nd <- dfL1L2TempArthropoda_2nd$IngroupDist_2nd * 0.5
dfL1L2TempArthropoda_3rd$startingBranchLength_3rd <- dfL1L2TempArthropoda_3rd$IngroupDist_3rd * 0.5

# Then determining half outgroup distance difference for each lineage
dfL1L2TempArthropoda_2nd$halfOutGroupDiffL1 <- (dfL1L2TempArthropoda_2nd$OugroupDistL1_2nd - dfL1L2TempArthropoda_2nd$OugroupDistL2_2nd) / 2
dfL1L2TempArthropoda_2nd$halfOutGroupDiffL2 <- (dfL1L2TempArthropoda_2nd$OugroupDistL2_2nd - dfL1L2TempArthropoda_2nd$OugroupDistL1_2nd) / 2

dfL1L2TempArthropoda_3rd$halfOutGroupDiffL1 <- (dfL1L2TempArthropoda_3rd$OugroupDistL1_3rd - dfL1L2TempArthropoda_3rd$OugroupDistL2_3rd) / 2
dfL1L2TempArthropoda_3rd$halfOutGroupDiffL2 <- (dfL1L2TempArthropoda_3rd$OugroupDistL2_3rd - dfL1L2TempArthropoda_3rd$OugroupDistL1_3rd) / 2

# Estimated Branch Lengths per Lineage - 2nd
dfL1L2TempArthropoda_2nd$branchLengthEstimationL1 <- dfL1L2TempArthropoda_2nd$startingBranchLength_2nd + dfL1L2TempArthropoda_2nd$halfOutGroupDiffL1
dfL1L2TempArthropoda_2nd$branchLengthEstimationL2 <- dfL1L2TempArthropoda_2nd$startingBranchLength_2nd + dfL1L2TempArthropoda_2nd$halfOutGroupDiffL2

# Estimated Branch Lengths per Lineage - 3rd
dfL1L2TempArthropoda_3rd$branchLengthEstimationL1 <- dfL1L2TempArthropoda_3rd$startingBranchLength_3rd + dfL1L2TempArthropoda_3rd$halfOutGroupDiffL1
dfL1L2TempArthropoda_3rd$branchLengthEstimationL2 <- dfL1L2TempArthropoda_3rd$startingBranchLength_3rd + dfL1L2TempArthropoda_3rd$halfOutGroupDiffL2

# Eliminate negative branch lengths 
negative_2nd <- which(dfL1L2TempArthropoda_2nd$branchLengthEstimationL1 > 0 & dfL1L2TempArthropoda_2nd$branchLengthEstimationL2 > 0)
dfL1L2TempArthropoda_2nd <- dfL1L2TempArthropoda_2nd[negative_2nd,]
negative_3rd <- which(dfL1L2TempArthropoda_3rd$branchLengthEstimationL1 > 0 & dfL1L2TempArthropoda_3rd$branchLengthEstimationL2 > 0)
dfL1L2TempArthropoda_3rd <- dfL1L2TempArthropoda_3rd[negative_3rd,]

# Standardization of temperature difference
dfL1L2TempArthropoda_2nd$standardDiffTemp2 <- dfL1L2TempArthropoda_2nd$TempDiff / sqrt(dfL1L2TempArthropoda_2nd$branchLengthEstimationL1 + dfL1L2TempArthropoda_2nd$branchLengthEstimationL2)
dfL1L2TempArthropoda_3rd$standardDiffTemp3 <- dfL1L2TempArthropoda_3rd$TempDiff / sqrt(dfL1L2TempArthropoda_3rd$branchLengthEstimationL1 + dfL1L2TempArthropoda_3rd$branchLengthEstimationL2)

# Standardization of lat difference
dfL1L2TempArthropoda_2nd$standardDiffLat2 <- dfL1L2TempArthropoda_2nd$latDelta.x / sqrt(dfL1L2TempArthropoda_2nd$branchLengthEstimationL1 + dfL1L2TempArthropoda_2nd$branchLengthEstimationL2)
dfL1L2TempArthropoda_3rd$standardDiffLat3 <- dfL1L2TempArthropoda_3rd$latDelta.x / sqrt(dfL1L2TempArthropoda_3rd$branchLengthEstimationL1 + dfL1L2TempArthropoda_3rd$branchLengthEstimationL2)

# Branch length diff (absolute)
dfL1L2TempArthropoda_2nd$branchLengthDiff2 <- (abs(dfL1L2TempArthropoda_2nd$branchLengthEstimationL1 + dfL1L2TempArthropoda_2nd$branchLengthEstimationL2))
dfL1L2TempArthropoda_3rd$branchLengthDiff3 <- (abs(dfL1L2TempArthropoda_3rd$branchLengthEstimationL1 + dfL1L2TempArthropoda_3rd$branchLengthEstimationL2))

# Standardization of branch length diff
dfL1L2TempArthropoda_2nd$standardDiffBranch2 <- dfL1L2TempArthropoda_2nd$branchLengthDiff2 / sqrt(dfL1L2TempArthropoda_2nd$branchLengthEstimationL1 + dfL1L2TempArthropoda_2nd$branchLengthEstimationL2)
dfL1L2TempArthropoda_3rd$standardDiffBranch3 <- dfL1L2TempArthropoda_3rd$branchLengthDiff3 / sqrt(dfL1L2TempArthropoda_3rd$branchLengthEstimationL1 + dfL1L2TempArthropoda_3rd$branchLengthEstimationL2)

# Signage of branch length diff 
foreach(i=1:nrow(dfL1L2TempArthropoda_2nd)) %do% 
  if(dfL1L2TempArthropoda_2nd$standardDiffBranch[i] < 0){
    dfL1L2TempArthropoda_2nd$standardDiffBranch2[i] <- dfL1L2TempArthropoda_2nd$standardDiffBranch2[i] * -1
  }

foreach(i=1:nrow(dfL1L2TempArthropoda_3rd)) %do% 
  if(dfL1L2TempArthropoda_3rd$standardDiffBranch[i] < 0){
    dfL1L2TempArthropoda_3rd$standardDiffBranch3[i] <- dfL1L2TempArthropoda_3rd$standardDiffBranch3[i] * -1
  }

# For second codon position
plotPublicationTemp2nd <- plot_ly(
  dfL1L2TempArthropoda_2nd, x = ~standardDiffTemp2) %>%
  add_markers(y = ~standardDiffBranch2, color = ~order_name.x.x, colors = c("#009e73", "#D55E00", "#0072B2", "#CC79A7","#7F7F7F")) %>% 
  add_lines(x = ~standardDiffTemp2, y = fitted(lm(standardDiffBranch2 ~ 0 + standardDiffTemp2, data = dfL1L2TempArthropoda_2nd)))

plotPublicationLat2nd <- plot_ly(
  dfL1L2TempArthropoda_2nd, x = ~standardDiffLat2) %>%
  add_markers(y = ~standardDiffBranch2, color = ~order_name.x.x, colors = c("#009e73", "#D55E00", "#0072B2", "#CC79A7","#7F7F7F")) %>% 
  add_lines(x = ~standardDiffLat2, y = fitted(lm(standardDiffBranch2 ~ 0 + standardDiffLat2, data = dfL1L2TempArthropoda_2nd)))

# For third codon position
plotPublicationTemp3rd <- plot_ly(
  dfL1L2TempArthropoda_3rd, x = ~standardDiffTemp3) %>%
  add_markers(y = ~standardDiffBranch3, color = ~order_name.x.x, colors = c("#009e73", "#D55E00", "#0072B2", "#CC79A7","#7F7F7F")) %>% 
  add_lines(x = ~standardDiffTemp3, y = fitted(lm(standardDiffBranch3 ~ 0 + standardDiffTemp3, data = dfL1L2TempArthropoda_3rd)))

plotPublicationLat3rd <- plot_ly(
  dfL1L2TempArthropoda_3rd, x = ~standardDiffLat3) %>%
  add_markers(y = ~standardDiffBranch3, color = ~order_name.x.x, colors = c("#009e73", "#D55E00", "#0072B2", "#CC79A7","#7F7F7F")) %>% 
  add_lines(x = ~standardDiffLat3, y = fitted(lm(standardDiffBranch3 ~ 0 + standardDiffLat3, data = dfL1L2TempArthropoda_3rd)))
