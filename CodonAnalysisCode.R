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
