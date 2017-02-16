##################
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# There is a copy of the GNU General Public License
# along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

###################
# Authors and Contributions:

# This program was authored primarily by Matthew Orton.
# Contributions by Jacqueline May for lines 347-378,
# editing, formatting and testing of this program.
# Contributions by David Lee for lines 1822-1843, 
# editing and formatting of this program.
# Contributions by Winfield Ly for editing and formatting of this program.
# Collaboration of Sarah Adamowicz (University of Guelph) 
# in designing the analyses for this program as well as editing, formatting 
# and testing of this program.

###################
# Program Purpose (Class-Level Analyses):

# This program will allow for the generation of latitudinally separated 
# sister pairings and associated outgroupings from nearly any taxa 
# (provided they have a suitable reference sequence 
# and are not too large) and geographical region found on the BOLD API in one 
# streamlined R pipeline! In this iteration of the pipeline, data is translated 
# directly from a BOLD tsv file of the user's choosing to a dataframe in R. 
# The generated sister pairs and outgroups can also be written to a csv or tsv, 
# and the file will appear in the current working directory of R. Additionally, 
# binomial and Wilcoxon tests can be performed on the relative outgroup 
# distances of the generated pairings, and a plot of the resultant relative 
# outgroup distances can be generated for each class. A world map visualizing 
# the latitudinally separated pairings for the taxon being run can also be 
# generated using using the R visualization tool plotly.

# Larger taxa including various phyla can be run and will get broken down into 
# classes, and sister pairs are sought within classes by default. 

# ***This version of the program is intended for the larger taxa over 10000
# BIN's in size.***

###################
# Guidelines and Tips on Using this Program:

# The versioning of R and R Studio is important and can influence the results 
# of this program. At this time we can only guarantee full function of this 
# program using R Studio versions 0.99.903 or 1.0.44 and R versions 3.3.1 
# or 3.2.4.

# It is highly recommended that you use the IDE R Studio to run this program 
# as this will provide a much more user friendly interface than just using R.

# This program is limited to using the BOLD API for download of taxonomic and 
# sequence data. The taxa being tested with this program must use BIN 
# identifiers as a means of identification; for instance, plants on the BOLD API
# cannot be run using this program since they do not use BIN identifiers.

# At this time, it is best not to run Arthropoda, Insecta, or Chordata 
# in their entirety since these are all very large and computationally demanding
# and should be broken down into smaller sub-taxa. As well be aware that the 
# larger insect orders including Diptera, Hymenoptera and Lepidoptera are large 
# enough that they require an additional code component that divides the dataset 
# into four separate geographical regions: North America, South America, 
# Australia and Eurasia/Africa. This code component will be located in 
# section 5 of the program.

# It is highly recommended that you save your workspaces frequently with large 
# taxa in case of unexpected crashes within R Studio.

# Larger taxa consume a great deal of working memory and R is very memory 
# intensive, so be mindful of the amount of available working memory you have 
# available on your PC/Mac. As a general rule, every 1000 BIN's processed using 
# this program requires roughly 1.5 Gb of working memory. 
# It's probably a good idea to ensure RStudio is the only memory intensive 
# application running when you run this program with a larger taxon. It is also
# important to note that some taxa are so large that running R Studio locally 
# isnt really possible due to computational demands. However there is a great
# online version of R Studio using an EC2 instance on Amazon Web 
# Services (cloud computing instance). There is a guide here that will show 
# how to use this service: 
# http://strimas.com/r/rstudio-cloud-1/

# Entries on BOLD which do not at least have class-level classification 
# will be filtered out since the program cannot properly categorize these. 

# You do not have to worry about entries missing certain pieces of data, 
# for example latitude coordinates, sequence data etc.
# The script will filter these records out automatically.

# There are two options for parsing the tsv. 
# You can either download the tsv directly from the BOLD API, 
# or you can use a tsv you have previously downloaded. You can also save your 
# workspace once a TSV is downloaded so you don't have to download it again.

# When testing a taxon for the first time, you will have to ensure that you have
# a suitable reference sequence for it and ensure that it is inserted into the 
# dfRefSeq dataframe in Section 5 of the code. Ensure that you keep the entire 
# sequence in one line; breaking it up on to separate lines will add a new line 
# character.

# The reference sequence is a high-quality sequence from a given taxon, which is
# used as the basis for the first alignment step. As well, the reference 
# sequence is used to trim all sequences to the same length for analysis. To 
# insert a sequence to the dfRefSeq dataframe, simply add another taxon and 
# sequence in quotation marks in the refseq dataframe commented in the 
# reference sequence section. 

# Some tips for using the BOLD API since this is what is used to grab 
# the relevant data needed from BOLD:
# To see details on how to use the BOLD API, 
# go to http://www.boldsystems.org/index.php/resources/api?type=webservices
# We can add additional restrictions to the url, for example 
# &instituiton=Biodiversity Institute of Ontario|York University or 
# &marker=COI-5P if you want to specifiy an institution.
# geo=all in the url means global but geo can be geo=Canada, for example.
# We can also use | modifier in the url. For example, &geo=Canada|Alaska 
# would give data for both Canada and Alaska, 
# or taxon=Aves|Reptilia would yield data for both Aves and Reptilia.

#################
# Important Dataframes:

# dfPVal gives the p-values for both the binomial test and Wilcoxon tests for 
# all classes (that have pairings) and each class separately.

# dfPairingResultsSummary shows a more user-friendly summation of the 
# pairing results.

# dfRelativeDist shows the relative distances to the outgroup for each pairing, 
# including pseudoreplicates.

# dfPairingResultsL1L2 is a more detailed finalized dataframe that contains all 
# of the finalized pairings and outgroupings and all relevant details for them.

# dfPairingResultsL1 and dfPairingResultsL2 represent dataframes for each 
# lineage of each pairing, respectively.

# dfInitial is the dataframe first produced by the import from BOLD and is 
# filtered by latitude, bin_uri, N content, gap content, and sequence length.

# dfOutGroupL1 and L2 contain the associated outgroupings only 
# (for each lineage), but each one does have a column indicating the pairing 
# with which it is associated. 

# dfCentroid contains centroid sequences for all BINs with more than one member. 

# dfNonCentroid contains all BINs that only have one member and thus do not need 
# centroid sequences.

# dfAllSeq is simply dfNonCentroid and dfCentroid combined together in one 
# dataframe.

# dfGeneticDistanceStack is all genetic distance matrices concatenated into one 
# long column of values. It is used to grab bin_uri's and distances for each 
# pairing and assist in determining outgroup bin_uri's and distances.

# dfRefSeq shows various taxa with a suitable reference sequence that has been 
# found for them.

# dfPseudoRep shows each separate ingroup pairing number of each set of 
# pseudoreplicates along with its respective relative outgroup distance.

# dfPseudoRepAverage shows the averaged distances for each set of 
# pseudoreplicates in dfPseudoRep.

##################
# Important Variables:

# alignment2 will show a truncated alignment of each class before divergent 
# sequences are removed. (For example, typing alignment2[1] will show the first 
# alignment performed.) Note that the console will not show the full alignment, 
# to view the full alignment you must convert and export to FASTA format. 
# To export to FASTA format and view the full alignment, please see the commands
# in section 6 of the code.

# alignmentFinal will show a truncated final alignment of each class after 
# divergent sequences are removed. (For example, typing alignmentFinal[1] will 
# show the first alignment performed.) Note that the console will not show the 
# full alignment, to view the full alignment you must convert and export to 
# FASTA format. To export to FASTA format and view the full alignment, please 
# see the commands in section 7 of the code.

# alignmentFinalTrim will show a truncated final alignment of each class that is 
# also trimmed by the reference sequence. (For example, typing 
# alignmentFinalTrim[1] will show the first trimmed alignment performed.)
# Note that the console will not show the full alignment, to view the full 
# alignment you must convert to FASTA format. To export to FASTA format and view 
# the full alignment, please see the commands in section 8 of the code.

# binomialTestOutGroup shows the results of the binomial test.

# wilcoxonTestOutgroup shows the results of the Wilcoxon test.

# mapLayout is a variable that will allow for customization of the world map 
# for map plotting using the visualization tool Plotly.

#################
# Packages Required:

# Note that once you have installed the packages (first time running program), 
# you only have to run the libraries again each time you open up R Studio or 
# restart R Studio.
# Therefore, remove the "#" symbol in front of the lines for installing the 
# packages the first time running the script. 

# We need the foreach package for several functions that require iteration 
# over dataframe rows.
# install.packages("foreach")
library(foreach)
# For genetic distance determination using the TN93 model, the ape package is
# required.
# install.packages("ape")
library(ape)
# The readr package will speed up parsing of the tsv file using the 
# read_tsv function.
# install.packages("readr")
library(readr)
# For sequence alignments we need the biostrings (DNAStringSet function) 
# and muscle libraries, as follows:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("muscle")
library(Biostrings) 
library(muscle)
# For the calculation of overlapping latitude regions we need the 
# Desctools package. 
# The Desctools package allows us to use the Overlap function which can 
# calculate overlap regions between BINs based on the maximum and minimum 
# latitudinal ranges of each BIN.
# install.packages("DescTools")
library(DescTools)
# Adding the data tables package for data table merging of outgroup BIN data 
# with pairing results BIN data in the pairing result section. 
# install.packages("data.table")
library(data.table)
# For plotting of relative outgroup distances between lineages we will also need 
# ggplot2.
# install.packages ("ggplot2") 
require(ggplot2)
# dplyr is required if this is not already installed.
# install.packages("dplyr")
library(dplyr)
# As a pre-requisite package to plotly we need the colorspace package.
# install.packages("colorspace")
library(colorspace)
# Prerequisite for the plotly package.
# install.packages("jsonlite")
library(jsonlite)
# plotly package used for map plotting functionality.
# install.packages("plotly")
library(plotly)

#################
# R Commands:

# Section 1: TSV Parsing
# In this section the intial TSV of the desired taxon is parsed from the BOLD 
# database API.

# First we download the TSV and convert it into a dataframe. The URL below is 
# what is modified by the user and will determine the taxon, geographic region, 
# etc. Example: taxon=Aves$geo=all, see above in Guidelines and Tips for more
# information on how to use the BOLD API.

# The read_tsv function has been modified to select only certain columns to save
# on downloading time:
dfInitial <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=&geo=all&format=tsv")[,
   c('recordID', 'bin_uri','phylum_taxID','phylum_name',
     'class_taxID','class_name','order_taxID',
     'order_name','family_taxID','family_name',
     'subfamily_taxID','subfamily_name','genus_taxID',
     'genus_name','species_taxID','species_name',
     'lat','lon','nucleotides','markercode')]

# If you want to run pre downloaded BOLD TSV's to avoid downloading of the same 
# tsv multiple times, this will let you choose a path to that TSV and parse. 
# Uncomment and run these commands:

# tsvParseDoc <- file.choose()
# dfInitial <- read_tsv(tsvParseDoc)[,
#  c('recordID', 'bin_uri','phylum_taxID','phylum_name',
#    'class_taxID','class_name','order_taxID',
#    'order_name','family_taxID','family_name',
#    'subfamily_taxID','subfamily_name','genus_taxID',
#    'genus_name','species_taxID','species_name',
#    'lat','lon','nucleotides','markercode')]

# ***Keep in mind you can also save your R workspace to avoid redownloading 
# BOLD TSV's.***

###################
# Section 2: Dataframe Filtering and Reorganization

# In this section, the dataframe is filtered according to latitude, 
# sequence quality, and sequence length. 
# The initial dataframe is also reorganized.

# Renaming record id column.
colnames(dfInitial)[1] <- "record_id"

# Removing sequences with marker codes other than COI-5P.
containCOI <- grep( "COI-5P", dfInitial$markercode)
dfInitial<-dfInitial[containCOI,]

# Removing sequences with no latitude values. We are filtering according to lat 
# since we only really need lat for the analysis.
containLat <- grep( "[0-9]", dfInitial$lat)
dfInitial<-dfInitial[containLat,]

# Next, we have to convert the lat column to num instead of chr type. This will 
# become important later on for median latitude determination.
latNum <- with(dfInitial, as.numeric(as.character(lat))) 
dfInitial$latNum <- latNum

# We will do lon as well for map plotting using plotly.
lonNum <- with(dfInitial, as.numeric(as.character(lon))) 
dfInitial$lonNum <- lonNum

# We next identify records missing a BIN assignment and eliminate rows with 
# missing bin_uri's since the presence of a BIN designation is an indicator of 
# sequence length and quality. As well, we use BINs later in the analysis.
# This is achieved by using grep by colon, since every record with a BIN 
# identifier will have this.
containBin <- grep( "[:]", dfInitial$bin_uri)
dfInitial <- dfInitial[containBin,]

# We next get rid of any records that don't have sequence data.
# Sometimes there are records that bear a BIN but for which the sequence was 
# subsequently deleted. This could occur if a record had a sequence and the 
# record holder later deleted the sequence, e.g. due to suspected contamination.
containNucleotides <- grep( "[ACGT]", dfInitial$nucleotides)
dfInitial <- dfInitial[containNucleotides,]

# Next, we filter out high gap content and N content to ensure high sequence 
# quality.

# First, we need to convert nucleotides to chr type.
dfInitial$nucleotides <- with(dfInitial, (as.character(nucleotides))) 
# Cut off starting Ns and gaps (large portions of Ns and gaps occur at the start 
# of a sequence).
startNGap <- sapply(regmatches(dfInitial$nucleotides, gregexpr("^[-N]", 
                    dfInitial$nucleotides)), length)
startNGap <- foreach(i=1:nrow(dfInitial)) %do%
  if (startNGap[[i]] > 0) {
    split <- strsplit(dfInitial$nucleotides[i], "^[-N]+")
    dfInitial$nucleotides[i] <- split[[1]][2]
  }
# Cut off ending Ns and gaps (large portions of Ns and gaps occur at the end of 
# a sequence).
endNGap <- sapply(regmatches(dfInitial$nucleotides, gregexpr("[-N]$",
                  dfInitial$nucleotides)), length)
endNGap <- foreach(i=1:nrow(dfInitial)) %do%
  if (endNGap[[i]] > 0) {
    split <- strsplit(dfInitial$nucleotides[i], "[-N]+$")
    dfInitial$nucleotides[i] <- split[[1]][1]
  }
# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalNGap <- sapply(regmatches(dfInitial$nucleotides, gregexpr("[-N]", 
                       dfInitial$nucleotides)), length)
# We then loop through each sequence to see if the number of Ns or gaps is 
# greater than 1% (0.01) of the total sequence length.
internalNGap <- foreach(i=1:nrow(dfInitial)) %do% 
  which((internalNGap[[i]] / nchar(dfInitial$nucleotides[i]) > 0.01 ))
nGapCheck <- sapply( internalNGap , function (x) length( x ) )
nGapCheck <- which(nGapCheck > 0)
# Subset out these higher gap and N content sequences.
dfInitial <- dfInitial[-nGapCheck, ]

# Filter out sequences less than 640 bp and greater than 1000 bp since these 
# sequence length extremes can interfere with the alignment,
# and this also helps to standardize sequence length against the reference 
# sequences, for consistency in subsequent analyses.
sequenceLengths <- nchar(gsub("-", "", dfInitial$nucleotides))
sequenceLengthCheck <- which(sequenceLengths > 1000 | sequenceLengths < 640)
dfInitial <- dfInitial[-sequenceLengthCheck,]

# Modifying BIN column slightly to remove "BIN:"
dfInitial$bin_uri <- substr(dfInitial$bin_uri, 6 , 13)

# Dataframe Reorganization
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum_taxID","phylum_name",
                           "class_taxID","class_name","order_taxID","order_name"
                           ,"family_taxID","family_name","subfamily_taxID",
                           "subfamily_name","genus_taxID","genus_name",
                           "species_taxID","species_name","nucleotides",
                           "latNum","lonNum")])
###################
# Section 3: BIN Stats and Median Latitude/Longitude Determination per BIN

# In this section, the median latitude is determined for each BIN as well as 
# other important pieces of information including number of record ids within 
# each BIN and the latitudinal min and max of each BIN.

# First we can make a smaller dataframe with the columns we want for each BIN: 
# bin_uri, latnum, lonnum, record_id, and nucleotides.
dfBinList <- 
  (dfInitial[,c("record_id","bin_uri","latNum","lonNum","nucleotides")])

# We convert latitudes to absolute values before median latitude values are 
# calculated on dfBinList.
dfBinList$latNumAbs <- abs(dfBinList$latNum) 

# Lon remains untouched as it is not needed for the analysis.

# Create groupings by BIN with each grouping representing a different bin_uri.
# Each element of this list represents a BIN with subelements representing the 
# various columns of the initial dataframe created, and the information is 
# grouped by BIN. 
binList <- lapply(unique(dfBinList$bin_uri), function(x) 
  dfBinList[dfBinList$bin_uri == x,])

# Now determine the median latitude for each BIN based on absolute values.
medianLatAbs <- sapply( binList , function(x) median( x$latNumAbs ) )

# Calculating median lat based upon the original latitude values, 
# for mapping purposes only.
medianLatMap <- sapply( binList , function(x) median( x$latNum ) )

# We also need a median longitude for map plotting.
medianLon <- sapply( binList , function(x) median( x$lonNum ) )

# We can also take a few other important pieces of data regarding each BIN using 
# sapply, including number of record_ids within each BIN and the latitudinal min 
# and max of each BIN.
latMin <- sapply( binList , function(x) min( x$latNum ) )
latMax <- sapply( binList , function(x) max( x$latNum ) )
binSize <- sapply( binList , function (x) length( x$record_id ) )

# Dataframe of our median lat values. This will be used in our final dataframe.
dfLatLon <- data.frame(medianLatAbs)

# Adding bin_uri, latMin, latMax, binSize and medianLatMap to dataframe. 
# medianLatMap is the median latitude without conversion to an absolute value 
# and is used purely for mapping purposes.
dfLatLon$bin_uri <- c(unique(dfInitial$bin_uri))
dfLatLon$medianLon <- c(medianLon)
dfLatLon$latMin <- c(latMin)
# Convert to absolute value for latMin and latMax. 
dfLatLon$latMin <- abs(dfLatLon$latMin)
dfLatLon$latMax <- c(latMax)
dfLatLon$latMax <- abs(dfLatLon$latMax)
dfLatLon$binSize <- c(binSize)
dfLatLon$medianLatMap <- c(medianLatMap)

# Merging LatLon to BinList for the sequence alignment step.
dfBinList <- merge(dfBinList, dfLatLon, by.x = "bin_uri", by.y = "bin_uri")

# Also reordering dfLatLon by bin_uri for a step later on.
dfLatLon <- dfLatLon[order(dfLatLon$bin_uri),]

###################
# Section 4: Selecting a Centroid Sequence Per BIN

# A single Barcode Index Number (BIN) can contain many record ids and sequences.
# In order to simplify the analyses, we select one sequence per BIN.
# To do this we are determining the centroid sequence for each BIN.
# Centroid Sequence: BIN sequence with minimum average pairwise distance to all 
# other sequences in a given BIN.

# First we have to subset dfBinList to find BINs with more than one member since 
# these BINs will need centroid sequences.
largeBin <- which(dfBinList$binSize > 1)
# If there is at least one BIN with more than one member, then a dataframe 
# dfCentroid will be created with those BINs.
if (length(largeBin) > 0) {
  dfCentroid <- dfBinList[largeBin,]
  
  # Also subset dfLatLon down to number of BINs in dfInitial
  dfLatLon <- subset(dfLatLon, dfLatLon$bin_uri %in% dfCentroid$bin_uri)
  
  # Also need to find the number of unique BINs in dfCentroid. 
  binNumberCentroid <- unique(dfCentroid$bin_uri)
  binNumberCentroid <- length(binNumberCentroid)
  
  # We also have to create another separate dataframe with BINs that only have 
  # one member called dfNonCentroid.
  dfNonCentroid <- dfBinList[-largeBin, ]
  
  # We then take the dfCentroid sequences and break it down into a list with 
  # each element being a unique BIN. 
  largeBinList <- lapply(unique(dfCentroid$bin_uri), 
                         function(x) dfCentroid[dfCentroid$bin_uri == x,])
  
  # Extract record id from each BIN.
  largeBinRecordId <- sapply( largeBinList , function(x) ( x$record_id ) )
  
  # Convert all of the sequences in the largeBinList to dnaStringSet format for 
  # the alignment step.
  dnaStringSet1 <- 
    sapply( largeBinList, function(x) DNAStringSet(x$nucleotides) )
  
  # name DNAStringSet with the record ids.
  for (i in seq(from = 1, to = binNumberCentroid, by = 1)) {
    names(dnaStringSet1[[i]]) <- largeBinRecordId[[i]]
  }
  
  # Multiple sequence alignment using the muscle package.
  # Refer here for details on package: 
  # https://www.bioconductor.org/packages/release/bioc/html/muscle.html
  # Run a multiple sequence alignment on each element of the dnaStringSet1 list
  # Using diags = TRUE with Muscle command to speed up alignment of each BIN.
  # ***For larger taxa, it is highly recommended to use a maxiter 
  # setting of 2 to reduce alignment times.***
  alignment1 <- foreach(i=1:binNumberCentroid) %do% 
    muscle(dnaStringSet1[[i]], maxiters = 2, diags = TRUE, gapopen = -3000)
  
  # We can then convert each alignment to DNAbin format.
  dnaBINCentroid <- 
    foreach(i=1:binNumberCentroid) %do% as.DNAbin(alignment1[[i]])
  
  # Then, we perform genetic distance determination with the TN93 model on each 
  # DNAbin list.
  geneticDistanceCentroid <- foreach(i=1:binNumberCentroid) %do% 
    dist.dna(dnaBINCentroid[[i]], model = "TN93", as.matrix = TRUE, 
             pairwise.deletion = TRUE)
  
  # The centroid sequence can be determined from the distance matrix alone.
  # It is the sequence in a bin with minimum average pairwise distance to all 
  # other sequences in its BIN.
  centroidSeq <- foreach(i=1:binNumberCentroid) %do% 
    which.min(rowSums(geneticDistanceCentroid[[i]]))
  centroidSeq <- unlist(centroidSeq)
  centroidSeq <- names(centroidSeq)
  centroidSeq <- as.numeric(centroidSeq)
  
  # Subset dfCentroid by the record ids on this list.
  dfCentroid <- subset(dfCentroid, record_id %in% centroidSeq)
  
  # Append our noncentroid to our centroid sequences. We will make a new 
  # dataframe for this - dfAllSeq. 
  # Now we have all of the sequences we need for the next alignment of all 
  # sequences, with one sequence representing each BIN.
  dfAllSeq <-  rbind(dfCentroid, dfNonCentroid)
  # We can then merge this with dfInitial to get all of the relevant data we 
  # need, merging to dfInitial to gain all the relevant taxanomic data.
  dfAllSeq <- merge(dfAllSeq, dfInitial, by.x = "record_id", by.y = "record_id")
  # Renaming and reorganizing the dataframe.
  dfAllSeq <- 
    (dfAllSeq[,c("bin_uri.x","binSize","record_id","phylum_taxID","phylum_name",
                 "class_taxID","class_name","order_taxID","order_name",
                 "family_taxID","family_name","subfamily_taxID",
                 "subfamily_name","genus_taxID","genus_name","species_taxID",
                 "species_name","nucleotides.x","medianLatAbs","medianLatMap",
                 "latMin","latMax","medianLon")])
  colnames(dfAllSeq)[1] <- "bin_uri"
  colnames(dfAllSeq)[18] <- "nucleotides"
  # Adding an index column to reference later with Pairing Results dataframe.
  dfAllSeq$ind <- row.names(dfAllSeq)
  # Deleting any possible duplicate entries.
  dfAllSeq <- (by(dfAllSeq, dfAllSeq["bin_uri"], head, n = 1))
  dfAllSeq <- Reduce(rbind, dfAllSeq)
  
} else {
  # Else if there are no BINs with more than one member then we would simply 
  # merge latlon with initial to get dfAllSeq and thus all of the sequences we 
  # want.
  dfAllSeq <- merge(dfLatLon, dfInitial, by.x = "bin_uri", by.y = "bin_uri")
  dfAllSeq <- 
    (dfAllSeq[,c("bin_uri.x","binSize","record_id","phylum_taxID","phylum_name",
                 "class_taxID","class_name","order_taxID","order_name",
                 "family_taxID","family_name","subfamily_taxID",
                 "subfamily_name","genus_taxID","genus_name","species_taxID",
                 "species_name","nucleotides.x","medianLatAbs","medianLatMap",
                 "latMin","latMax","medianLon")])
  dfAllSeq$ind <- row.names(dfAllSeq)
}

# To simplify the number of dataframes, can remove dfLatLon and dfBinList since 
# these are redundant.
rm(dfBinList) 
rm(dfLatLon)
gc()

###################
# Section 5: Creation of Reference Sequence Dataframe

# In this section we establish reference sequences that will be used in the 
# final alignment for appropriate trimming of the alignment to a standard 
# sequence length of 620 bp.

# Currently, the script is finetuned for reference sequences of 658 bp of 
# the barcode region of the cytochrome c oxidase subunit I (COI) gene. Minor 
# modifications to the alignment settings or trimming amount would be needed to 
# use sequences of different lengths or different genes.

# We selected reference sequences that fell between the primers for COI 
# by Folmer et al. (1994): 

# Folmer O, M, Black WH, Lutz R, Vrijenhoek R (1994) DNA primers for 
# amplification of mitochondrial cytochrome C oxidase subunit I from metazoan 
# invertebrates. Molecular Marine Biology and Biotechnology 3: 294-299.

# In DNA barcoding studies, many taxonomic groups are amplified by these primers
# or by variants of these primers that bind at the same position. For most taxa,
# the reference sequences were exactly 658 bp in length. For certain taxa 
# (e.g. Bivalvia), the length differed from 658 bp due to amino acid indels, but
# the same start and end point, as verified using an amino acid alignment, 
# was used.

# Manual input of reference sequences into dfRefSeq dataframe
dfRefSeq <- data.frame(taxa = c("",""),
                       nucleotides = c("",
                                       ""))
colnames(dfRefSeq)[2] <- "nucleotides"
dfRefSeq$nucleotides <- as.character(dfRefSeq$nucleotides)

# Symmetrical trimming of the references to a standard 620 bp from 658 bp.
# A different trimming length could be used, depending upon the distribution of 
# sequence lengths in a particular taxon.
dfRefSeq$nucleotides <- 
  substr(dfRefSeq$nucleotides, 20, nchar(dfRefSeq$nucleotides) - 19)

# Check of sequence length.
dfRefSeq$seqLength <- nchar(dfRefSeq$nucleotides)

# Subset dfAllSeq by entries in the reference sequence dataframe.
dfAllSeq <- subset(dfAllSeq, dfAllSeq$class_name %in% dfRefSeq$taxa)

# Breakdown of a taxon by geographical region, please uncomment one set of these
# commands (ex: North American Region) to breakdown by a specific geographical 
# region. Breaking down by region will require 4 runthroughs after section 5 to 
# complete a taxon.

# ***For the large insect orders including Hymenoptera, Diptera and Lepidoptera,
# running by region will most likely be required due to their large size.***

# North American Region

# northAmericaFilter <- which(dfAllSeq$medianLatMap > 7 & 
#   dfAllSeq$medianLon > (-170) & dfAllSeq$medianLon < (-15) & 
#   dfAllSeq$medianLatMap < 90)
# dfAllSeq <- dfAllSeq[northAmericaFilter, ]

# South American Region

# southAmericaFilter <- which(dfAllSeq$medianLatMap < 7 & 
#   dfAllSeq$medianLon > (-135) & dfAllSeq$medianLon < (-30) & 
#   dfAllSeq$medianLatMap > (-90))
# dfAllSeq <- dfAllSeq[southAmericaFilter, ]

# Australasian Region

# australasianFilter1 <- which(dfAllSeq$medianLatMap > (-90) & 
#   dfAllSeq$medianLon < (-135) & dfAllSeq$medianLon > (-180) & 
#   dfAllSeq$medianLatMap < 0)
# australasianFilter2 <- which(dfAllSeq$medianLatMap > (-90) & 
#   dfAllSeq$medianLon < (180) & dfAllSeq$medianLon > (90) &
#   dfAllSeq$medianLatMap < 0)
# australasianFilter3 <- append(australasianFilter1, australasianFilter2)
# dfAllSeq <- dfAllSeq[australasianFilter3, ]

# Eurasia/African Region

# northAmericaFilter <- which(dfAllSeq$medianLatMap > 7 & 
#   dfAllSeq$medianLon > (-170) & dfAllSeq$medianLon < (-15) & 
#   dfAllSeq$medianLatMap < 90)
# southAmericaFilter <- which(dfAllSeq$medianLatMap < 7 & 
#   dfAllSeq$medianLon > (-135) & dfAllSeq$medianLon < (-30) & 
#   dfAllSeq$medianLatMap > (-90))
# australasianFilter1 <- which(dfAllSeq$medianLatMap > (-90) & 
#   dfAllSeq$medianLon < (-135) & dfAllSeq$medianLon > (-180) & 
#   dfAllSeq$medianLatMap<0)
# australasianFilter2 <- which(dfAllSeq$medianLatMap > (-90) & 
#   dfAllSeq$medianLon < (180) & dfAllSeq$medianLon > (90) &
#   dfAllSeq$medianLatMap < 0 )
# australasianFilter3 <- append(australasianFilter1, australasianFilter2)
# eurasiaAfricaFilter <- append(australasianFilter3, northAmericaFilter)
# eurasiaAfricaFilter2 <- append(southAmericaFilter, eurasiaAfricaFilter)
# dfAllSeq <- dfAllSeq[-eurasiaAfricaFilter2, ]

# Break dfAllSeq down into the various classes.
taxaListComplete <- lapply(unique(dfAllSeq$class_taxID), 
                           function(x) dfAllSeq[dfAllSeq$class_taxID == x,])

# Revise dfRefSeq dataframe to reflect this.
classList <- foreach(i=1:nrow(dfRefSeq)) %do% 
  unique(taxaListComplete[[i]]$class_name)
classList <- unlist(classList)
# This command will ensure the reference sequence dataframe is in the same order
# as the dfAllSeq dataframe.
dfRefSeq <- dfRefSeq[match(classList, dfRefSeq$taxa), ]

###################
# Section 6: Preliminary Alignment with Muscle and Pairwise Distance Matrix with
# the TN93 Model to Identify and Remove Divergent Sequences

# In this section we want to eliminate sequences that are too divergent in 
# pairwise distance by doing an initial alignment and distance matrix.
# Highly divergent sequences can cause major gaps in the alignment and 
# inaccurate pairwise distances between sequences.
# As well, the analysis further below compares latitudinally separated BINs that
# are relatively closely related.
# Using the Muscle alignment algorithm and the TN93 model.
# This ensures the final alignment will only be using sequences that are closely
# related to some others in pairwise distance.

# Extraction of Sequences from taxaListComplete for the alignment.
classSequences <- foreach(i=1:nrow(dfRefSeq)) %do% 
  taxaListComplete[[i]]$nucleotides
# Extraction of BINs from taxalist for the alignment.
classBin <- foreach(i=1:nrow(dfRefSeq)) %do% taxaListComplete[[i]]$bin_uri
# Conversion to DNAStringSet format.
dnaStringSet2 <- 
  foreach(i=1:nrow(dfRefSeq)) %do% DNAStringSet(classSequences[[i]])

# Naming the DNAStringSet with the appropriate bin_uri's.
for (i in seq(from = 1, to = nrow(dfRefSeq), by = 1)) {
  names(dnaStringSet2[[i]]) <- classBin[[i]]
}

# Multiple sequence alignment using muscle package on the dnaStringSet2 list for
# each class.
# Refer here for details on package: 
# https://www.bioconductor.org/packages/release/bioc/html/muscle.html
# This could take several minutes to hours depending on the taxa.
# Using default settings of package, the muscle command can detect DNA is being 
# used and align accordingly.
# The settings can be modified, if needed, depending upon the size of the 
# data set to be aligned and the patterns of sequence divergence in given data 
# set.
# ***For larger taxa, it is highly recommended to use a maxiter 
# setting of 2 to reduce alignment times.***
alignment2 <- foreach(i=1:length(classBin)) %do% 
  muscle(dnaStringSet2[[i]], maxiters = 2, diags = TRUE, gapopen = -3000)

# To check each preliminary alignment (each class) and output the preliminary 
# alignments to FASTA format, uncomment and run these three commands:

# classFileNames <- foreach(i=1:nrow(dfRefSeq)) %do% 
#   paste("alignmentPrelim",dfRefSeq$taxa[i],".fas",sep="")

# alignmentPrelimFasta <- foreach(i=1:nrow(dfRefSeq)) %do% 
#   DNAStringSet(alignment2[[i]])

# foreach(i=1:nrow(dfRefSeq)) %do% 
#   writeXStringSet(alignmentPrelimFasta[[i]], 
#   file=classFileNames[[i]], format = "fasta", width = 1500)

# ***FASTA files will show up in working directory of R and be named according 
# to the appropriate class being aligned.
# For instance the class Polychaeta would show up as the file 
# alignmentPrelimPolychaeta.fas.***

# Conversion to DNAbin format before using pairwise distance matrix.
dnaBin <- foreach(i=1:length(taxaListComplete)) %do% as.DNAbin(alignment2[[i]])

# Details on each model can be found here: 
# https://cran.rstudio.com/web/packages/ape/index.html
# We are using the TN93 model for our data, but the model selection can be 
# altered.
matrixGeneticDistance <- foreach(i=1:length(taxaListComplete)) %do% 
  dist.dna(dnaBin[[i]], model = "TN93", 
           as.matrix = TRUE, pairwise.deletion = TRUE)

# Conversion to dataframe format. 
geneticDistanceMatrixList <- foreach(i=1:length(taxaListComplete)) %do% 
  as.data.frame(matrixGeneticDistance[i])

# Putting it into a stack (each column concatenated into one long column of 
# indexes and values)so it can be easily subsetted.
geneticDistanceStackList <- foreach(i=1:length(taxaListComplete)) %do% 
  stack(geneticDistanceMatrixList[[i]])
# Identification and removal of divergent sequences, a divergent sequence being 
# one that is greater than 0.15 in pairwise distance to all other sequences.

# Finding all sequences between 0.15 and 0 distance.
# This command will remove divergent sequences which do not have pairwise 
# distances between 0.15 and 0.
divergentSequencesRemove <- foreach(i=1:length(taxaListComplete)) %do% 
  which(geneticDistanceMatrixList[[i]] <= 0.15 & 
        geneticDistanceMatrixList[[i]] > 0)
# Subsetting the distance stack by BINs meeting this criteria.
divergentSequencesRemove <- foreach(i=1:length(taxaListComplete)) %do% 
  geneticDistanceStackList[[i]][c(divergentSequencesRemove[[i]]), ]
# Grabbing all BINs from this stack that meet this criteria and converting to 
# character type.
divergentSequencesRemove <- foreach(i=1:length(taxaListComplete)) %do% 
  as.character(unique(divergentSequencesRemove[[i]]$ind))

###################
# Section 7: Finalized Alignment with Addition of Reference Sequence

# In this section we can perform the final alignment after divergent sequences 
# have been removed and the addition of a reference sequence.

# Subset taxalist by BINs meeting 0.15 criteria so no divergent sequences are 
# present.
taxaListComplete <- foreach(i=1:length(taxaListComplete)) %do% 
  subset(taxaListComplete[[i]], taxaListComplete[[i]]$bin_uri %in% 
         divergentSequencesRemove[[i]])

# Extract sequences and bin_uri from each class.
classBin <- foreach(i=1:nrow(dfRefSeq)) %do% taxaListComplete[[i]]$bin_uri
classSequences <- 
  foreach(i=1:nrow(dfRefSeq)) %do% taxaListComplete[[i]]$nucleotides
classSequencesNames <- classBin

# Take our reference sequences.
alignmentRef <- as.character(dfRefSeq$nucleotides)
dfRefSeq$reference <- "reference"
# Name our reference as reference for each class so it can be identified as such
# in the alignment.
alignmentRefNames <- dfRefSeq$reference

# Merge our reference sequences with each of our class sequences.
alignmentSequencesPlusRef <- foreach(i=1:nrow(dfRefSeq)) %do% 
  append(classSequences[[i]], alignmentRef[[i]])
# Merge the names together.
alignmentNames <- foreach(i=1:nrow(dfRefSeq)) %do% 
  append(classSequencesNames[[i]], alignmentRefNames[[i]])

# Converting all sequences in dfAllSeq plus reference to DNAStringSet format, 
# the format required for the alignment.
dnaStringSet3 <- foreach(i=1:nrow(dfRefSeq)) %do% 
  DNAStringSet(alignmentSequencesPlusRef[[i]])

# Name the DNAStringSet List with the appropriate BIN uri's.
for (i in seq(from=1, to=nrow(dfRefSeq), by = 1)){
  names(dnaStringSet3[[i]]) <- alignmentNames[[i]]
}

# Multiple sequence alignment using muscle package on each element of the 
# dnaStringSet2 list for each class.
# Refer here for details on package: 
# https://www.bioconductor.org/packages/release/bioc/html/muscle.html
# This could take several minutes to hours depending on the taxa. 
# Using default parameters, muscle will detect DNA is being used and will 
# align accordingly.
# ***For larger taxa, it is highly recommended to use a maxiter 
# setting of 2 to reduce alignment times.***
alignmentFinal <- foreach(i=1:length(classBin)) %do% 
  muscle(dnaStringSet3[[i]], maxiters = 2, diags = TRUE, gapopen = -3000)

# To check each final alignment (each class) and output the final alignments to 
# FASTA format, uncomment and run these three commands:

# classFileNames2 <- foreach(i=1:nrow(dfRefSeq)) %do% 
#   paste("alignmentFinal", dfRefSeq$taxa[i], ".fas", sep="")

# alignmentFinalFasta <- foreach(i=1:nrow(dfRefSeq)) %do% 
#   DNAStringSet(alignmentFinal[[i]])

# foreach(i=1:nrow(dfRefSeq)) %do% 
#   writeXStringSet(alignmentFinalFasta[[i]], file=classFileNames2[[i]], 
#                   format="fasta", width=1500)

# ***FASTA files will show up in working directory of R and be named according 
# to the appropriate class being aligned. For instance the class Polychaeta 
# would show up as the file alignmentFinalPolychaeta.fas.***

# Removal of unecessary large variables to minimize workspace size and 
# free up working memory
rm(geneticDistanceStackList)
rm(geneticDistanceMatrixList)
rm(matrixGeneticDistance)
rm(geneticDistanceCentroid)
rm(dfInitial)
rm(end_N_gap)
rm(start_N_gap)
rm(internal_N_gap)
rm(binList)
gc()

#***Save workspace before proceeding further*** 

###################
# Section 8: Sequence Trimming of Finalized Alignment According to the Reference 
# Sequence

# In this section, the final alignment generated is now trimmed according to the 
# position of the reference in the alignment. This step ensures that a 
# consistent genetic region is used for the subsequent distance calculations.

# For trimming of the sequences we have to determine where in the alignment the 
# reference sequence is and determine its start and stop positions relative to 
# the other sequences. We can then use these positions to trim the rest of the 
# sequences in the alignment.
refSeqPos <- foreach(i=1:nrow(dfRefSeq)) %do% 
  which(alignmentFinal[[i]]@unmasked@ranges@NAMES == "reference")
refSeqPos <- foreach(i=1:nrow(dfRefSeq)) %do% 
  alignmentFinal[[i]]@unmasked[refSeqPos[[i]]]

# Finding start position by searching for the first nucleotide position of the 
# reference sequence.
refSeqPosStart <- foreach(i=1:nrow(dfRefSeq)) %do% 
  regexpr("[ACTG]", refSeqPos[[i]])
refSeqPosStart <- as.numeric(refSeqPosStart)

# Finding last nucleotide position of the reference sequence.
refSeqPosEnd <- foreach(i=1:nrow(dfRefSeq)) %do% 
  (nchar(dfRefSeq$nucleotides[i]) + refSeqPosStart[i])
refSeqPosEnd <- as.numeric(refSeqPosEnd)

# Then we can substr the alignment by these positions to effectively trim the 
# alignment.
alignmentFinalTrim <- foreach(i=1:nrow(dfRefSeq)) %do% 
  substr(alignmentFinal[[i]], refSeqPosStart[i] + 1, refSeqPosEnd[i])

# To check each final trimmed alignment (each class) and output the final 
# trimmed alignments to FASTA format, uncomment and run these three commands:

# classFileNames3 <- foreach(i=1:nrow(dfRefSeq)) %do% 
#   paste("alignmentFinalTrim", dfRefSeq$taxa[i], ".fas", sep="")

# alignmentFinalTrimFasta <- foreach(i=1:nrow(dfRefSeq)) %do% 
#   DNAStringSet(alignmentFinalTrim[[i]])

# foreach(i=1:nrow(dfRefSeq)) %do% 
#   writeXStringSet(alignmentFinalTrimFasta[[i]], file=classFileNames3[[i]], 
#   format="fasta", width=658)

# ***FASTA files will show up in working directory of R and be named according 
# to the appropriate class being aligned. For instance the class Polychaeta 
# would show up as the file alignmentFinalTrimPolychaeta.fas.***

# Again, convert to dnaStringSet format.
dnaStringSet4 <- foreach(i=1:nrow(dfRefSeq)) %do% 
  DNAStringSet(alignmentFinalTrim[[i]])

# Establish where the reference sequence is in each alignment for removal of 
# the reference from further analysis.
refSeqRemove <- foreach(i=1:nrow(dfRefSeq)) %do% 
  which(dnaStringSet4[[i]]@ranges@NAMES == "reference")

# We also have to check dnaStringSet4 for alignments that only contain a 
# reference sequence and remove these.
alignmentCheck <- foreach(i=1:nrow(dfRefSeq)) %do% 
  which(length(dnaStringSet4[[i]]) == 1)
alignmentCheck <- which(alignmentCheck > 0)
# Subset dnaStringSet4, taxalistcomplete, alignmentFinalTrim and refSeqRemove by 
# the alignment check.
if (length(alignmentCheck > 0)) {
  dnaStringSet4 <- dnaStringSet4[-alignmentCheck]
  taxaListComplete <- taxaListComplete[-alignmentCheck]
  alignmentFinalTrim <- alignmentFinalTrim[-alignmentCheck]
  refSeqRemove <- refSeqRemove[-alignmentCheck]
}

# Removal of references from both dnaStringSet4 and alignmentFinalTrim.
dnaStringSet4 <- foreach(i=1:length(taxaListComplete)) %do% 
  subset(dnaStringSet4[[i]][-refSeqRemove[[i]]])
alignmentFinalTrim <- foreach(i=1:length(taxaListComplete)) %do% 
  alignmentFinalTrim[[i]][-refSeqRemove[[i]]]

###################
# Section 9: Pairwise Distance Determination with the TN93 Model of Finalized 
# Alignment

# In this section we do the final pairwise distance determination step using the 
# TN93 model, before initial pairings can be generated.

# Conversion to DNAbin format before using genetic distance matrix.
dnaBin2 <- foreach(i=1:length(taxaListComplete)) %do% 
  as.DNAbin(dnaStringSet4[[i]])

# Once again using the TN93 model for pairwise distance computation, this time 
# from the 3rd alignment after divergent sequences have been removed and after 
# trimming according to the reference sequences.
matrixGeneticDistance2 <- foreach(i=1:length(taxaListComplete)) %do% 
  dist.dna(dnaBin2[[i]], model = "TN93", as.matrix = TRUE, 
           pairwise.deletion = TRUE)

# Convert to dataframe format.
geneticDistanceMatrixList2 <- foreach(i=1:length(taxaListComplete)) %do% 
  as.data.frame(matrixGeneticDistance2[i])

# Putting it into a stack (each column concatenated into one long column of 
# indexes and values) so it can be easily subsetted.
geneticDistanceStackList2 <- foreach(i=1:length(taxaListComplete)) %do% 
  stack(geneticDistanceMatrixList2[[i]])

###################
# Section 10: Finding Appropriate Pairings According to Genetic Distance Criteria 
# and Latitude Criteria

# In this section, we establish preliminary pairings of BINs showing up to a 
# maximum of 0.15 genetic divergence and 20 degrees latitude difference
# and divide preliminary pairings into distinct lineages.

# These values can easily be edited to add more or less stringency to the 
# matches. Will produce lists with indexes (BINs, row and column names) of each
# match according to the maximum genetic distance criterion of 0.15.
pairingResultCandidates <- foreach(i=1:length(taxaListComplete)) %do% 
  which(geneticDistanceMatrixList2[[i]] <= 0.15)

# Generate a matrix of latitude differences between bins.
matrixLatitudeDistance <- foreach(i=1:length(taxaListComplete)) %do% 
  as.matrix( dist(taxaListComplete[[i]]$medianLatAbs) )
# Establish which bins meet a min lat difference of 20 degrees.
latitudeResultCandidates <- foreach(i=1:length(taxaListComplete)) %do% 
  which(matrixLatitudeDistance[[i]] > 20)
# Find pairs which have both a min lat difference of 20 degrees and also have a 
# divergence under 0.15.
pairingResultCandidates <- foreach(i=1:length(taxaListComplete)) %do%
  intersect(pairingResultCandidates[[i]], latitudeResultCandidates[[i]])

# Next we can delete the classes that don't have any candidate pairings.
pairingResultCheck <- foreach(i=1:length(taxaListComplete)) %do% 
  which(length(pairingResultCandidates[[i]]) == 0)
pairingResultCheck <- which(pairingResultCheck > 0)
# Subset taxalistcomplete, geneticDistance, pairingResultCandidate, 
# dnaStringSet lists according to pairingResultCheck.
# so that we aren't keeping classes without pairings in these lists.
if (length(pairingResultCheck > 0)) {
  taxaListComplete <- taxaListComplete[-pairingResultCheck]
  geneticDistanceStackList2 <- geneticDistanceStackList2[-pairingResultCheck]
  pairingResultCandidates <- pairingResultCandidates[-pairingResultCheck]
  dnaStringSet4 <- dnaStringSet4[-pairingResultCheck]
}

# Use these pairing results and reference against the genetic distance stack 
# list to subset it.
pairingResultCandidates2 <- foreach(i=1:length(taxaListComplete)) %do% 
  geneticDistanceStackList2[[i]][c(pairingResultCandidates[[i]]),]
# This command will grab the precise row and column numbers of each pairing so 
# they can be used to determine a unique identifier for each pairing.
pairingResultCandidates3 <- foreach(i=1:length(taxaListComplete)) %do%
  data.table(which(geneticDistanceMatrixList2[[i]] <= 0.15 & 
                   matrixLatitudeDistance[[i]] > 20, arr.ind = TRUE))

# Merge pairingResultCandidate list 2 and 3 together and combine into the 
# pairing results dataframe.
dfPairingResultsL1L2 <- do.call(rbind, Map(data.frame, pairingResultCandidates2,
                                           pairingResultCandidates3))

# Combine AllSeq with pairingResultCandidates to get all relevant taxonomic 
# and latitudinal data.
dfPairingResultsL1L2 <- suppressWarnings(merge(dfPairingResultsL1L2, dfAllSeq, 
                                               by.x = "ind", by.y = "bin_uri"))

# Creating a key that can be used to uniquely identify each pairing based on its 
# precise position in its respective pairwise distance matrix from the row and 
# column number. This key is generated using a pairing function called
# called the Cantor pairing function which will generate a unique value
# value for each pairing.

# First adding min and max columns to the dfPairingResultsL1L2 dataframe
# so that this key can be calculated. Min refers to the min of value from
# each row/col and max refers to the max value from each row/col
dfPairingResultsL1L2 <- transform(dfPairingResultsL1L2, max = pmax(row, col)) 
dfPairingResultsL1L2 <- transform(dfPairingResultsL1L2, min = pmin(row, col))

# Adding the class_taxID to the max and min value of each lineage so that
# no pairing duplicates can be present across classes
dfPairingResultsL1L2$max <- dfPairingResultsL1L2$max + 
  dfPairingResultsL1L2$class_taxID
dfPairingResultsL1L2$min <- dfPairingResultsL1L2$min + 
  dfPairingResultsL1L2$class_taxID

# Inputting min and max into the Cantor pairing function so that a unique
# key can be created for each pairing.
dfPairingResultsL1L2$pairingKey <- 0.5 *
  (dfPairingResultsL1L2$max + dfPairingResultsL1L2$min) *
  (dfPairingResultsL1L2$max + dfPairingResultsL1L2$min + 1) +
   dfPairingResultsL1L2$min

# order by pairingKey so all pairings are ordered correctly.
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$pairingKey), ]

# Changing significant digits back.
options(digits =5)
# Get rid of zero ingroupdistance entries.
zeroDistance <- which(dfPairingResultsL1L2$values == 0)
if (length(zeroDistance) > 0) {
  dfPairingResultsL1L2 <- dfPairingResultsL1L2[-zeroDistance, ]
}

# Then we can multiply these ingroupdistance values by 1.3 to determine the
# minimum outgroup distance from the pairings. 
# We put this in another column in the pairing results dataframe.
# This minimum outgroup distance could also be a user adjustable parameter and 
# can be easily modified.
dfPairingResultsL1L2$inGroupDistx1.3 <- dfPairingResultsL1L2$values * 1.3

# Also grouping dataframe every 2 rows to reflect each unique pairing/matching, 
# pairing column will give a number value for each pairing ordered.
dfPairingResultsL1L2$inGroupPairing <- 
  rep(1:(nrow(dfPairingResultsL1L2) / 2), each = 2)

# Reorganizing and renaming some columns in pairing result candidate dataframe 
# to make more easily readable.
colnames(dfPairingResultsL1L2)[1] <- "bin_uri"
dfPairingResultsL1L2 <- 
  (dfPairingResultsL1L2[,c("inGroupPairing","record_id","bin_uri","values",
                           "inGroupDistx1.3","medianLatAbs","medianLatMap",
                           "latMin","latMax","binSize","phylum_taxID",
                           "phylum_name","class_taxID","class_name",
                           "order_taxID","order_name","family_taxID",
                           "family_name","subfamily_taxID","subfamily_name",
                           "genus_taxID","genus_name","species_taxID",
                           "species_name","nucleotides","ind.1","medianLon",
                           "pairingKey")])
colnames(dfPairingResultsL1L2)[4] <- "inGroupDist"
colnames(dfPairingResultsL1L2)[26] <- "indexNo"
dfPairingResultsL1L2$index <- seq.int(nrow(dfPairingResultsL1L2))
                              
# Order by ingrouppairing between pairings.
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupPairing),]

# Dividing Pairing Candidates into two lineages:

# First lineage, this is done by picking even-numbered rows of pairing results 
# dataframe: 2,4,6...
dfPairingResultsL1 <- 
  dfPairingResultsL1L2[dfPairingResultsL1L2$index %% 2 == 0, ]
# Second lineage, this is done with odd-numbered rows.
dfPairingResultsL2 <- 
  dfPairingResultsL1L2[dfPairingResultsL1L2$index %%2 > 0, ]

###################
# Section 11: Ensuring Every Pairing is a Unique Set of BINs and Addition of 
# Trimmed Sequences to Pairing Results

# In this section we want to eliminate duplicate BINs in the pairing results.
# If a BIN is found in multiple pairings that meet the latitude criterion, 
# then we will only retain the pairing with the smallest ingroup distance.
# Pairings are also numbered based on ordering from lowest ingroup distance to 
# highest.

# Once again, order by ingroup distance.
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupDist), ]
# Take the first instance of each BIN, the first instance of it being the 
# pairing in which it exhibits the lowest ingroup distance.
dfPairingResultsL1L2 <- 
  by(dfPairingResultsL1L2, dfPairingResultsL1L2["bin_uri"], head, n=1)
dfPairingResultsL1L2 <- Reduce(rbind, dfPairingResultsL1L2)
dfPairingResultsL1 <- 
  subset(dfPairingResultsL1L2, duplicated(dfPairingResultsL1L2$inGroupPairing))
#Subsetting against L1L2 to get finalized pairings in L1L2.
dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, 
                               dfPairingResultsL1L2$inGroupPairing %in% 
                               dfPairingResultsL1$inGroupPairing)
# Now rebuilding L1 and L2 based on a subsetted L1L2.
dfPairingResultsL1 <- 
  dfPairingResultsL1L2[dfPairingResultsL1L2$index %% 2 == 0, ]
dfPairingResultsL2 <- 
  dfPairingResultsL1L2[dfPairingResultsL1L2$index %% 2 > 0, ]

# Merge aligned and trimmed nucleotide sequences with pairing candidates, 
# replacing raw sequences from intial tsv.
alignmentFinalTrimUnlist <- unlist(alignmentFinalTrim)
# New dataframe with trimmed and aligned sequences.
dfTrimmedSeq <- data.frame(alignmentFinalTrimUnlist)
dfTrimmedSeq$bin_uri <- names(alignmentFinalTrimUnlist)
colnames(dfTrimmedSeq)[1] <- "trimmedNucleotides"

# Merge sequences to dfPairingResultsL1L2, L1 and L2 dataframes.
dfPairingResultsL1L2 <- merge(dfPairingResultsL1L2, dfTrimmedSeq, 
                              by.x = "bin_uri", by.y = "bin_uri")
dfPairingResultsL1 <- merge(dfPairingResultsL1, dfTrimmedSeq, 
                            by.x = "bin_uri", by.y = "bin_uri")
dfPairingResultsL2 <- merge(dfPairingResultsL2, dfTrimmedSeq, 
                            by.x = "bin_uri", by.y = "bin_uri")

# Again reorganization of pairingresults to have trimmedNucleotides added 
# instead of nucleotides.
dfPairingResultsL1L2 <- 
  (dfPairingResultsL1L2[,c("inGroupPairing","record_id","bin_uri","inGroupDist",
                           "inGroupDistx1.3","medianLatAbs","medianLatMap",
                           "latMin","latMax","binSize","phylum_taxID",
                           "phylum_name","class_taxID","class_name",
                           "order_taxID","order_name","family_taxID",
                           "family_name","subfamily_taxID","subfamily_name",
                           "genus_taxID","genus_name","species_taxID",
                           "species_name","trimmedNucleotides","indexNo",
                           "medianLon","index","pairingKey")])
dfPairingResultsL1 <- 
  (dfPairingResultsL1[,c("inGroupPairing","record_id","bin_uri","inGroupDist",
                         "inGroupDistx1.3","medianLatAbs","medianLatMap",
                         "latMin","latMax","binSize","phylum_taxID",
                         "phylum_name","class_taxID","class_name","order_taxID",
                         "order_name","family_taxID","family_name",
                         "subfamily_taxID","subfamily_name","genus_taxID",
                         "genus_name","species_taxID","species_name",
                         "trimmedNucleotides","indexNo","medianLon","index",
                         "pairingKey")])
dfPairingResultsL2 <- 
  (dfPairingResultsL2[,c("inGroupPairing","record_id","bin_uri","inGroupDist",
                         "inGroupDistx1.3","medianLatAbs","medianLatMap",
                         "latMin","latMax","binSize","phylum_taxID",
                         "phylum_name","class_taxID","class_name","order_taxID",
                         "order_name","family_taxID","family_name",
                         "subfamily_taxID","subfamily_name","genus_taxID",
                         "genus_name","species_taxID","species_name",
                         "trimmedNucleotides","indexNo","medianLon","index",
                         "pairingKey")])

# Order and renumber pairings starting at 1 according to pairingKey.
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$pairingKey), ]
dfPairingResultsL1L2$inGroupPairing <- 
  rep(1:(nrow(dfPairingResultsL1L2) / 2), each = 2)
dfPairingResultsL1 <- 
  dfPairingResultsL1[order(dfPairingResultsL1$pairingKey), ]
dfPairingResultsL1$inGroupPairing <- 1:nrow(dfPairingResultsL1)
dfPairingResultsL2 <- 
  dfPairingResultsL2[order(dfPairingResultsL2$pairingKey), ]
dfPairingResultsL2$inGroupPairing <- 1:nrow(dfPairingResultsL2)

# Removal of variables to free up working memory and reduce workspace size.
rm(matrixLatitudeDistance)
rm(matrixGeneticDistance2)
rm(pairingResultCandidates)
rm(pairingResultCandidates3)
rm(pairingResultCandidates2)
rm(geneticDistanceMatrixList2)
rm(latitudeResultCandidates)
gc()

# ***Save workspace before proceeding***

###################
# Section 12: Eliminating Pairings based on Overlapping Latitudinal Range
# In this section we eliminate pairings which have latitudinal range overlap of 
# greater than 25%.
# This section may take a few mins to process.

# Before eliminating based on overlapping latitudinal range, another check of 
# latitude difference is performed.
# This is to address an issue where some pairings may still have latitudes not 
# meeting the 20 degrees criteria.

# Determining latitude difference between pairings:
latitudeDiffCheck <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  abs(dfPairingResultsL1$medianLatAbs[i] - dfPairingResultsL2$medianLatAbs[i])

# Ensuring all pairings meet the latitude difference of at least 20 degrees 
# (this choice can be modified by the user).
latitudeDiffCheck <- which(latitudeDiffCheck<20)
if (length(latitudeDiffCheck) > 0) {
  # Subsetting Lineage 1 by this latitude check.
  dfPairingResultsL1 <- dfPairingResultsL1[-latitudeDiffCheck, ]
  # Now subsetting Lineage 2.
  dfPairingResultsL2 <- dfPairingResultsL2[-latitudeDiffCheck, ]
  # Now subsetting for dataframe with both lineages.
  dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, 
                                 dfPairingResultsL1L2$inGroupPairing %in% 
                                 dfPairingResultsL1$inGroupPairing)
}

# Next we can work on establishing latitudinal ranges for each pairing.
# If the two lineages of a pairing have overlapping latitude regions of greater 
# than 25% then we would not consider that pairing as a viable pairing.
# This is because we want each lineage of a pairing to meet an appropriate 
# difference in latitude.

# We can define the overlap range threshold as 25% of the latitude range of L1.
# If an overlap is greater than this value we would discard this pairing.
# Of course this value could be easily modified to add more or less stringency 
# to the script.
rangeThresholdL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  (abs(dfPairingResultsL1$latMax[i] - dfPairingResultsL1$latMin[i]) * 0.25) 

# Same process for L2, as we apply the overlap criterion to each sister lineage 
# in a pairing.
rangeThresholdL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  (abs(dfPairingResultsL2$latMax[i] - dfPairingResultsL2$latMin[i]) * 0.25) 

# Define our latitude ranges for each lineage of a pairing.
rangeL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  range(dfPairingResultsL1$latMax[i], dfPairingResultsL1$latMin[i])
rangeL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  range(dfPairingResultsL2$latMax[i], dfPairingResultsL2$latMin[i])

# Then we can determine the overlap region between them using the Overlap 
# function from the Desctools package.
# Overlap will return an absolute value so we don't have to worry about 
# negatives for overlap values.
overlapValueL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  Overlap(rangeL1[[i]], 
          rangeL2[[i]])
overlapValueL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  Overlap(rangeL2[[i]], 
          rangeL1[[i]])

# Then if there is a range overlap between two lineages in a pairing, we can 
# determine if this overlap is actually larger than the 25% value of 
# rangeThreshold for each individual pairing.
rangeOverlapCheckL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(rangeThresholdL1[[i]]<overlapValueL1[[i]])

# Same process for L2.
rangeOverlapCheckL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(rangeThresholdL2[[i]]<overlapValueL2[[i]])

# Then overlap values exceeding that 25% value we set should be returned as an 
# integer of 1. If an overlap does not exceed this value, then it will return an 
# integer of 0. If there is a value that meets these criteria, its associated 
# pairing will be removed from the pairing results dataframe. This will name 
# each element of rangeOverlapCheck with the inGroupPairing number to identify 
# the pairing we need to eliminate.
names(rangeOverlapCheckL1) <- paste0(dfPairingResultsL1$inGroupPairing)
names(rangeOverlapCheckL2) <- paste0(dfPairingResultsL2$inGroupPairing)
# Identify which pairing in rangeOverlapCheck is greater than 0; do this with 
# respect to both lineages.
overlapIndL1 <- which(rangeOverlapCheckL1>0)
overlapIndL2 <- which(rangeOverlapCheckL2>0)
overlapIndTotal <- append(overlapIndL1,overlapIndL2)
overlapIndTotal <- unique(overlapIndTotal)
# If overlapIndTotal is not empty:
if (length(overlapIndTotal) > 0) {
  # Eliminate based on that pairing number(s) for each pairing result dataframe
  dfPairingResultsL1 <- dfPairingResultsL1[-overlapIndTotal, ]
  dfPairingResultsL2 <- dfPairingResultsL2[-overlapIndTotal, ]
  dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, 
                                 dfPairingResultsL1L2$inGroupPairing %in% 
                                 dfPairingResultsL1$inGroupPairing)
}

# Reorder and renumber pairings again
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$pairingKey),]
dfPairingResultsL1L2$inGroupPairing <- 
  rep(1:(nrow(dfPairingResultsL1L2) / 2), each = 2)
dfPairingResultsL1 <- 
  dfPairingResultsL1[order(dfPairingResultsL1$pairingKey), ]
dfPairingResultsL1$inGroupPairing <- 1:nrow(dfPairingResultsL1)
dfPairingResultsL2 <- 
  dfPairingResultsL2[order(dfPairingResultsL2$pairingKey), ]
dfPairingResultsL2$inGroupPairing <- 1:nrow(dfPairingResultsL2)

###################
# Section 13: Outgroup Determination for Each Pairing

# In this section, we use the 1.3x divergence criterion to establish a distance 
# range for outgroups. Outgroup BINs are then determined based on this 
# divergence criterion.

# First we can take all genetic distance values from each distance matrix and 
# concatenate them into one stack of values to be easily subsetted from.
dfGeneticDistanceStack <- do.call(rbind, geneticDistanceStackList2)

# Then we can search for outgroupings relative to lineage 1.
# We can use the bin_uri's to subset the dfGeneticDistanceStackList according to 
# bin_uri's represented in lineage 1. This will be called dfOutGroupL1.
# This will limit our outgroup distances to those associated with lineage1.
dfOutGroupL1 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% 
                         dfPairingResultsL1$bin_uri)

# Order outgroups by Lineage1.
dfOutGroupL1 <- dfOutGroupL1[order(match(dfOutGroupL1[, 2],
                                         dfPairingResultsL1[, 3])), ]

# Find which BINs in the outgroup genetic distances match the BINs in the 
# pairing results.
outGroupCandidatesL1a <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfOutGroupL1$ind == dfPairingResultsL1$bin_uri[i])

# Find which outgroups meet the 1.3x criterion.
outGroupCandidatesL1b <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfOutGroupL1$values >= dfPairingResultsL1$inGroupDistx1.3[i])

# Intersection of the two using mapply to find the correct outgroupings for each
# pairing.
outGroupCandidatesL1c <- 
  mapply(intersect, outGroupCandidatesL1a, outGroupCandidatesL1b)
# Unlist to make into one vector.
outGroupCandidatesL1c <- unlist(outGroupCandidatesL1c)
# Adding an rownum column to dfBestOutGroup. This represents the second index of 
# each outgroup candidate.
dfOutGroupL1$rownum <- seq.int(nrow(dfOutGroupL1))
# Then we can subset to rownum column based on the outgroup candidates.
dfOutGroupL1 <- dfOutGroupL1[dfOutGroupL1$rownum %in% 
                               outGroupCandidatesL1c, ]
# Same process for lineage 2.
dfOutGroupL2 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% 
                       dfPairingResultsL2$bin_uri)
dfOutGroupL2 <- dfOutGroupL2[order(match(dfOutGroupL2[, 2],
                                         dfPairingResultsL2[, 3])), ]

# Find which BINs in the outgroup genetic distances match the BINs in the 
# pairing results.
outGroupCandidatesL2a <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfOutGroupL2$ind == dfPairingResultsL2$bin_uri[i])

# Find which outgroups meet the 1.3x criterion.
outGroupCandidatesL2b <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfOutGroupL2$values >= dfPairingResultsL2$inGroupDistx1.3[i])

# Intersection of the two using mapply to find the correct outgroupings for 
# each pairing.
outGroupCandidatesL2c <- 
  mapply(intersect, outGroupCandidatesL2a, outGroupCandidatesL2b)
# Unlist to make into one vector.
outGroupCandidatesL2c <- unlist(outGroupCandidatesL2c)
# Adding a rownum column to dfBestOutGroup. This represents the second index of 
# each outgroup candidate.
dfOutGroupL2$rownum <- seq.int(nrow(dfOutGroupL2))
# Then we can subset to rownum column based on the outgroup candidates.
dfOutGroupL2 <- dfOutGroupL2[dfOutGroupL2$rownum %in% 
                               outGroupCandidatesL2c, ]

# This next command will determine which outgroup candidates are shared between 
# both lineages.
dfOutGroupL1 <- dfOutGroupL1[dfOutGroupL1$rownum %in% dfOutGroupL2$rownum, ]
dfOutGroupL2 <- dfOutGroupL2[dfOutGroupL2$rownum %in% dfOutGroupL1$rownum, ]

# Then we can determine averages between outgroup distances and find the 
# outgroup with the min average distance to ensure the outgroup is still 
# relatively close in distance to each lineage but at least at the 1.3x 
# divergence value determined for each pairing.
dfOutGroupL1$distAverage <- (dfOutGroupL1$values + dfOutGroupL2$values) / 2
dfOutGroupL2$distAverage <- (dfOutGroupL1$values + dfOutGroupL2$values) / 2
dfOutGroupL1 <- dfOutGroupL1[!duplicated(dfOutGroupL1$distAverage), ]
dfOutGroupL2 <- dfOutGroupL2[!duplicated(dfOutGroupL2$distAverage), ]
outGroupList <- lapply(unique(dfOutGroupL1$ind), 
                       function(x) dfOutGroupL1[dfOutGroupL1$ind == x,])
minOutGroupAverage <- sapply( outGroupList , function(x) min( x$distAverage ) )

# Subsetting each outgroup dataframe by our min outgroup averages.
dfOutGroupL1 <- subset(dfOutGroupL1, distAverage %in% minOutGroupAverage)
dfOutGroupL2 <- subset(dfOutGroupL2, distAverage %in% minOutGroupAverage)

# Next find the bin_uri's for each outgroup candidate by subsetting the 
# dfGeneticDistanceStack dataframe with values from the L1 and L2 dataframes.
# L1.
dfOutGroupL1 <- merge(dfOutGroupL1, dfGeneticDistanceStack, by.x = "values", 
                      by.y = "values", all.x = TRUE)
dfOutGroupL1 <- dfOutGroupL1[order(match(dfOutGroupL1[,2], 
                                         dfPairingResultsL1[,3])), ]
outGroupL1Check <- which(dfOutGroupL1$ind.x != dfOutGroupL1$ind.y)
dfOutGroupL1 <- dfOutGroupL1[outGroupL1Check, ]

# L2 (same process).
dfOutGroupL2 <- merge(dfOutGroupL2, dfGeneticDistanceStack, by.x = "values", 
                      by.y = "values", all.x = TRUE)
dfOutGroupL2 <- dfOutGroupL2[order(match(dfOutGroupL2[,2], 
                                         dfPairingResultsL2[,3])), ]
outGroupL2Check <- which(dfOutGroupL2$ind.x != dfOutGroupL2$ind.y)
dfOutGroupL2 <- dfOutGroupL2[outGroupL2Check, ]

# Removal of variables not used anymore to reduce R workspace size and free 
# memory before dfOutGroupMerge line
rm(outGroupCandidatesL2b)
rm(outGroupCandidatesL1b)
rm(geneticDistanceStackList2)
rm(dfGeneticDistanceStack)
rm(outGroupList)
gc()

#***Save workspace before proceeding further

# For faster merging of Outgroup dataframes, converting to data tables.
dfOutGroupL1 <- data.table(dfOutGroupL1)
dfOutGroupL2 <- data.table(dfOutGroupL2)
# Setting keys for merge.
setkey(dfOutGroupL1, ind.y)
setkey(dfOutGroupL2, ind.y)

# Merging of outgroup dataframes so final filtering of outgroups can be done.
dfOutGroupMerge <- merge(dfOutGroupL1, dfOutGroupL2, by.x = "ind.y", 
                         by.y = "ind.y", all.x = TRUE)
outGroupMergeCheck <- 
  which(dfOutGroupMerge$rownum.x == dfOutGroupMerge$rownum.y)
dfOutGroupMerge <- dfOutGroupMerge[outGroupMergeCheck,]
dfOutGroupMerge <- dfOutGroupMerge[!duplicated(dfOutGroupMerge$distAverage.x), ]

# Reorganization and merging of outgroup dataframes with dfAllSeq to get all 
# taxonomic information.
# L1.
dfOutGroupL1 <- (dfOutGroupMerge[,c("ind.y","values.x","ind.x.x")])
dfOutGroupL1 <- merge(dfOutGroupL1, dfAllSeq, by.x = "ind.y", by.y = "bin_uri")
# Then merge with dfTrimmedSeq to replace raw sequences with trimmed and 
# aligned ones.
dfOutGroupL1 <- merge(dfOutGroupL1, dfTrimmedSeq, by.x = "ind.y", 
                      by.y = "bin_uri")
dfOutGroupL1 <- subset(dfOutGroupL1, select = -c(nucleotides) )
# Renaming a few columns.
colnames(dfOutGroupL1)[1] <- "outGroupBin"
colnames(dfOutGroupL1)[2] <- "outGroupDistance"
colnames(dfOutGroupL1)[3] <- "associatedInGroupBin"

# Same process for L2.
dfOutGroupL2 <- (dfOutGroupMerge[,c("ind.y","values.y","ind.x.y")])
dfOutGroupL2 <- merge(dfOutGroupL2, dfAllSeq, by.x = "ind.y", by.y = "bin_uri")
dfOutGroupL2 <- merge(dfOutGroupL2, dfTrimmedSeq, by.x = "ind.y", 
                      by.y = "bin_uri")
dfOutGroupL2 <- subset(dfOutGroupL2, select = -c(nucleotides) )
# Renaming a few columns.
colnames(dfOutGroupL2)[1] <- "outGroupBin"
colnames(dfOutGroupL2)[2] <- "outGroupDistance"
colnames(dfOutGroupL2)[3] <- "associatedInGroupBin"

# Remove outGroupMerge and dfTrimmedSeq since these are redundant and 
# no longer needed.
rm(dfOutGroupMerge)
rm(dfTrimmedSeq)

# Some pairings may not have any viable outgroupings, we can filter those out.
noOutGroupCheck <- setdiff(dfPairingResultsL1$bin_uri, 
                           dfOutGroupL1$associatedInGroupBin)

# If there is at least one pairing without an outgroup, then we subset pairing 
# results dataframes by those pairings.
if (length(noOutGroupCheck) > 0) {
  noOutGroupCheck <- foreach(i=1:length(noOutGroupCheck)) %do% 
    which(dfPairingResultsL1$bin_uri == noOutGroupCheck[[i]])
  noOutGroupCheck <- unlist(noOutGroupCheck)
  # Subsetting Lineage 1 by this outgroup check.
  dfPairingResultsL1 <- dfPairingResultsL1[-noOutGroupCheck, ]
  # Now subsetting Lineage 2 by this outgroup check.
  dfPairingResultsL2 <- dfPairingResultsL2[-noOutGroupCheck, ]
  # Now subsetting for dataframe with both lineages.
  dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, 
                                 dfPairingResultsL1L2$inGroupPairing %in% 
                                 dfPairingResultsL1$inGroupPairing)
}

# Renumber pairings again.
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$pairingKey),]
dfPairingResultsL1L2$inGroupPairing <- 
  rep(1:(nrow(dfPairingResultsL1L2) / 2), each = 2)
dfPairingResultsL1 <- dfPairingResultsL1[order(dfPairingResultsL1$pairingKey),]
dfPairingResultsL1$inGroupPairing <- 1:nrow(dfPairingResultsL1)
dfPairingResultsL2 <- dfPairingResultsL2[order(dfPairingResultsL2$pairingKey),]
dfPairingResultsL2$inGroupPairing <- 1:nrow(dfPairingResultsL2)

# Removal of more variables
rm(dfOutGroupMerge)
rm(dfTrimmedSeq)
rm(outGroupCandidatesL1a)
rm(outGroupCandidatesL2a)
rm(outGroupCandidatesL1c)
rm(outGroupCandidatesL2c)
gc()

# ***Save workspace before proceeding further***

###################
# Section 14: Finalizing the Pairing Results Dataframes

# In this section the final pairing results are established with outgroups.
# Pairing result dataframes also undergo some reorganization.

# First order outgroup lineage dataframes again by BIN.
dfOutGroupL1 <- data.frame(dfOutGroupL1)
dfOutGroupL1 <- dfOutGroupL1[order(match(dfOutGroupL1[,3],
                                         dfPairingResultsL1[,3])), ]
dfOutGroupL2 <- data.frame(dfOutGroupL2)
dfOutGroupL2 <- dfOutGroupL2[order(match(dfOutGroupL2[,3],
                                         dfPairingResultsL2[,3])), ]

# Adding an additional column for the delta lat difference of each pairing.
latDelta <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  abs(dfPairingResultsL1$medianLatAbs[i] - dfPairingResultsL2$medianLatAbs[i])

# Finally merge outgroups to both L1, L2 and L1L2 dataframes.
# Using the data tables package for this for reliable merging.
# Making these dataframes data tables.
dfPairingResultsL1 <- data.table(dfPairingResultsL1)
dfOutGroupL1 <- data.table(dfOutGroupL1)
dfOutGroupL1$inGroupPairing <- seq.int(nrow(dfOutGroupL1))
# Setting keys for merge.
setkey(dfOutGroupL1, inGroupPairing)
setkey(dfPairingResultsL1, inGroupPairing)
dfPairingResultsL1 <- merge(dfPairingResultsL1, dfOutGroupL1, 
                            by.x = "inGroupPairing", by.y = "inGroupPairing")
dfPairingResultsL1 <- data.frame(dfPairingResultsL1)
dfPairingResultsL1$latDelta <- latDelta

# L2, same process.
dfPairingResultsL2 <- data.table(dfPairingResultsL2)
dfOutGroupL2 <- data.table(dfOutGroupL2)
dfOutGroupL2$inGroupPairing <- seq.int(nrow(dfOutGroupL2))
# Setting keys for merge.
setkey(dfOutGroupL2, inGroupPairing)
setkey(dfPairingResultsL2, inGroupPairing)
dfPairingResultsL2 <- merge(dfPairingResultsL2, dfOutGroupL2, 
                            by.x = "inGroupPairing", by.y = "inGroupPairing")
dfPairingResultsL2 <- data.frame(dfPairingResultsL2)
dfPairingResultsL2$latDelta <- latDelta

# Find pairings where each lineage is equal in distance to the outgroup and 
# remove these.
equalOutGroupDistCheck <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  (dfPairingResultsL1$outGroupDistance[i] - 
     dfPairingResultsL2$outGroupDistance[i])
equalOutGroupDistCheck <- which(equalOutGroupDistCheck == 0)
# If there is at least one pairing with equal outgroup distance from both 
# lineages.
if (length(equalOutGroupDistCheck) > 0) {
  # Subsetting Lineage 1 by this outgroup check.
  dfPairingResultsL1 <- dfPairingResultsL1[-equalOutGroupDistCheck, ]
  # Now subsetting Lineage 2 by this outgroup check.
  dfPairingResultsL2 <- dfPairingResultsL2[-equalOutGroupDistCheck, ]
}
# Renumber pairing again.
dfPairingResultsL1$inGroupPairing <- seq.int(nrow(dfPairingResultsL1))
dfPairingResultsL2$inGroupPairing <- seq.int(nrow(dfPairingResultsL2))

# Now some dataframe reorganization and ordering.
dfPairingResultsL1 <- 
  (dfPairingResultsL1[,c("inGroupPairing","associatedInGroupBin","inGroupDist",
                         "inGroupDistx1.3","medianLatAbs.x","latDelta",
                         "latMin.x","latMax.x","record_id.x","binSize.x",
                         "phylum_taxID.x","phylum_name.x","class_taxID.x",
                         "class_name.x","order_taxID.x","order_name.x",
                         "family_taxID.x","family_name.x","subfamily_taxID.x",
                         "subfamily_name.x","genus_taxID.x","genus_name.x",
                         "species_taxID.x","species_name.x",
                         "trimmedNucleotides.x","medianLatMap.x","medianLon.x",
                         "outGroupBin","outGroupDistance","record_id.y",
                         "binSize.y","phylum_taxID.y","phylum_name.y",
                         "class_taxID.y","class_name.y","order_taxID.y",
                         "order_name.y","family_taxID.y","family_name.y",
                         "subfamily_taxID.y","subfamily_name.y","genus_taxID.y",
                         "genus_name.y","species_taxID.y","species_name.y",
                         "trimmedNucleotides.y","pairingKey")])
colnames(dfPairingResultsL1)[2] <- "inGroupBin"
dfPairingResultsL1 <- 
  dfPairingResultsL1[order(dfPairingResultsL1$inGroupPairing), ]

# L2.
dfPairingResultsL2 <- 
  (dfPairingResultsL2[,c("inGroupPairing","associatedInGroupBin","inGroupDist",
                         "inGroupDistx1.3","medianLatAbs.x","latDelta",
                         "latMin.x","latMax.x","record_id.x","binSize.x",
                         "phylum_taxID.x","phylum_name.x","class_taxID.x",
                         "class_name.x","order_taxID.x","order_name.x",
                         "family_taxID.x","family_name.x","subfamily_taxID.x",
                         "subfamily_name.x","genus_taxID.x","genus_name.x",
                         "species_taxID.x","species_name.x",
                         "trimmedNucleotides.x","medianLatMap.x","medianLon.x",
                         "outGroupBin","outGroupDistance","record_id.y",
                         "binSize.y","phylum_taxID.y","phylum_name.y",
                         "class_taxID.y","class_name.y","order_taxID.y",
                         "order_name.y","family_taxID.y","family_name.y",
                         "subfamily_taxID.y","subfamily_name.y","genus_taxID.y",
                         "genus_name.y","species_taxID.y","species_name.y",
                         "trimmedNucleotides.y","pairingKey")])  
colnames(dfPairingResultsL2)[2] <- "inGroupBin"
dfPairingResultsL2 <- 
  dfPairingResultsL2[order(dfPairingResultsL2$inGroupPairing), ]

# Rebuilding of L1L2 dataframe with outgroups.
dfPairingResultsL1L2 <- rbind(dfPairingResultsL1, dfPairingResultsL2)
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupPairing), ]

###################
# Section 15: Identifying PseudoReplicates

# In this section, we check for the phylogenetics problem of pseudoreplication.
# We want each section of branch length to be used just once.
# This involves generating another smaller distance matrix with our selected 
# pairings only for each class.
# If a BIN from one pairing is actually closer to a BIN from another pairing 
# (as opposed to its paired sister lineage) then we would have to average the 
# results of those two pairings in the statistics section.

# To do this we can generate pairwise distance matrices with the pairing 
# BINs of each lineage only for each class.

# First, we need the BINs associated with each ingroup pairing.

# We can break pairing results down by class first in a list.
pairingResultBreakdownClass <- lapply(unique(dfPairingResultsL1L2$class_name.x), 
  function(x) dfPairingResultsL1L2[dfPairingResultsL1L2$class_name.x == x,])
# Extract BINs from each list based on class.
binClass <-  foreach(i=1:length(pairingResultBreakdownClass)) %do% 
  as.character(pairingResultBreakdownClass[[i]]$inGroupBin)
# Divide by lineage.
binClassL1 <- foreach(i=1:length(pairingResultBreakdownClass)) %do% 
  binClass[[i]][seq(1, length(binClass[[i]]), 2)]
binClassL2 <- foreach(i=1:length(pairingResultBreakdownClass)) %do% 
  binClass[[i]][seq(2, length(binClass[[i]]), 2)]
# Extract sequences from each list based on class.
sequenceClass <- foreach(i=1:length(pairingResultBreakdownClass)) %do% 
  as.character(pairingResultBreakdownClass[[i]]$trimmedNucleotides.x)
# Convert to DNAStringSet format
sequenceClassStringSet <- foreach(i=1:length(sequenceClass)) %do% 
  DNAStringSet(sequenceClass[[i]])
# Naming of DNAStringSet.
for (i in seq(from=1, to=length(sequenceClassStringSet), by = 1)) {
  names(sequenceClassStringSet[[i]]) <- binClass[[i]]
}

# Conversion to DNAbin format (no need to align since the trimmed sequences 
# should be aligned already).
pseudoRepDNABin <- foreach(i=1:length(sequenceClassStringSet)) %do% 
  as.DNAbin(sequenceClassStringSet[[i]])
# Distance Matrix Creation with BINs of each pairing lineage only, 
# using the TN93 model.
pseudoRepMatrixList <- foreach(i=1:length(sequenceClassStringSet)) %do% 
  dist.dna(pseudoRepDNABin[[i]], model = "TN93", as.matrix = TRUE, 
           pairwise.deletion = TRUE)
# Conversion to dataframe format.
pseudoRepMatrixList <- foreach(i=1:length(sequenceClassStringSet)) %do% 
  as.data.frame(pseudoRepMatrixList[[i]])
# Order Matrix by ordering of BINs in pairing results.
pseudoRepMatrixList <- foreach(i=1:length(sequenceClassStringSet)) %do%
  pseudoRepMatrixList[[i]][binClassL1[[i]], binClassL2[[i]]]
# Then we have to remove matrices (from each class) in this list which either 
# only have one pairing (2 BINs) or no pairings present since these cannot have 
# pseudoreplicates present.
pseudoRepMatrixLengthCheck <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  which(length(pseudoRepMatrixList[[i]]) < 2)
pseudoRepMatrixLengthCheck <- 
  sapply( pseudoRepMatrixLengthCheck , function(x) length( x ) )
pseudoRepMatrixLengthCheck <- which(pseudoRepMatrixLengthCheck > 0)
if (length(pseudoRepMatrixLengthCheck) > 0) {
  pseudoRepMatrixList <- pseudoRepMatrixList[-pseudoRepMatrixLengthCheck]
}

# For each column and row of each distance matrix (each matrix belonging to a 
# specific class), if a distance value is lower than the pairwise distances of 
# each pairing then we know there is another BIN which is actually closer.
# We check for this by using our pairwise distances from the pseudoRepMatrixList
# and finding the minimums of each column and row.

# First we check the minimum values for each row of each matrix.
pseudoRepMatrixMinRow <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  apply(pseudoRepMatrixList[[i]], 1, which.min)

# Then we check the minimum values for each column of each matrix.
pseudoRepMatrixMinColumn <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  apply(pseudoRepMatrixList[[i]], 2, which.min)

# If BIN distance values intersect between these two minima (row and column) 
# then we know there are no pseudoreplicates present for these BINs.
# If BIN distance values do not intersect, then we know there is a 
# pseudoreplicate BIN that is actually closer in distance than one of the
# pairing BINs.

# This command identifies indices of BINs which do not intersect, thus 
# identifying pseudoreplicates.
pseudoRepMatrixCandidates <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  which(pseudoRepMatrixMinRow[[i]] != pseudoRepMatrixMinColumn[[i]])

# Identification of BINs not intersecting on columns of each matrix.
pseudoRepMatrixCandidates <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  pseudoRepMatrixMinColumn[[i]][pseudoRepMatrixCandidates[[i]]]

# Identification of BINs not intersecting on rows of each matrix.
pseudoRepMatrixCandidates2 <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  pseudoRepMatrixMinRow[[i]][pseudoRepMatrixCandidates[[i]]]

# Naming and unlisting each list of candidates.
# First set of candidates (from columns of matrix).
pseudoRepMatrixCandidatesNames1 <- foreach(i=1:length(pseudoRepMatrixList)) %do%
  names(pseudoRepMatrixCandidates[[i]])
pseudoRepMatrixCandidates <- unlist(pseudoRepMatrixCandidatesNames1)
# Second set of candidates (from rows of matrix).
pseudoRepMatrixCandidatesNames2 <- foreach(i=1:length(pseudoRepMatrixList)) %do%
  names(pseudoRepMatrixCandidates2[[i]])
pseudoRepMatrixCandidates2 <- unlist(pseudoRepMatrixCandidatesNames2)

# Now we can get the inGroupPairing numbers associated with the candidates.

# To do this we have to create a named vector containing BINs named with pairing
# numbers of all pairing results.
pairingNumVector <- as.character(dfPairingResultsL1L2$inGroupBin)
names(pairingNumVector) <- dfPairingResultsL1L2$inGroupPairing

# Find pairing numbers of pseudoRepMatrixCandidates
if (length(pseudoRepMatrixCandidates) > 1) {
  pseudoRepPairingNum <- foreach(i=1:length(pseudoRepMatrixCandidates)) %do% 
    which(pairingNumVector == pseudoRepMatrixCandidates[[i]])
  pseudoRepPairingNum <- unlist(pseudoRepPairingNum)
  pseudoRepPairingNum <- names(pseudoRepPairingNum)
  pseudoRepPairingNum2 <- foreach(i=1:length(pseudoRepMatrixCandidates)) %do% 
    which(pairingNumVector == pseudoRepMatrixCandidates2[[i]])
  pseudoRepPairingNum2 <- unlist(pseudoRepPairingNum2)
  pseudoRepPairingNum2 <- names(pseudoRepPairingNum2)
}

# If there is at least 1 pseudoreplicate found in any of the distance matrices
if (length(pseudoRepMatrixCandidates) > 1) {
  # Creation of a dataframe for pseudoreplicates
  dfPseudoRep <- data.frame(pseudoRepPairingNum, stringsAsFactors=FALSE)
  dfPseudoRep$pseudoRepPairingNum2 <- pseudoRepPairingNum2
  # Find and remove instances where pairing numbers equal each other
  duplicatePseudoRep <- 
    which(dfPseudoRep$pseudoRepPairingNum == dfPseudoRep$pseudoRepPairingNum2)
  dfPseudoRep <- dfPseudoRep[-duplicatePseudoRep, ]
  # Order by pseudoRepPairingNum2 and remove duplicates
  dfPseudoRep <- dfPseudoRep[order(dfPseudoRep$pseudoRepPairingNum2), ]
  dfPseudoRep <- (by(dfPseudoRep, dfPseudoRep["pseudoRepPairingNum2"], 
                     head, n=1))
  dfPseudoRep <- Reduce(rbind, dfPseudoRep)
  # converting to numeric type for easy ordering
  dfPseudoRep$pseudoRepPairingNum <- as.numeric(dfPseudoRep$pseudoRepPairingNum)
  dfPseudoRep$pseudoRepPairingNum2 <- 
    as.numeric(dfPseudoRep$pseudoRepPairingNum2)
  # Remove instances of duplicates across columns, ex: 2,1 and 1,2
  dfPseudoRep <- unique(t(apply(dfPseudoRep, 1, sort)))
  dfPseudoRep <- data.frame(dfPseudoRep)
  colnames(dfPseudoRep)[1] <- "pseudoRepPairingNum"
  colnames(dfPseudoRep)[2] <- "pseudoRepPairingNum2"
  dfPseudoRep <- dfPseudoRep[order(dfPseudoRep$pseudoRepPairingNum), ]
}

# These identified pseudoreplicates will now have their signed 
# relative outgroup distances averaged in the statistics section.

###################
# Section 16: Relative OutGroup Distance Determination for Pairings
# In this section, signed relative outgroup distance is determined for each 
# individual pairing.

# Positive relative distance indicates a pairing where the lower latitude 
# lineage also has the further outgroup distance, thus supporting the hypothesis
# of higher evolutionary rates in the tropics.

# Negative would be where the lower latitude lineage has the smaller outgroup 
# distance.

# Create variables that hold outgroup distances.
testOutGroupDist <- as.numeric(dfPairingResultsL1L2$outGroupDistance)
# Stored as a vector.
# Set trueTestOutGroupDist to null to ensure empty variable.
trueTestOutGroupDist<-NULL
# Counter for trueTestOutGroupDist
x = 1
# For loop that increments by 2 to divide adjacent outgroup distances found in
# outGroupDist. This will give our relative genetic distances in the vector 
# testOutGroupDist.
for (i in seq(from = 1, to = length(testOutGroupDist), by = 2)){
  if (testOutGroupDist[i] > testOutGroupDist[i + 1]) {
    trueTestOutGroupDist[x] <- testOutGroupDist[i] / testOutGroupDist[i + 1]
  }
  else if (testOutGroupDist[i + 1] > testOutGroupDist[i]) {
    trueTestOutGroupDist[x] = testOutGroupDist[i + 1] / testOutGroupDist[i]
  }
  else if (testOutGroupDist[i + 1] == testOutGroupDist[i]) {
    trueTestOutGroupDist[x] = testOutGroupDist[i + 1] / testOutGroupDist[i]
  }
  x <- x + 1
}

# Checking to see which latitudes are lower for lineage 1 compared to lineage2.
lowerLatL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfPairingResultsL1$medianLatAbs.x[i] < 
          dfPairingResultsL2$medianLatAbs.x[i]) 
# Checking to see which genetic distances are greater for lineage 1 compared to 
# lineage2.
greaterDistanceL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfPairingResultsL1$outGroupDistance[i] > 
          dfPairingResultsL2$outGroupDistance[i])
# Then an integer of 1 will be assigned to those in Lineage1 that meet those 
# criteria. Meeting these criteria is defined as a success.
successValL1 <- mapply(intersect, lowerLatL1, greaterDistanceL1)
# This will give the index of each pairing relative to Lineage 1 that met 
# those criteria.
successValL1 <- which(successValL1 > 0) 

# Do the same process for lineage 2 compared to lineage1.
lowerLatL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfPairingResultsL2$medianLatAbs.x[i] < 
          dfPairingResultsL1$medianLatAbs.x[i]) 
greaterDistanceL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfPairingResultsL2$outGroupDistance[i] > 
          dfPairingResultsL1$outGroupDistance[i])
successValL2 <- mapply(intersect, lowerLatL2 ,greaterDistanceL2)
successValL2 <- which(successValL2 > 0) 

# Append both success lists together, representing our total number of successes.
successOverall <- append(successValL1, successValL2)
# Sort the list.
successOverall <- sort(successOverall)

# Also making a vector for pairings that were not successes.
failureOverall <- dfPairingResultsL1$inGroupPairing[-successOverall]

# Can then subset our trueTestOutGroupDist vector containing relative distances 
# with our successes.
relativeDistPos <- trueTestOutGroupDist[successOverall]
# Then subset based on which relative distances were not successes.
relativeDistNeg <- trueTestOutGroupDist[failureOverall]
# These values will also be assigned a negative value.
relativeDistNeg <- relativeDistNeg * -1
# Also identifying which pairing was positive and which was negative.
names(relativeDistPos) <- paste0(successOverall)
names(relativeDistNeg) <- paste0(failureOverall)

# Creating a new named vector with positive and negative values appended 
# together but with pairing numbers as names.
relativeDistOverall <- append(relativeDistNeg, relativeDistPos)
# Make relativeDistOverall into a dataframe.
dfRelativeDist <- data.frame(relativeDistOverall)
# Adding rownames variable column of dfRelativeDistOverall.
dfRelativeDist$variable <- rownames(dfRelativeDist)
# Change to numeric.
dfRelativeDist$variable <- as.numeric(dfRelativeDist$variable)
# Order by this column.
dfRelativeDist <- dfRelativeDist[order(dfRelativeDist$variable), ]
# Adding relative distances to a new column called value.
colnames(dfRelativeDist)[1] <- "value"
# Creating another value called sign to determine which is positive and which is
# negative.
dfRelativeDist[["sign"]] = ifelse(dfRelativeDist[["value"]] >= 0, "positive", 
                                  "negative")

# Appending relative outgroup distance results to ResultsSummary dataframe and 
# ResultsL1L2 dataframe.
dfPairingResultsL1$relativeOutGroupDist <- dfRelativeDist$value
dfPairingResultsL2$relativeOutGroupDist <- dfRelativeDist$value
dfPairingResultsL1L2 <- rbind(dfPairingResultsL1, dfPairingResultsL2)
dfPairingResultsL1L2 <- 
  dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupPairing), ]

# Also creating another dataframe called pairing results summary which will give
# a quick idea of pairing results and make the results easier to read. This will
# include relative outgroup distance.
dfPairingResultsSummary <- 
  (dfPairingResultsL1L2[,c("inGroupPairing","relativeOutGroupDist","inGroupBin",
                           "inGroupDist","inGroupDistx1.3","outGroupDistance",
                           "outGroupBin","medianLatAbs.x","latDelta","latMin.x",
                           "latMax.x","class_name.x","order_name.x",
                           "family_name.x","genus_name.x","species_name.x",
                           "trimmedNucleotides.x","class_name.x","order_name.y",
                           "family_name.y","genus_name.y","species_name.y",
                           "trimmedNucleotides.y","pairingKey")])
colnames(dfPairingResultsSummary)[8] <- "inGroupMedianLatAbs"
colnames(dfPairingResultsSummary)[10] <- "inGroupMinLat"
colnames(dfPairingResultsSummary)[11] <- "inGroupMaxLat"
colnames(dfPairingResultsSummary)[12] <- "inGroupClass"
colnames(dfPairingResultsSummary)[13] <- "inGroupOrder"
colnames(dfPairingResultsSummary)[14] <- "inGroupFamily"
colnames(dfPairingResultsSummary)[15] <- "inGroupGenus"
colnames(dfPairingResultsSummary)[16] <- "inGroupSpecies"
colnames(dfPairingResultsSummary)[17] <- "inGroupNucleotides"
colnames(dfPairingResultsSummary)[18] <- "outGroupClass"
colnames(dfPairingResultsSummary)[19] <- "outGroupOrder"
colnames(dfPairingResultsSummary)[20] <- "outGroupFamily"
colnames(dfPairingResultsSummary)[21] <- "outGroupGenus"
colnames(dfPairingResultsSummary)[22] <- "outGroupSpecies"
colnames(dfPairingResultsSummary)[23] <- "outGroupNucleotides"
colnames(dfPairingResultsSummary)[24] <- "pairingKey"
# Reordering.
dfPairingResultsSummary <- 
  dfPairingResultsSummary[order(dfPairingResultsSummary$inGroupPairing),]

# Merge RelativeDist dataframe to a pairing results dataframe to grab the class 
# names so results can be broken down by class.
dfRelativeDist <- merge(dfRelativeDist, dfPairingResultsL1, 
                        by.x = "variable", by.y = "inGroupPairing")
dfRelativeDist <- (dfRelativeDist[,c("variable","value","sign","class_name.x")])

###################
# Section 17: PseudoReplicate Relative Outgroup Distance Averaging
# In this section we average pseudoreplicate relative outgroup distances so they
# can be added to the relative distance results and get included in the 
# statistics section.

# Again, if at least 1 pseudoreplicate is present in the dfPseudoRep dataframe, 
# then:
if (length(pseudoRepMatrixCandidates) > 1) {
  # First we have to add the signed relative outgroup distances to dfPseudoRep, 
  # merging dfRelativeDist to do this:
  dfPseudoRep <- merge(dfPseudoRep, dfRelativeDist, 
                       by.x = "pseudoRepPairingNum", by.y = "variable")
  dfPseudoRep <- merge(dfPseudoRep, dfRelativeDist, 
                       by.x = "pseudoRepPairingNum2", by.y = "variable")
  dfPseudoRep <- 
    (dfPseudoRep[,c("pseudoRepPairingNum","value.x",
                    "pseudoRepPairingNum2","value.y","class_name.x.x")])
  colnames(dfPseudoRep)[2] <- "relativeDist1"
  colnames(dfPseudoRep)[4] <- "relativeDist2"
  
  # Also ordering dfPseudoRep.
  dfPseudoRep <- dfPseudoRep[order(dfPseudoRep$pseudoRepPairingNum), ]
  # Rounding.
  dfPseudoRep$relativeDist1 <- round(dfPseudoRep$relativeDist1, 6)
  dfPseudoRep$relativeDist2 <- round(dfPseudoRep$relativeDist2, 6)
  
  # Now subtracting 1 from positive values and adding one to negative values for 
  # averaging. This is to address an issue where averages of pairings differing 
  # in sign were producing decimal results that were difficult to interpret 
  # biologically.
  # RelativeDist1:
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if (dfPseudoRep$relativeDist1[i] > 0) {
      dfPseudoRep$relativeDist1[i] <- dfPseudoRep$relativeDist1[i] - 1
    }
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if (dfPseudoRep$relativeDist1[i] < 0) {
      dfPseudoRep$relativeDist1[i] <- dfPseudoRep$relativeDist1[i] + 1
    }
  # RelativeDist2:
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if (dfPseudoRep$relativeDist2[i] > 0) {
      dfPseudoRep$relativeDist2[i] <- dfPseudoRep$relativeDist2[i] - 1
    }
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if (dfPseudoRep$relativeDist2[i] < 0) {
      dfPseudoRep$relativeDist2[i] <- dfPseudoRep$relativeDist2[i] + 1
    }
  
  # Breaking pseudoreplicates down to a list.
  pseudoRepList <- 
    lapply(unique(dfPseudoRep$pseudoRepPairingNum), function(x) 
      dfPseudoRep[dfPseudoRep$pseudoRepPairingNum == x,])
  
  # Number of pseudoreplicates per pairing.
  pseudoRepLength <- 
    sapply( pseudoRepList , function (x) length( x$relativeDist2 ) )
  # Distances and names associated with pseudoreplicates only:
  pseudoRepRelativeDist2 <- 
    sapply( pseudoRepList , function (x) ( x$relativeDist2 ) )
  pseudoRepRelativeDist2Names <- 
    sapply( pseudoRepList , function (x) ( x$pseudoRepPairingNum2 ) )
  pseudoRepRelativeDist2Names <- pseudoRepRelativeDist2Names
  # Distances and names associated with ingroup pairings only.
  pseudoRepRelativeDist1 <- unique(dfPseudoRep$relativeDist1)
  pseudoRepRelativeDist1Names <- 
    sapply( pseudoRepList , function (x) as.character( x$pseudoRepPairingNum ) )
  pseudoRepRelativeDist1Names <- 
    sapply( pseudoRepRelativeDist1Names , function (x) unique( x ) )
  # Append these distances together for averaging using map.
  # Map will append to a list format.
  pseudoRepAllRelativeDist <- 
    (Map(c, pseudoRepRelativeDist2, pseudoRepRelativeDist1))
  pseudoRepAllRelativeDistNames <- 
    Map(c, pseudoRepRelativeDist2Names, pseudoRepRelativeDist1Names)
  pseudoRepAllRelativeDistNames <- unname(pseudoRepAllRelativeDistNames)
  
  # Now we can finally average the relative distances based on the values in the 
  # pseudoRepAllRelativeDist list.
  pseudoRepAverage <- sapply( pseudoRepAllRelativeDist , function(x) mean( x ) )
  
  # Making another dataframe with the averages for the pseudoreplicates called 
  # dfPseudoRepAverage.
  dfPseudoRepAverage <- data.frame(pseudoRepAverage)
  
  # Now subtracting 1 if negative or adding 1 if positive to averages.
  # Once again this will ensure that all pseudoreplicate average values are 
  # above 1 or below -1.
  foreach(i=1:nrow(dfPseudoRepAverage)) %do% 
    if (dfPseudoRepAverage$pseudoRepAverage[i] < 0) {
      dfPseudoRepAverage$pseudoRepAverage[i] <- 
        dfPseudoRepAverage$pseudoRepAverage[i] - 1
    }
  foreach(i=1:nrow(dfPseudoRepAverage)) %do% 
    if (dfPseudoRepAverage$pseudoRepAverage[i]>0) {
      dfPseudoRepAverage$pseudoRepAverage[i] <- 
        dfPseudoRepAverage$pseudoRepAverage[i] + 1
    }
  
  # Adding another column for the pairings associated with each average.
  dfPseudoRepAverage$variable <- pseudoRepAllRelativeDistNames
  TrimFunction <- function (x) sub("[c(]","", x)
  dfPseudoRepAverage$variable <- TrimFunction(dfPseudoRepAverage$variable)
  colnames(dfPseudoRepAverage)[1] <- "value"
  # Also adding the signs of each average for plotting later on.
  dfPseudoRepAverage[["sign"]] = ifelse(dfPseudoRepAverage[["value"]] >= 0, 
                                        "positive", "negative")
  # Now let's subset dfRelativeDist to remove pairings that were averaged.
  # First, making a vector with all pairings averaged.
  pseudoRepPairings <- append(dfPseudoRep$pseudoRepPairingNum, 
                              dfPseudoRep$pseudoRepPairingNum2)
  pseudoRepPairings <- pseudoRepPairings[!duplicated(pseudoRepPairings)]
  # Then, subsetting dfRelativeDist based on this.
  dfRelativeDist <- 
    dfRelativeDist[setdiff(dfRelativeDist$variable, pseudoRepPairings),]  
  # Adding class names to dfPseudoRepAverage.
  dfPseudoRepAverage$className <- 
    sapply( pseudoRepList, function(x) unique( x$class_name.x.x ) )
  # Column name change for dfRelativeDist.
  colnames(dfRelativeDist)[4] <- "className"
  # Now dfPseudoRepAverage will be appended to the relativeDist dataframe.
  dfRelativeDist <- rbind(dfRelativeDist,dfPseudoRepAverage)
} else {
  # Column name change for dfRelativeDist.
  colnames(dfRelativeDist)[4] <- "className"
}

###################
# Section 18: Statistical Analysis of Pairings

# In this section we are performing Binomial Tests and Wilcoxon Tests on the 
# relative distance results of the pairings to investigate the hypothesis of 
# higher rates of molecular evolution at lower latitudes.

# Performing a binomial test on our pairings to test to see if pairings with 
# lower latitude AND greater outgroup distance are more prevalent than the null 
# expectation of 50%, the null expectation being that pairings with high 
# lat/high outgroup dist are equally prevalent compared with low lat/high 
# outgroup dist pairings.

# Binomial test on relative outgroup distances for sister pairs.
# Number of successes defined as number of positive values in the dfRelativeDist
# dataframe:
successOverall <- length(which(dfRelativeDist$sign == "positive"))
# Total number of trials is number of rows of dfRelativeDist dataframe.
# Null expectation set to 50%.

# Binomial test for all classes.
binomialTestOutGroup <- binom.test(successOverall, nrow(dfRelativeDist), p= 0.5)

# Binomial test for each class separately.

# Break relative outgroup distances down by class.
classBinomList <- lapply(unique(dfRelativeDist$className), function(x) 
  dfRelativeDist[dfRelativeDist$className == x,])

# Total number of observations for each class.
classNumObservations <- foreach(i=1:length(classBinomList)) %do% 
  length(classBinomList[[i]]$variable)

# Number of successes for each class.
classNumSuccesses <- foreach(i=1:length(classBinomList)) %do% 
  which(classBinomList[[i]]$sign == "positive" ) 
classNumSuccesses <- foreach(i=1:length(classBinomList)) %do% 
  length(classNumSuccesses[[i]])

# Binom test for each class.
binomClass <- foreach(i=1:length(classNumObservations)) %do% 
  binom.test(classNumSuccesses[[i]], classNumObservations[[i]], p=0.5)

# p-value extraction for binomial test from each class.
pvalBinomialTotal <- binomialTestOutGroup$p.value
pvalBinomialClass <- foreach(i=1:length(classBinomList)) %do% 
  binomClass[[i]]$p.value

# Creation of a dataframe for p-values.
allClassVariable <- "AllClasses"
classNames1 <- foreach(i=1:length(classBinomList)) %do% 
  unique(classBinomList[[i]]$className)
classNames1 <- unlist(classNames1)
classNames2 <- append(allClassVariable, classNames1)
pValBinomial <- append(pvalBinomialTotal, pvalBinomialClass)
dfPVal <- data.frame(classNames2)
dfPVal$pValueBinomial <- round(as.numeric(pValBinomial), digits = 6)

# Wilcoxon Test

# Next we do a Wilcoxon test on all of the signed relative distances 
# (pos or neg) from the value column of the relativeDist dataframe to compare 
# the median for a significant difference from the null expectation of zero.
# This test will consider both the magnitude and direction but is also 
# non-parametric.

# Wilcoxon test for all classes together.
wilcoxTestOutGroup <- wilcox.test(dfRelativeDist$value, mu=0)

# Wilcoxon test for each class separately:

# Relative outgroup values for each class to be fed into test.
classValue <-  foreach(i=1:length(classBinomList)) %do% 
  (classBinomList[[i]]$value)

# Wilcoxon test for each class.
wilcoxonClass <- foreach(i=1:length(classValue)) %do% 
  wilcox.test(classValue[[i]], mu=0) 

# p-value extraction for Wilcoxon test.
pvalWilcoxon <- wilcoxTestOutGroup$p.value
pvalWilcoxonClass <- foreach(i=1:length(wilcoxonClass)) %do% 
  (wilcoxonClass[[i]]$p.value)
pValWilcoxonTotal <- append(pvalWilcoxon, pvalWilcoxonClass)

# Addition of Wilcoxon p-values to dfPVal dataframe.
dfPVal$pValueWilcoxon <- round(as.numeric(pValWilcoxonTotal), digits = 6)

###################
# Section 19: Plotting of Relative Outgroup Distance Results

# In this section, we do plotting of our relative distances based on pairing 
# number. Points will plot red if the value is below 0 (meaning not a success) 
# and blue if above 0 (success!).

# Make the variable column a factor for dfRealtiveDistOverall so ggplot classes 
# pairings correctly.
dfRelativeDist$variable <- factor(dfRelativeDist$variable, 
                                  levels = dfRelativeDist$variable)

# Breaking dfRelativeDist down to a list by class so plots can be generated 
# according to class.
relativeDistClass <- lapply(unique(dfRelativeDist$className), function(x) 
  dfRelativeDist[dfRelativeDist$className == x,])

# Variables for the title of each plot; this will include name of the class and 
# associated p-values.
classTitle <- foreach(i=1:length(classNames1)) %do% 
  paste("Relative Outgroup Distances of Each Pairing of", classNames1[i])
pValBinomialTitle <- foreach(i=1:length(pvalBinomialClass)) %do%
  paste("Binomial Test p-value:", round(pvalBinomialClass[[i]], digits = 6))
pvalWilcoxonTitle <- foreach(i=1:length(pvalBinomialClass)) %do%
  paste("Wilcoxon Test p-value:", round(pvalWilcoxonClass[[i]], digits = 6))

# Plot of signed relative distance per pairing using ggplot:

# Separate plots will be generated for each class.

foreach(i = 1:length(relativeDistClass)) %do%
  (print(ggplot(relativeDistClass[[i]], aes(x = variable, y = value, 
                                            color = sign))
         + geom_point(stat = "identity", size = 2.5)
         + theme(
           text = element_text(size = 13),
           axis.text.x = element_text(
             face = "bold",
             angle = 90,
             vjust = 1
           )
         )
         + ggtitle(
           paste0(
             classTitle[i],
             "\n",
             pValBinomialTitle[i],
             "\n",
             pvalWilcoxonTitle[i]
           )
         )
         + labs(x = "Pairing Number", y = "Signed Relative OutGroup Distance", 
                color = "sign")
  ))


# Click on zoom icon in the bottom right hand corner to see plot more clearly.

# Plots can be scrolled through using the arrow icons in the bottom right 
# hand corner.

# Plots that are important to keep can then be exported to any directory using 
# the Export function in the viewer (bottom right hand of RStudio).

# If "Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) 
# : invalid graphics state" message appears, use this command:
# dev.off()

###################
# Section 20: Map Plotting of Latitude Separated Pairings using Plotly

# In this section, we can plot each pairing result on a customizable world map 
# where each pairing is color coded and located on the map according to its 
# respective median latitude and longitude.

# Info can be found here on this: https://plot.ly/r/scatter-plots-on-maps/

# Code adapted from Plotly website.
# Note there are more options for map customization that can be found on plotly.
# This just represents a basic world map representation.

# New variable for basic map layout characteristics:
# (using eckert4 projection which will reduce area near poles and expand area 
# near equator.)

mapLayout <- list(
  showland = TRUE,
  showlakes = TRUE,
  showcountries = TRUE,
  showocean = TRUE,
  countrywidth = 0.5,
  landcolor = toRGB("grey90"),
  lakecolor = toRGB("white"),
  oceancolor = toRGB("white"),
  projection = list(type = 'eckert4'),
  lonaxis = list(
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  ),
  lataxis = list(
    showgrid = TRUE,
    gridcolor = toRGB("gray40"),
    gridwidth = 0.5
  )
)

# New dataframe column with data for hovering over points on the map.
# Can add more columns to hover if you want more detail on the map for each
# pairing.
dfPairingResultsL1L2$hover <- 
  paste("PairNum:",dfPairingResultsL1L2$inGroupPairing,
        "Ingroup:",dfPairingResultsL1L2$inGroupBin, 
        dfPairingResultsL1L2$species_name.x, 
        round(dfPairingResultsL1L2$inGroupDist, 6), 
        "Outgroup:", dfPairingResultsL1L2$outGroupBin,
        dfPairingResultsL1L2$species_name.y,
        sep = "<br>")

# Title for the map:

# First taking all unique class names to name the map.
uniqueClass <- unique(dfPairingResultsL1L2$class_name.x)
# Pasting them together.
classConcatenate <- paste(uniqueClass , collapse = ', ')

# Map title command.
mapTitle <- paste("Map of Latitude Separated Sister Pairings for", 
                  classConcatenate)

# This command will ensure the pairing results dataframe can be read by plotly.
attach(dfPairingResultsL1L2)

# This command will show a map organized by the pairing number with each pairing
# number having its own distinctive color on a spectrum of purple 
# (shortest distance) to red (longest distance). The map can become crowded with
# large numbers of pairings, so individual maps may need to be manually created
# for subsets of the pairing results.

plot_ly(dfPairingResultsL1L2, lat = medianLatMap.x, lon = medianLon.x, 
        text = hover, color = as.ordered(inGroupPairing), colors = "Spectral", 
        mode = "markers+lines", type = 'scattergeo') %>%
  layout(title = paste0(mapTitle) , geo = mapLayout)

# A map can be generated by class, with each class having its own distinctive 
# color. Due to limitations with plotly, lines will not connect pairings:

plot_ly(dfPairingResultsL1L2, lat = medianLatMap.x, lon = medianLon.x, 
        text = hover, color = class_name.x, mode = "markers", 
        type = 'scattergeo') %>%
  layout(title = paste0(mapTitle) , geo = mapLayout)

# ***You will need to click on the icon in the viewer (bottom right corner) that
# says "show in new window" (little box with arrow beside the refresh icon). 
# Unfortunately, this does not show the actual map directly in Rstudio.
# The map will appear in a web browser window, though you don't have to be 
# online to do this.***

###################
# Section 21: Posting of Map to Plotly Server for Online Viewing
# ***Only important if you want to post maps online***

# For uploading to plotly server for online viewing of map.

# You will first have to create a plotly account to do this:
# https://plot.ly/

# Note there is a limit of one plot for the free version of the Plotly account,
# and the plot is public, meaning other people on plotly can view the plot, 
# though it is not easily found on the website without the direct link.

# To obtain additional plots on the server, you have to pay for a package.

# Run these commands for uploading user details, enter username and API key in 
# the empty quotations to run commands:
# (obtained from making an account and in settings of account details) 

# Sys.setenv("plotly_username"="") 
# Sys.setenv("plotly_api_key"="")

# Run these commands to make plot as a variable:

# plot <- plot_ly(dfPairingResultsL1L2, lat = medianLatMap.x, lon = medianLon.x,
# text = hover, color = as.ordered(inGroupPairing), colors = "Spectral", 
# mode = "markers+lines", type = 'scattergeo') %>%
# layout(title = 'Latitude Separated Sister Pairings', geo = mapLayout)

# Run this command for posting of map to plotly server 
# (can rename in quotations to a name you prefer):

# plotly_POST(plot, filename = "LatitudeSeparatedPairingsMap")

###################
# Section 22: Output of Pairing Results to TSV or CSV
# ***Only important if you want to export results outside of R.***

# Uncomment the appropriate commands for either TSV or CSV output.

# First create this dataframe for dfPairingResultsL1L2, this will allow output 
# of this file:

# dfPairingResultsOut <- data.frame(lapply(dfPairingResultsL1L2, as.character),
#   stringsAsFactors=FALSE)

# Alternatively for results summary:

# dfPairingResultsOut <- data.frame(lapply(dfPairingResultsSummary, 
#   as.character), stringsAsFactors=FALSE)

# Defining another variable to give a unique name to the CSV.
# Name would be what you specify it to be. R will prompt you to insert a name 
# for the file in the console once this command is uncommented:
# filename <- readline(prompt="")

# Then you can uncomment one of these write.table commands to output.
# ***will output file to current working directory of R***

# To output results to CSV format:

# write.table(dfPairingResultsOut, file=paste(filename, ".csv", sep=""), 
# quote=FALSE, sep=',', col.names = NA)

# To output results to TSV format:

# write.table(dfPairingResultsOut, file=paste(filename, ".tsv", sep=""), 
# quote=FALSE, sep='\t', col.names = NA)  
