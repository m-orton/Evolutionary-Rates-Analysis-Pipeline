#################
#Authored by Matthew Orton

#This program will allow for the generation of latitudinally separated sister pairings and associated outgroupings from ANY taxa 
#(provided they have a suitable reference sequence) and geographical region found on BOLD in one streamlined R pipeline!
#In this new and improved iteration, data is translated directly from a BOLD tsv file of the users choosing to a dataframe in R
#The generated sister pairs and outgroups can then be written to a csv or tsv and the file will appear in the current working 
#directory of R

##################
#A few important tips:

#There are two options for parsing the tsv, you can either download the tsv directly from the BOLD API OR you can use a tsv 
#you have previously downloaded

#Larger taxa, ex: aves can potentially take from several minutes to hours in order to process the tsv and run the sequence alignment. 
#This can potentially consume a lot of working memory. I wouldnt run any taxa that is very large, (insecta for example) until 
#we know what kind of memory resources it will consume.

#When testing a taxa for the first time you will have to ensure that you have a suitable reference sequence for it 
#and ensure that it is inserted into the dfRefSeq dataframe
#Once it is inserted you will then have to modify which row of the dfRefSeq dataframe is being used further 
#ex: alignmentRef <- as.character(dfRefSeq$nucleotides[1]) would correspond to first row of dfRefSeq
#alignmentRef <- as.character(dfRefSeq$nucleotides[2]) would correspond to the second row of dfRefSeq

#Some tips for using the BOLD API since this is what is used to grab the relevant data we need: 
#To see details on how to use the bold API, go to http://www.boldsystems.org/index.php/resources/api?type=webservices
#Can add additional restrictions to url, for example &instituiton=Biodiversity Institute of Ontario|York University or 
#&marker=COI-5P if you want to specifiy an institution or specific genetic marker in the url
#geo=all in the url means global but geo can be geo=Canada for example
#Can use | modifier in the url, for example &geo=Canada|Alaska would give data for both Canada and alaska, 
#or taxon=Aves|Reptilia would yield give data for both Aves and Reptilia

#I recommend using RStudio: https://www.rstudio.com/home/, provides a nicer interface than just using R
#Large matrices and dataframes may not show all columns, you can overcome this by typing the command: utils::View() 
#where you would insert the dataframe or matrix you want in the brackets of that command

#Important dataframes:
#dfMatchOverallBest is the finalized dataframe that contains all of the finalized pairings and outgroupings
#dfMatchOverall is the dataframe of ingroup pairings before bins with multiple pairings are eliminated based 
#on distance and before the problems of pseudoreplication or latitudinal range is addressed, it is used as a reference to 
#see which pairings are intially created
#dfMatchOverallLineage1 and Lineage2 represent dataframes for each lineage of each pairing
#dfInitial is the dataframe first produced by the import from BOLD and is trimmed by lat, bin_uri etc.
#dfLatLon just contains relevant information for each bin: bin size, maximum lat, median lat, minimum lat, median lon
#dfBestOutGroupL1 and L2 contain the associated outgroupings only (for each lineage) but each one does have a column for 
#which pairing its associated with
#dfLatitudeDistance and dfGeneticDistance are matrices converted into dataframes that show all possible distances between bins
#dfAllSeq is the dataframe that contains sequence data for both consensus sequences and nonconsensus (bins with one member to them) 
#sequences
#dfGeneticDistanceStack is dfGeneticDistance with all columns concatenated into one long column, it is used to grab index numbers 
#for each pairing
#dfDistancePair represents pairwise distances between lineages in the final pairings only, it used to determine if 
#pseduoreplication is present in some pairings
#dfDistancePairTotal represents the indexes and bins of the pseudoreplicate bins
#dfConsensus and dfNonconsensus are dataframes with the consensus and nonconcensus sequences respectively but is also 
#similar to dfLatLon in the data it contains
#dfRelativeDist shows the relative distances to the outgroup for each pairing
#dfRefSeq shows vairous taxa with a suitable reference sequence that has been found for them

#################
#Packages required

#We need the foreach package for several functions that require iteration over dataframe rows
install.packages("foreach")
library(foreach)
#For genetic distance determination using the TN93 model, we use the ape package
install.packages("ape")
library(ape)
#Speeds up parsing of the tsv file with read_tsv function
install.packages("readr")
library(readr)
#For sequence alignments we need the biostrings (DNAStringSet function) and msa packages, run each of these 
#commands individually, sometimes it skips the libraries
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")
biocLite("msa")
library("msa")
#For overlapping latitude regions we need the Desctools package
install.packages("DescTools")
library(DescTools)
#Also adding data tables for table merging in the outgrouping section
install.packages("data.table")
library(data.table)
#For plotting of relative outgroup distances between lineages we will also need ggplot2
require(ggplot2)

#################
#R Commands:

#TSV Parsing

#First we download the TSV and convert it into a dataframe, this URL is what is modified by the user and will 
#determine the taxa, geographic region etc.
dfInitial <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Sphingidae&geo=all&format=tsv")

#If you want to run pre downloaded BOLD tsv's to avoid downloading of the same tsv multiple times, this will 
#let you choose a path to that tsv and parse
#tsvParseDoc <- file.choose()
#dfInitial <- read_tsv(tsvParseDoc)

##############
#Dataframe Filtering and Reorganization

#Filtering this df according to the relevant columns we need
dfInitial <- (dfInitial[,c("recordID","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","lat","lon","nucleotides")])
colnames(dfInitial)[1] <- "record_id"

#Removing sequences with no latitude values, filtering according to lat since we only really need lat for the analysis
containLat <- grep( "[0-9]", dfInitial$lat)
dfInitial<-dfInitial[containLat,]

#Next we have to convert lat column to num instead of chr type, this will become important later on for median latitude determination
latNum <- with(dfInitial, as.numeric(as.character(lat))) 
dfInitial$latNum <- latNum

#Can do lon as well to get numeric values instead of characters
lonNum <- with(dfInitial, as.numeric(as.character(lon))) 
dfInitial$lonNum <- lonNum

#First identifying missing bins and eliminating rows with missing bin_uri's since bin is a big indicator of sequence quality
#Grep by colon since every record with a bin identifier will have this
containBin <- grep( "[:]", dfInitial$bin_uri)
dfInitial<-dfInitial[containBin,]

#Getting rid of any records that dont have sequence data (sometimes there are a few)
containNucleotides <- grep( "[ACGT]", dfInitial$nucleotides)
dfInitial<-dfInitial[containNucleotides,]

#Really high dash content will affect the alignment so we can get rid of sequences with greater than 5% dashes
#This will give the number of positions where an N is found for each sequence
containDash <- gregexpr( "[-]", dfInitial$nucleotides)
#We then go through each sequence and see if the number of dashes is greater than 5% of total sequence length 
#(0.05 can easily modified to add more or less stringency)
containDash <- foreach(i=1:nrow(dfInitial)) %do% which((containDash[[i]]/nchar(dfInitial$nucleotides[i])>0.01))
dashcheck <- sapply( containDash , function (x) length( x ) )
dashcheck <- which(dashcheck>0)
#Subset out these higher N content sequences
dfInitial<-dfInitial[-dashcheck,]

#Modifying Bin column slightly to remove "BIN:"
dfInitial$bin_uri <- substr(dfInitial$bin_uri, 6 , 13)

#Dataframe Reorganization
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","latNum","lonNum")])

############
#Bin Stats and Median Latitude/Longitude Determination per bin

#First we can make a smaller dataframe with the columns we want for each bin - bin_uri, latnum, lonnum, record_id, 
#if we didnt do this, the binList would consume a huge amount of memory
dfBinList <- (dfInitial[,c("record_id","bin_uri","latNum","lonNum","nucleotides")])
#Create groupings by bin with each grouping representing a different bin_uri
#Each element of this list represents a bin with subelements representing the various columns of the initial 
#dataframe created and the information is grouped by bin 
binList <- lapply(unique(dfBinList$bin_uri), function(x) dfBinList[dfBinList$bin_uri == x,])

#Now to determine a median latitude for each bin
medianLat <- sapply( binList , function(x) median( x$latNum ) )

#We also need a median longitude for each if we are going to plot on a map for a visual interface
medianLon <- sapply( binList , function(x) median( x$lonNum ) )

#we can also take a few other important pieces of data regarding each bin using sapply including number 
#of record_ids to a bin and latitudinal min and max of each bin
latMin <- sapply( binList , function(x) min( x$latNum ) )
latMax <- sapply( binList , function(x) max( x$latNum ) )
binSize <- sapply( binList , function (x) length( x$record_id ) )

#Dataframe of our median lat values, this will be used in our final dataframe
dfLatLon <- data.frame(medianLat)

#Adding bin_uri, median longitude, latMin, latMax and binSize to dataframe with medianLat
dfLatLon$bin_uri <- c(unique(dfInitial$bin_uri))
dfLatLon$medianLon <- c(medianLon)
dfLatLon$latMin <- c(latMin)
dfLatLon$latMax <- c(latMax)
dfLatLon$binSize <- c(binSize)

#Can also convert to 180 degree latitude scale, and longitude too
dfLatLon$medianLat <- dfLatLon$medianLat + 90
dfLatLon$latMax <- dfLatLon$latMax + 90
dfLatLon$latMin <- dfLatLon$latMin + 90
dfLatLon$medianLon <- dfLatLon$medianLon + 180

#Merging LatLon to BinList for the sequence alignment step
dfBinList <- merge(dfBinList, dfLatLon, by.x = "bin_uri", by.y = "bin_uri")

###############
#Selecting One Sequence Per Bin

#First filtering out any bins with high N content before we pick a bin
#To avoid using sequences with high N content we can eliminate sequences with greater than 1% total N content

#This will give the number of positions where an N is found for each sequence
containN <- gregexpr( "[N]", dfInitial$nucleotides)
#We then go through each sequence and see if the number of N's is greater than 1% of total sequence length 
#(0.01 can easily modified to add more or less stringency)
containN <- foreach(i=1:nrow(dfInitial)) %do% which((containN[[i]]/nchar(dfInitial$nucleotides[i])>0.01))
ncheck <- sapply( containN , function (x) length( x ) )
ncheck <- which(ncheck>0)
#Subset out these higher N content sequences
dfInitial<-dfInitial[-ncheck,]

#Then picking one sequence per bin
dfInitial <- by(dfInitial, dfInitial["bin_uri"], head, n=1)
dfInitial <- Reduce(rbind, dfInitial)

#Creation of a new dataframe dfAllSeq with all data needed for the alignment and all further steps
dfAllSeq <- merge(dfLatLon, dfInitial, by.x = "bin_uri", by.y = "bin_uri", all.Y=TRUE)
dfAllSeq <- (dfAllSeq[,c("bin_uri","binSize","record_id","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLat","latMin","latMax","medianLon")])
dfAllSeq$ind <- row.names(dfAllSeq)

#############
#Multiple Sequence Alignment of All Sequences with Reference Sequence

#Lets first start off by identifying our reference sequence
#we can make a smaller dataframe with the name of the taxa as one column and the sequence as another column
#This dataframe can be populated further with more reference sequences as we find more in the literature
#Can call this dfRefSeq, for the sake of testing just using a standard length (658 bp) COI-5P sequence from Sphingidae for now
dfRefSeq <- data.frame(taxa = c("Sphingidae"),
                       nucleotides = c("GGTATCTGAGCTGGAATAGTTGGTACTTCATTAAGATTATTAATTCGGGCAGAATTAGGAAATCCAGGGTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACAGCCCATGCATTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGATTAGTACCCTTAATATTAGGGGCACCTGATATAGCTTTCCCACGAATAAATAATATAAGATTTTGGATATTACCCCCTTCTTTAACTCTTCTTATCTCAAGAAGAATTGTAGAAAATGGGGCTGGAACTGGATGAACAGTTTATCCCCCTTTATCCTCTAATATTGCTCATAGAGGGAGATCTGTTGATTTAGCTATTTTTTCTCTCCATTTAGCTGGGATTTCCTCCATTTTAGGAGCTATTAACTTTATTACTACTATTATTAACATGCGAGTTAATAATATATCATTTGATCAAATACCATTATTTGTTTGAGCAGTTGGGATTACTGCATTTTTACTTTTACTATCTTTACCAGTTTTAGCTGGAGCTATTACCATACTATTAACTGACCGAAATTTAAATACATCTTTTTTTGATCCGGCAGGAGGAGGAGATCCTATTCTCTATCAACATTTATTT"))
colnames(dfRefSeq)[2] <- "nucleotides"
dfRefSeq$nucleotides <- as.character(dfRefSeq$nucleotides)

#We also have to convert to type character for this to work and add all of the sequences plus the 
#reference into a vector, reference sequence is added as the first sequence
alignmentSequences <- as.character(dfAllSeq$nucleotides)
#Also name our sequences according to bin
names(alignmentSequences) <- dfAllSeq$bin_uri
alignmentRef <- as.character(dfRefSeq$nucleotides[1])
#Name our reference as reference
names(alignmentRef) <- "reference"
#Append our sequences together
alignmentSequencesPlusRef <- append(alignmentRef, alignmentSequences)

#Converting all sequences in dfAllSeq plus reference to DNAStringSet format, this is the format required for the alignment
dnaStringSet2 <- DNAStringSet(alignmentSequencesPlusRef)

#Run a multiple sequence alignment of all sequences
#Using default settings of ClustalOmega of the msa package to speed up the alignment
#This could take several minutes to even hours in some cases depending on the taxa
alignment2 <- msaClustalOmega(dnaStringSet2)

##############
#Sequence Trimming according to the Reference Sequence

#For trimming of the sequences we have to determine where in the alignment the reference sequence is and 
#determine its start and stop positions relative to the other sequences
#we can then use these positions to trim the rest of the sequences in the alignment
refSeqPos <- which(alignment2@unmasked@ranges@NAMES == "reference")
refSeqPos <- alignment2@unmasked[refSeqPos]

#Finding start position by searching for the first nucleotide position of the reference sequence
refSeqPosStart <- regexpr("[ACTG]", refSeqPos)
refSeqPosStart <- as.numeric(refSeqPosStart)

#Finding last nucleotide position of the reference sequence, the regex may need to be changed depending on the reference sequence
#I use the last 20 bp taken from the reference to identify where the sequence ends in the alignment
last20bp <- substr(dfRefSeq$nucleotides[1], (nchar(dfRefSeq$nucleotides[1])-20), nchar(dfRefSeq$nucleotides[1]))
refSeqPosEnd <- regexpr(last20bp, refSeqPos)
refSeqPosEnd <- as.numeric(refSeqPosEnd)

#Then we can substr the alignment by these positions to effectively trim the alignment
alignment2Trim <- substr(alignment2, refSeqPosStart, refSeqPosEnd)

#Again convert to dnaStringSet format
dnaStringSet3 <- DNAStringSet(alignment2Trim)

#Remove our reference sequence from this as we dont want this to be included in further analysis
refSeqRemove <- which(dnaStringSet3@ranges@NAMES == "reference")
dnaStringSet3 <- subset(dnaStringSet3[-refSeqRemove])

#Reorganization of the AllSeq dataframe

#Also reordering dfAllSeq according to the ordering of bins produced in the alignment (now contained in the dnaStringSet object)
#First make a variable with the ordering
alignmentOrder <- dnaStringSet3@ranges@NAMES
#Then order dfAllSeq according to this
dfAllSeq <- dfAllSeq[match(alignmentOrder, dfAllSeq$bin_uri),]
#Then renumber indices of dfAllSeq accordingly
dfAllSeq$ind <- 1:nrow(dfAllSeq)
#Have to change to character as well
dfAllSeq$ind <- as.character(dfAllSeq$ind)
#Also will repopulate dfAllSeq with the newly trimmed sequences instead of the raw sequences in dfInitial
trimmedSeq <- as.character(dnaStringSet3)
dfAllSeq$nucleotides <- trimmedSeq

#Remnaming dnastringset3 by dfAllSeq index (important later on for dataframe merging)
dnaStringSet3@ranges@NAMES = dfAllSeq$ind

####################
#Pairwise Distance Determination with TN93

#Conversion to DNAbin format before using genetic distance matrix
DNAbin <- as.DNAbin(dnaStringSet3)

#Now for the computation of genetic distance, several models can be used - "raw", "N", "TS", "TV", "JC69", "K80" (the default), 
#"F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#Details on each model can be found here: http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Using the TN93 model for our data
matrixGeneticDistance <- dist.dna(DNAbin, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
#convert to dataframe
dfGeneticDistance <-as.data.frame(matrixGeneticDistance)
#Putting it into a stack (each column concatenated into one long column of indexes and values) so it can be easily subsetted
dfGeneticDistanceStack <-stack(dfGeneticDistance)

#############
#Latitude Distance Determination

#Can be determined easily with the dist function which will determine latitudinal differences between all bins
distLat <- dist(dfAllSeq$medianLat)
#Then we convert this to a matrix
matrixLatitudeDistance <- as.matrix( dist(dfAllSeq$medianLat) )
#Then convert to dataframe since its easier to manipulate, this dataframe will be used further down 
#(in testing I found matrices tend to be less reliable for further manipulations)
dfLatitudeDistance <-as.data.frame(matrixLatitudeDistance)

##############
#Finding appropriate pairings according to latitude and distance criteria

#So now we have two dataframes, we can find ideal matchings based on less than 15% divergence and 20 degrees latitude separation
#These values can easily be edited to add more or less stringency to the matches
#Will produce lists with indexes of each match according to our set criteria
geneticDistanceMatchI <- which(dfGeneticDistance<=0.15)
latitudeDistanceMatch <- which(dfLatitudeDistance>=20)

#Match the dfs against each other to find the sister pairs using intersect
matchOverall <- intersect(geneticDistanceMatchI,latitudeDistanceMatch)

#Use these index matches and reference against the genetic distance stack dataframe to subset it according to matches
dfMatchOverall <- dfGeneticDistanceStack[c(matchOverall), ]
#Each duplication of genetic distance values corresponds to a pairing 
#Will now merge the Allseq dataframe to the pairings generated 
dfMatchOverall <- merge(dfAllSeq, dfMatchOverall, by.x = "ind", by.y = "ind")

#Pairings will be ordered according to genetic distance, least divergent pairing to most divergent pairing
dfMatchOverall <- dfMatchOverall[order(dfMatchOverall$values),] 

#Then we can multiply these distance values by 1.3 to determine the minimum outgroup distance from the pairings, 
#we put this in another column in the matchOverall dataframe
#This minimum outgroup distance could also be a user adjustable parameter and can be easily modified
dfMatchOverall$inGroupDistx1.3 <- dfMatchOverall$values * 1.3

#Also grouping dataframe every 2 rows to reflect each unique pairing/matching, pairing column will give a 
#number value for each pairing ordered 
dfMatchOverall$inGroupPairing <- rep(1:(nrow(dfMatchOverall)/2), each = 2)

#Reorganizing and renaming some columns in MatchOverall dataframe to make more easily readable 
dfMatchOverall <- (dfMatchOverall[,c("inGroupPairing","record_id","bin_uri","values","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","ind")])
colnames(dfMatchOverall)[4] <- "inGroupDist"
colnames(dfMatchOverall)[26] <- "indexNo"

##############
#Taking the Best Possible Pairings (no bin duplicated in any pairing) and Creating Dataframes for each 
#Pairing Lineage and each Complete Pairing (both lineages)

#For bins that have multiple pairings associated with them we can select which pairing has the smallest divergence 
#to its other ingroup member since the pairings are ordered according to distance
#Basically we are taking the best possible pairing per bin based on divergence such that a bin would never be 
#duplicated in multiple pairings
#This dataframe will be dfMatchOverallLineage1
dfMatchOverallLineage1 <- by(dfMatchOverall, dfMatchOverall["bin_uri"], head, n=1)
#This dataframe represents the first lineage of each sister pairing
dfMatchOverallLineage1 <- Reduce(rbind, dfMatchOverallLineage1)
dfMatchOverallLineage1 <- subset(dfMatchOverallLineage1, duplicated(dfMatchOverallLineage1$inGroupPairing))

#To get both lineages for each pairing we do another subset against dfMatchOverall, we can call this dfMatchOverallBest
dfMatchOverallBest <- subset(dfMatchOverall, dfMatchOverall$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)

#To get the second lineage of each sister pairing we can subtract dfMatchOverallLineage1 by dfMatchOverall to 
#get the bin uri's for the alternative lineage
dfMatchOverallLineage2 <- setdiff(dfMatchOverallBest$bin_uri, dfMatchOverallLineage1$bin_uri)
dfMatchOverallLineage2 <- as.data.frame(dfMatchOverallLineage2)
colnames(dfMatchOverallLineage2)[1] <- "bin_uri"
#merge with match overallbest to get all the relevant data for the second lineage
dfMatchOverallLineage2 <- merge(dfMatchOverallLineage2, dfMatchOverallBest, all.x = TRUE)
#Reorganization of columns for second lineage
dfMatchOverallLineage2 <- (dfMatchOverallLineage2[,c("inGroupPairing","record_id","bin_uri","inGroupDist","inGroupDistx1.3","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo")])

##############
#Eliminating Pairings based on Overlapping Latitudinal Range

#Next we can work on establishing latitudinal ranges for each pairing
#If the two lineages of a pairing have overlapping latitude regions of greater than 25% then we would not 
#consider that pairing as a viable pairing
#This is because we want each lineage of a pairing to meet an appropriate difference in latitude

#First lets order our lineages
dfMatchOverallLineage1 <- dfMatchOverallLineage1[order(dfMatchOverallLineage1$inGroupPairing),] 
dfMatchOverallLineage2 <- dfMatchOverallLineage2[order(dfMatchOverallLineage2$inGroupPairing),] 

#we can define the overlap range threshold as 25% of the latitude range of L1
#If an overlap is greater than this value we would discard with this pairing
#Of course this value could be easily modified to add more or less stringency to the script
rangeThreshold <- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% ((dfMatchOverallLineage1$latMax[l] - dfMatchOverallLineage1$latMin[l]) * 0.25) 

#Define our latitude ranges for each lineage of a pairing
rangeL1 <- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% range(dfMatchOverallLineage1$latMax[l], dfMatchOverallLineage1$latMin[l])
rangeL2 <- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% range(dfMatchOverallLineage2$latMax[l], dfMatchOverallLineage2$latMin[l])

#Then we can determine the overlap region between them using the Overlap function from the Desctools package
#Overlap will return an absolute value so we dont have to worry about negatives for overlap values
overlapValue <- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% Overlap(rangeL1[[l]], rangeL2[[l]])

#Then if there is a range overlap between two lineages in a pairing, we can determine if this overlap is actually larger 
#than the 25% value of rangeThreshold for each individual pairing
rangeOverlapCheck <- foreach(l=1:nrow(dfMatchOverallLineage1)) %do% which(rangeThreshold[[l]]<overlapValue[[l]])

#Then overlaps values exceeding that 25% value we set should be returned as an integer of 1, if an overlap does 
#not exceed this value, then it will return a value of 0
#If there is a value that meets these criteria, its associated pairing will be removed from the MatchOverallBest 
#and lineage dataframes (we will retain it in dfMatchOverall for reference)
#This will name each element of rangeOverlapCheck with the inGroupPairing number to identify the pairing we need to eliminate
names(rangeOverlapCheck) <- paste0(dfMatchOverallLineage1$inGroupPairing)
#Identify which pairing in rangeOverlapCheck is greater than 0
overlapInd <- which(rangeOverlapCheck>0)
#if overlapInd is not empty:
if(length(overlapInd)>0){
  #Eliminate based on that pairing number(s) for MatchOverallLineage1
  dfMatchOverallLineage1 <- dfMatchOverallLineage1[-overlapInd,]
  #Now subset to both dfMatchOverallBest and dfMatchOverallLineage2
  dfMatchOverallBest <- subset(dfMatchOverallBest, dfMatchOverallBest$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)
  dfMatchOverallLineage2 <- subset(dfMatchOverallLineage2, dfMatchOverallLineage2$inGroupPairing %in% dfMatchOverallLineage1$inGroupPairing)
}

################
#Outgroup determination for each Pairing

#Now we can search for the best possible outgroupings for each pairing
#This will involve searching for outgroups relative to each pairing and finding one that is far enough away from both lineages

#First we can search for outgroupings relative to lineage 1
#We can use the indexNo to subset the dfGeneticDistanceStackDataframe according to indexes represented 
#in lineage 1, this will be called dfBestOutGroupL1
#This will essentially limit our outgroup distances to those associated with lineage1
dfBestOutGroupL1 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallLineage1$indexNo)
#We also order by Lineage1
dfBestOutGroupL1 <- dfBestOutGroupL1[order(match(dfBestOutGroupL1[,2],dfMatchOverallLineage1[,26])),]

#Then we can find which indices match lineage1 with bestoutgroup l1
outGroupCandidatesL1a <- foreach(i=1:nrow(dfMatchOverallLineage1)) %do% which(dfBestOutGroupL1$ind == dfMatchOverallLineage1$indexNo[i])
#Determining indices with correct outgroup distance, setting a small range between minimum outgroup distance and +0.1 distance 
#to narrow results, range could be modified if outgroups cannot be found
outGroupCandidatesL1b <- foreach(i=1:nrow(dfMatchOverallLineage1)) %do% which(dfBestOutGroupL1$values >= dfMatchOverallLineage1$inGroupDistx1.3[i] & dfBestOutGroupL1$values < dfMatchOverallLineage1$inGroupDistx1.3[i]+0.15)
#Intersection of the two using mapply to find the correct outgroupings for each pairing
outGroupCandidatesL1c <- mapply(intersect,outGroupCandidatesL1a,outGroupCandidatesL1b)
#Unlist to make into one vector 
outGroupCandidatesL1c <- unlist(outGroupCandidatesL1c)
#Adding an rownum column to dfBestOutGroup, this represents the second index of each outgroup candidate
dfBestOutGroupL1$rownum<-seq.int(nrow(dfBestOutGroupL1))
#Then we can subset to rownum column based on the outgroup candidates
dfBestOutGroupL1 <- dfBestOutGroupL1[dfBestOutGroupL1$rownum %in% outGroupCandidatesL1c, ]

#Now for the second lineage, we go through the same process
dfBestOutGroupL2 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% dfMatchOverallLineage2$indexNo)
dfBestOutGroupL2 <- dfBestOutGroupL2[order(match(dfBestOutGroupL2[,2],dfMatchOverallLineage2[,26])),]
outGroupCandidatesL2a <- foreach(i=1:nrow(dfMatchOverallLineage2)) %do% which(dfBestOutGroupL2$ind == dfMatchOverallLineage2$indexNo[i])
outGroupCandidatesL2b <- foreach(i=1:nrow(dfMatchOverallLineage2)) %do% which(dfBestOutGroupL2$values >= dfMatchOverallLineage2$inGroupDistx1.3[i] & dfBestOutGroupL2$values < dfMatchOverallLineage2$inGroupDistx1.3[i]+0.15)
outGroupCandidatesL2c <- mapply(intersect,outGroupCandidatesL2a,outGroupCandidatesL2b)
outGroupCandidatesL2c <- unlist(outGroupCandidatesL2c)
dfBestOutGroupL2$rownum<-seq.int(nrow(dfBestOutGroupL2))
dfBestOutGroupL2 <- dfBestOutGroupL2[dfBestOutGroupL2$rownum %in% outGroupCandidatesL2c, ]

#Now we find the intersection between outGroupListL1 and outGroupListL2, this will determine which outgroup
#candidates are shared between both lineages
dfBestOutGroupL1 <- dfBestOutGroupL1[dfBestOutGroupL1$rownum %in% dfBestOutGroupL2$rownum, ]
dfBestOutGroupL2 <- dfBestOutGroupL2[dfBestOutGroupL2$rownum %in% dfBestOutGroupL1$rownum, ]

#Taking the first entry of each index
dfBestOutGroupL1 <- by(dfBestOutGroupL1, dfBestOutGroupL1["ind"], head, n=1)
dfBestOutGroupL1 <- Reduce(rbind, dfBestOutGroupL1)
#Making sure that dfBestOutGroupOverall is ordered by rownum, this will be important for the next step
dfBestOutGroupL1 <- dfBestOutGroupL1[order(dfBestOutGroupL1$rownum),] 
#Lineage2
dfBestOutGroupL2 <- by(dfBestOutGroupL2, dfBestOutGroupL2["ind"], head, n=1)
dfBestOutGroupL2 <- Reduce(rbind, dfBestOutGroupL2)
dfBestOutGroupL2 <- dfBestOutGroupL2[order(dfBestOutGroupL2$rownum),] 

#this will give the correct index value for the rownum column, doing this for both L1 and L2
j=0:(nrow(dfBestOutGroupL1)-1)
dfBestOutGroupL1$indexNo <- (dfBestOutGroupL1$rownum / (nrow(dfAllSeq))-j) * (nrow(dfAllSeq))
#for merging, need to change ind column in dfAllSeq to numeric first
indNum <- with(dfAllSeq, as.numeric(as.character(ind)))
dfAllSeq$ind <- indNum
#using data.tables and setting keys for the merge
dfAllSeq <- data.table(dfAllSeq)
dfBestOutGroupL1 <- data.table(dfBestOutGroupL1)
setkey(dfBestOutGroupL1,indexNo)
setkey(dfAllSeq,ind)
dfBestOutGroupL1 <- merge(dfBestOutGroupL1, dfAllSeq, by.x = "indexNo", by.y = "ind", all.x = TRUE)
#Converting back to dataframe for further use
dfBestOutGroupL1 <- as.data.frame(dfBestOutGroupL1)
dfBestOutGroupL1 <- dfBestOutGroupL1[order(match(dfBestOutGroupL1[,3],dfMatchOverallLineage1[,26])),]

#Lineage2
j=0:(nrow(dfBestOutGroupL2)-1)
dfBestOutGroupL2$indexNo <- (dfBestOutGroupL2$rownum / (nrow(dfAllSeq))-j) * (nrow(dfAllSeq))
dfBestOutGroupL2 <- data.table(dfBestOutGroupL2)
setkey(dfBestOutGroupL2,indexNo)
dfBestOutGroupL2 <- merge(dfBestOutGroupL2, dfAllSeq, by.x = "indexNo", by.y = "ind", all.x = TRUE)
dfBestOutGroupL2 <- as.data.frame(dfBestOutGroupL2)
dfBestOutGroupL2 <- dfBestOutGroupL2[order(match(dfBestOutGroupL2[,3],dfMatchOverallLineage2[,26])),]

#Can rename certain columns to more closely resemble dfMatchOverallLineage dataframes 
#(except the outGroupDist column which would be unique to each lineage of course)
dfBestOutGroupL1 <- dfBestOutGroupL1[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo","ind")]
colnames(dfBestOutGroupL1)[3] <- "outGroupDist"
dfBestOutGroupL2 <- dfBestOutGroupL2[,c("bin_uri","record_id","values","medianLat","latMin","latMax","binSize","phylum_taxID","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name","species_taxID","species_name","nucleotides","medianLon","indexNo","ind")]
colnames(dfBestOutGroupL2)[3] <- "outGroupDist"

#adding an ingroup pairing column to each bestoutgroup dataframes
dfBestOutGroupL1$inGroupPairing <- dfMatchOverallLineage1$inGroupPairing
dfBestOutGroupL2$inGroupPairing <- dfMatchOverallLineage2$inGroupPairing

#Can now merge our outgroups to our associated lineages
#Each outgroup will be duplicated for each member of the pairing however the distances should be unique 
#to represent each distance from each outgroup to each lineage

#.x beside headings are the ingroup pairing data (except unique columns), scroll far enough to the right 
#and headings with .y are the outgroup related columns
dfMatchOverallLineage1 <- merge(dfMatchOverallLineage1, dfBestOutGroupL1, by.x = "inGroupPairing", by.y = "inGroupPairing")
dfMatchOverallLineage2 <- merge(dfMatchOverallLineage2, dfBestOutGroupL2, by.x = "inGroupPairing", by.y = "inGroupPairing")

#Then we can revise dfMatchOverallBest to correctly reflect both pairings with correct outgroup distances
#We do this by using rbind to combine both lineages together into the dfMatchOverallBest dataframe
dfMatchOverallBest <-  rbind(dfMatchOverallLineage1, dfMatchOverallLineage2)
#Then order by inGroupPairing once again for a better organization of the dataframe
dfMatchOverallBest <- dfMatchOverallBest[order(dfMatchOverallBest$inGroupPairing),]

#As a last step for outgroup determination, if a suitable outgroup for both lineages could not be found and 
#is NA for a pairing, then we can filter these pairings out
#If a suitable outgroup for a pairing cannot be found, I suggest going back and adjusting the criteria for the outgroupings
noOutGroup<-which(is.na(dfMatchOverallBest$bin_uri.y))
if(length(noOutGroup) >0){
  dfMatchOverallBest<-dfMatchOverallBest[-noOutGroup,]
}
noOutGroup<-which(is.na(dfMatchOverallLineage1$bin_uri.y))
if(length(noOutGroup) >0){
  dfMatchOverallLineage1<-dfMatchOverallLineage1[-noOutGroup,]
}
noOutGroup<-which(is.na(dfMatchOverallLineage2$bin_uri.y))
if(length(noOutGroup) >0){
  dfMatchOverallLineage2<-dfMatchOverallLineage2[-noOutGroup,]
}

#Also renumber pairings starting from 1 again
dfMatchOverallBest$inGroupPairing <- rep(1:(nrow(dfMatchOverallBest)/2), each = 2)
dfMatchOverallLineage1$inGroupPairing <- 1:nrow(dfMatchOverallLineage1)
dfMatchOverallLineage2$inGroupPairing <- 1:nrow(dfMatchOverallLineage2)

################
#Ouput to CSV or TSV

#The user can output the pairing results to CSV/TSV if they wish

#First create this dataframe for dfMatchOverallBest, this will allow output of this file
#dfMatchOverallBestOut <- data.frame(lapply(dfMatchOverallBest, as.character), stringsAsFactors=FALSE)

#Defining another variable to give a unique name to the CSV
#Name would be what you specify it to be, R will prompt you to insert a name for the file in the console:
#filename <- readline(prompt="")

#Then you can uncomment one of these write.table commands to output

#CSV
#write.table(dfMatchOverallBestOut, file=paste(filename, ".csv", sep=""), quote=FALSE, sep=',', col.names = NA)

#TSV
#write.table(dfMatchOverallBestOut, file=paste(filename, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)
