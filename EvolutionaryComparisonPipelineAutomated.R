#################
#Authored primarily by Matthew Orton
#Contribution by Jacqueline May for lines 220-246 of script 
#Contributions by David Lee for lines 1148-1169, editing and formatting of script
#Contributions by Winfield Ly for editing and formatting of script
#Collaboration of Sarah Adamowicz (University of Guelph) in designing the analyses


#Pipeline Purpose (Class-Level Analyses):

#This pipeline will allow for the generation of latitudinally separated sister pairings and
#associated outgroupings from nearly any taxa (provided they have a suitable reference sequence and are not too large)
#and geographical region found on the BOLD API in one streamlined R pipeline!
#In this new and improved iteration, data are translated directly from a BOLD tsv file of
#the user's choosing to a dataframe in R. The generated sister pairs and outgroups can then
#be written to a csv or tsv, and the file will appear in the current working directory of R.
#Additionally, binomial and Wilcoxon tests can be performed on the relative outgroup 
#distances of the generated pairings, and a plot of the resultant relative outgroup distances can be 
#generated for each order.
#A world map visualizing the latitudinally separated pairings for each taxon can also be run using plotly.

#Larger taxa including various phyla can be run and will get broken down into classes, and sister pairs are sought
#within classes by default. 
#Lowever-level taxonomic groups (below the class level) can be run. However, the analysis will still be performed at the class level. (Matt a clarification is needed here. I do not follow the explanation on this line, as currently phrased.)

##################
#Some important tips:

#Make sure to use RStudio since this will provide a better interface to run the script than just R.

#The taxa being tested must use BIN identifiers as a means of identification; 
#for instance, plants on the BOLD API cannot be run using this pipeline since they do not use BIN identifiers.

#At this time, it is best not to run Arthropoda, Insecta, or Chordata since these are all very large and computationally demanding 
#and should be broken down into smaller sub-taxa.

#Larger taxa consume a great deal of working memory, so be mindful of the amount of available working memory you have available on your 
#PC/Mac.
#As a benchmark, the class Actinopterygii which contains roughly 250000 unique entries on BOLD consumes close to 3.4 Gb of working memory
#to run.
#It's probably a good idea to ensure RStudio is the only major application running when you run the script

#Entries on BOLD which do not at least have order-level classification will be filtered out (Matt - should this read "class-level" in this sentence, for this version of the pipeline?)
#since the pipeline cannot properly categorize these. 

#You do not have to worry about entries missing certain pieces of data, for example latitude coordinates, sequence data etc.
#The script will filter these records out automatically.

#There are two options for parsing the tsv. You can either download the tsv directly from the
#BOLD API OR you can use a tsv you have previously downloaded. Don't forget you can also
#save your workspace once a TSV is downloaded so you dont have to download it again.

#When testing a taxon for the first time, you will have to ensure that you have a suitable 
#reference sequence for it and ensure that it is inserted into the dfRefSeq dataframe.
#The reference sequence is a high-quality sequence from a given taxon, which is used as the basis for the first alignment step.
#As well, the reference sequence is used to trim all sequences to the same length for analysis.
#To insert a sequence to the dfRefSeq dataframe, simply add another taxon and sequence in quotations in the refseq dataframe command
#in the reference sequence section.
#Make sure that you keep the entire sequence in one line; breaking it up on to separate lines will add a new line character.

#Some tips for using the BOLD API since this is what is used to grab the relevant data we need: 
#To see details on how to use the bold API, 
#go to http://www.boldsystems.org/index.php/resources/api?type=webservices
#We can add additional restrictions to the url, for example 
#&instituiton=Biodiversity Institute of Ontario|York University or 
#&marker=COI-5P if you want to specifiy an institution or specific genetic marker in the url.
#geo=all in the url means global but geo can be geo=Canada, for example.
#Can use | modifier in the url, for example &geo=Canada|Alaska would give data for 
#both Canada and Alaska, 
#or taxon=Aves|Reptilia would yield give data for both Aves and Reptilia.

#################
#Important dataframes:

#dfPVal gives the p-values for both the binomial test and wilcoxon tests for all classes (that have pairings)
#and each order separately.

#dfPairingResultsSummary shows a more user-friendly summation of the pairing results.

#dfRelativeDist shows the relative distances to the outgroup for each pairing, including pseudoreplicates.

#dfPairingResultsL1L2 is the more detailed finalized dataframe that contains all of the finalized pairings
#and outgroupings.

#dfPairingResultsL1 and dfPairingResultsL2 represent dataframes for each lineage of each pairing respectively.

#dfInitial is the dataframe first produced by the import from BOLD and is filtered by lat, 
#bin_uri, N content, Gap content and sequence length.

#dfOutGroupL1 and L2 contain the associated outgroupings only (for each lineage), but 
#each one does have a column indicating the pairing with which it is associated. (Matt - please verify this edit.)

#dfCentroid contains centroid sequences for all bins with more than one member. 

#dfNonCentroid contains all bins that only have one member and thus do not need centroid sequences.

#dfAllSeq is simply dfNonCentroid and dfCentroid combined together in one dataframe.

#dfGeneticDistanceStack is all genetic distance matrices concatenated into one long 
#column of values. It is used to grab bin_uri's and distances for each pairing and assist in determining outgroup 
#bin_uri's and distances.

#dfRefSeq shows various taxa with a suitable reference sequence that has been found for them.
#This dataframe may be broken up into dfRefSeqOrder and dfRefSeqLargeOrder if there are very large orders that need to broken
#down to families.

#dfPseudoRepAverage shows the averaged distances for each set of pseudoreplicates in dfPseudoRep.

##################
#Important Variables

#alignment2 will show the alignment of class
#(ex: typing alignment2[1] will show the first alignment performed)
#It is important to always check this variable to ensure the alignment is producing a 
#reasonable result

#alignment2Trim will show the trimmed sequences after alignment 
#with the reference and trimming with the reference
#Again, good to check this to make sure the result is reasonable

#binomialTestOutGroup shows the results of the binomial test

#wilcoxonTestOutgroup shows the results of the wilcoxon test

#mapLayout will allow for customization of the world map for map plotting using Plotly

#################
#Packages required
#Note that once you have installed the packages (first time running the script), 
#you only have to run the libraries again each time you open up RStudio

#We need the foreach package for several functions that require iteration over dataframe rows
#install.packages("foreach")
library(foreach)
#For genetic distance determination using the TN93 model, we use the ape package
#install.packages("ape")
library(ape)
#Speeds up parsing of the tsv file with read_tsv function
#install.packages("readr")
library(readr)
#For sequence alignments we need the biostrings (DNAStringSet function) and msa packages, 
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library("Biostrings")
#biocLite("msa")
library("msa")
#For overlapping latitude regions we need the Desctools package
#install.packages("DescTools")
library(DescTools)
#Also adding data tables for table merging in the outgrouping section
#install.packages("data.table")
library(data.table)
#For plotting of relative outgroup distances between lineages we will also need ggplot2
require(ggplot2)
#plotly for map plotting functionality
#install.packages("plotly")
library(plotly)

#################
#R Commands:

#TSV Parsing

#First we download the TSV and convert it into a dataframe, the URL below is what is modified 
#by the user and will determine the taxa, geographic region, etc.
#ex: taxon=Aves$geo=all

#read_tsv has been modified to select only certain columns to save on downloading time 
dfInitial <- read_tsv(
  "http://www.boldsystems.org/index.php/API_Public/combined?taxon=Annelida&geo=all&format=tsv")[ ,
                                              c('recordID', 'bin_uri','phylum_taxID','phylum_name','class_taxID',
                                                'class_name','order_taxID','order_name','family_taxID','family_name',
                                                'subfamily_taxID','subfamily_name','genus_taxID','genus_name',
                                                'species_taxID','species_name','lat','lon','nucleotides')]

#If you want to run pre downloaded BOLD TSV's to avoid downloading of the same tsv multiple times, 
#this will let you choose a path to that TSV and parse:

#tsvParseDoc <- file.choose()
#dfInitial <- read_tsv(tsvParseDoc)[ ,
#                          c('recordID', 'bin_uri','phylum_taxID','phylum_name','class_taxID',
#                          'class_name','order_taxID','order_name','family_taxID','family_name',
#                          'subfamily_taxID','subfamily_name','genus_taxID','genus_name',
#                          'species_taxID','species_name','lat','lon','nucleotides')]

#***Keep in mind you can also save your R workspace to avoid redownloading TSV's***

##############
#Dataframe Filtering and Reorganization

colnames(dfInitial)[1] <- "record_id"

#Removing sequences with no latitude values, filtering according to lat since we only really 
#need lat for the analysis
containLat <- grep( "[0-9]", dfInitial$lat)
dfInitial<-dfInitial[containLat,]

#Next we have to convert the lat column to num instead of chr type, this will become important 
#later on for median latitude determination
latNum <- with(dfInitial, as.numeric(as.character(lat))) 
dfInitial$latNum <- latNum

#Can do lon as well for map plotting using plotly
lonNum <- with(dfInitial, as.numeric(as.character(lon))) 
dfInitial$lonNum <- lonNum

#First identifying missing bins and eliminating rows with missing bin_uri's since bin is 
#a big indicator of sequence quality
#Grep by colon since every record with a bin identifier will have this
containBin <- grep( "[:]", dfInitial$bin_uri)
dfInitial <- dfInitial[containBin,]

#Getting rid of any records that dont have sequence data (sometimes there are a few)
containNucleotides <- grep( "[ACGT]", dfInitial$nucleotides)
dfInitial <- dfInitial[containNucleotides,]

#Gap content and N content will affect the Clustal Omega alignment and the alignment will give warning messages
#so we need to filter out sequences with high gap and N content

# First, we need to convert nucleotides to chr type.
dfInitial$nucleotides <- with(dfInitial, (as.character(nucleotides))) 
# Cut off starting Ns and gaps (large portions of Ns and gaps at the start of a sequence).
start_N_gap <- sapply(regmatches(dfInitial$nucleotides, gregexpr("^[-N]", dfInitial$nucleotides)), length)
start_N_gap <- foreach(i=1:nrow(dfInitial)) %do%
  if (start_N_gap[[i]] > 0) {
    split <- strsplit(dfInitial$nucleotides[i], "^[-N]+")
    dfInitial$nucleotides[i] <- split[[1]][2]
  }
# Cut off ending Ns and gaps (large portions of Ns and gaps at the end of a sequence).
end_N_gap <- sapply(regmatches(dfInitial$nucleotides, gregexpr("[-N]$", dfInitial$nucleotides)), length)
end_N_gap <- foreach(i=1:nrow(dfInitial)) %do%
  if (end_N_gap[[i]] > 0) {
    split <- strsplit(dfInitial$nucleotides[i], "[-N]+$")
    dfInitial$nucleotides[i] <- split[[1]][1]
  }
# This will give the number of positions where an *internal* N or gap is found for each sequence.
internal_N_gap <- sapply(regmatches(dfInitial$nucleotides, gregexpr("[-N]", dfInitial$nucleotides)), length)
# We then go loop through each sequence to see if the number of Ns or gaps is greater than 1% (0.01) of
# the total sequence length.
internal_N_gap <- foreach(i=1:nrow(dfInitial)) %do% 
  which((internal_N_gap[[i]]/nchar(dfInitial$nucleotides[i])>0.01))
N_gap_check <- sapply( internal_N_gap , function (x) length( x ) )
N_gap_check <- which(N_gap_check>0)
# Subset out these higher gap and N content sequences.
dfInitial <- dfInitial[-N_gap_check,]

#Filter out sequences less than 620 bp and greater than 1000 bp since these sequence length extremes can interfere with the alignment
#and this also helps to standardize sequence length for the reference sequences.
sequenceLengths <- nchar(dfInitial$nucleotides)
sequenceLengthCheck <- which(sequenceLengths>1000 | sequenceLengths<620)
dfInitial <- dfInitial[-sequenceLengthCheck,]

#Modifying Bin column slightly to remove "BIN:"
dfInitial$bin_uri <- substr(dfInitial$bin_uri, 6 , 13)

#Dataframe Reorganization
dfInitial <- (dfInitial[,c("record_id","bin_uri","phylum_taxID","phylum_name","class_taxID",
                           "class_name","order_taxID","order_name","family_taxID",
                           "family_name","subfamily_taxID","subfamily_name","genus_taxID",
                           "genus_name","species_taxID","species_name","nucleotides",
                           "latNum","lonNum")])
############
#Bin Stats and Median Latitude/Longitude Determination per bin

#First we can make a smaller dataframe with the columns we want for each bin - bin_uri,
#latnum, lonnum, record_id, 
dfBinList <- (dfInitial[,c("record_id","bin_uri","latNum","lonNum","nucleotides")])

#Conversion to absolute values before median latitude values are calculated on dfBinList
dfBinList$latNumAbs <- abs(dfBinList$latNum) 

#Lon remains untouched as it is not needed for the analysis

#Create groupings by bin with each grouping representing a different bin_uri
#Each element of this list represents a bin with subelements representing the various 
#columns of the initial 
#dataframe created and the information is grouped by bin 
binList <- lapply(unique(dfBinList$bin_uri), function(x) dfBinList[dfBinList$bin_uri == x,])

#Now to determine a median latitude for each bin based on absolute values
medianLatAbs <- sapply( binList , function(x) median( x$latNumAbs ) )

#median lat for mapping purposes only
medianLatMap <- sapply( binList , function(x) median( x$latNum ) )

#We also need a median longitude for map plotting
medianLon <- sapply( binList , function(x) median( x$lonNum ) )

#we can also take a few other important pieces of data regarding each bin using sapply 
#including number of record_ids to a bin and latitudinal min and max of each bin
latMin <- sapply( binList , function(x) min( x$latNum ) )
latMax <- sapply( binList , function(x) max( x$latNum ) )
binSize <- sapply( binList , function (x) length( x$record_id ) )

#Dataframe of our median lat values, this will be used in our final dataframe
dfLatLon <- data.frame(medianLatAbs)

#Adding bin_uri, latMin, latMax and binSize to dataframe with medianLat
dfLatLon$bin_uri <- c(unique(dfInitial$bin_uri))
dfLatLon$medianLon <- c(medianLon)
dfLatLon$latMin <- c(latMin)
#Convert to absolute value for latMin and Max
dfLatLon$latMin <- abs(dfLatLon$latMin)
dfLatLon$latMax <- c(latMax)
dfLatLon$latMax <- abs(dfLatLon$latMax)
dfLatLon$binSize <- c(binSize)
dfLatLon$medianLatMap <- c(medianLatMap)

#Merging LatLon to BinList for the sequence alignment step
dfBinList <- merge(dfBinList, dfLatLon, by.x = "bin_uri", by.y = "bin_uri")

#Also reordering dfLatLon by bin_uri for a step later on
dfLatLon <- dfLatLon[order(dfLatLon$bin_uri),]

###############
#Selecting a Centroid Sequence Per Bin
#Centroid Sequence: bin sequence with minimum average pairwise distance to all other bin sequences in a given bin

#First we have to subset dfBinList to find bins with more than one member since these bins will need centroid sequences
largeBin <-which(dfBinList$binSize > 1)
#If there is at least one bin with more than one member, 
#then a dataframe dfCentroid will be created with those bins
if(length(largeBin) >0){
  dfCentroid <- dfBinList[largeBin,]
  
  #Also subset dfLatLon down to number of bins in dfIntiial
  dfLatLon <- subset(dfLatLon, dfLatLon$bin_uri %in%
                       dfCentroid$bin_uri)
  
  #Also need to find the number of unique bins in dfConsensus
  binNumberCentroid <- unique(dfCentroid$bin_uri)
  binNumberCentroid <- length(binNumberCentroid)
  
  #We also have to create another separate dataframe with bins that only have one member called dfNonCentroid
  dfNonCentroid <- dfBinList[-largeBin,]
  
  #We then take the dfConsensus sequences and break it down into a list with each element being a unique bin
  largeBinList <- lapply(unique(dfCentroid$bin_uri), function(x) dfCentroid[dfCentroid$bin_uri == x,])
  
  #Extract record id from each bin
  largeBinRecordId <- sapply( largeBinList , function(x) ( x$record_id ) )
  
  #Convert all of the sequences in the largeBinList to dnaStringSet format for the alignment step
  dnaStringSet1 <- sapply( largeBinList, function(x) DNAStringSet(x$nucleotides) )
  
  #name DNAStringSet with the record ids
  for (i in seq(from=1, to=binNumberCentroid, by = 1)){
    names(dnaStringSet1[[i]]) <- largeBinRecordId[[i]]
  }
  
  #Multiple sequence alignment using msa package - 3 different algorithms can be used for this: ClustalW, ClustalOmega, MUSCLE 
  #Sequence alignment package can be found here: 
  #http://master.bioconductor.org/packages/3.2/bioc/vignettes/msa/inst/doc/msa.pdf
  #Using the default parameters of ClustalOmega, 
  #the default of Clustal Omega is known as "Gonnet", you will see this as an output in the console
  #Sometimes you may see this message during this alignment: 
  #Transfer:hhalign/hhalignment-C.h:2968: profile has no leading and/or trailing residues (h=-1:t=0:#=1)
  #Those messages simply indicates a particular sequence being aligned that 
  #is flanked by gaps hence "no leading or trailing residues"
  #Run a multiple sequence alignment on each element of the dnaStringSet1 list
  alignment1 <- foreach(i=1:binNumberCentroid) %do% msaClustalOmega(dnaStringSet1[[i]])
  
  #We can then convert each alignment to DNAbin format
  dnaBINCentroid <- foreach(i=1:binNumberCentroid) %do% as.DNAbin(alignment1[[i]])
  
  #Then performing genetic distance determination with the TN93 model on each DNAbin list
  geneticDistanceCentroid <- foreach(i=1:binNumberCentroid) %do% 
    dist.dna(dnaBINCentroid[[i]], model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
  
  #centroid sequence can be determined from the distance matrix alone, it is the sequence in 
  #a bin with minimum average pairwise distance to all other bin sequences.
  centroidSeq <- foreach(i=1:binNumberCentroid) %do% which.min(rowSums(geneticDistanceCentroid[[i]]))
  centroidSeq <- unlist(centroidSeq)
  centroidSeq <- names(centroidSeq)
  centroidSeq <- as.numeric(centroidSeq)
  
  #subset dfCentroid by the record ids on this list
  dfCentroid <- subset(dfCentroid, record_id %in% centroidSeq)
  
  #Append our nonconsensus to our consensus, will make a new dataframe for this - dfAllSeq
  #now we have all of the sequences we need for the next alignment of all sequences
  dfAllSeq <-  rbind(dfCentroid, dfNonCentroid)
  #Can then merge this with dfInitial to get all of the relevant data we need
  #merging to dfInitial to gain all the relevant taxanomic data
  dfAllSeq <- merge(dfAllSeq, dfInitial, by.x = "record_id", by.y = "record_id")
  #Renaming and reorganizing the dataframe
  dfAllSeq <- (dfAllSeq[,c("bin_uri.x","binSize","record_id","phylum_taxID","phylum_name",
                           "class_taxID","class_name","order_taxID","order_name","family_taxID",
                           "family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name",
                           "species_taxID","species_name","nucleotides.x","medianLatAbs","medianLatMap","latMin","latMax",
                           "medianLon")])
  colnames(dfAllSeq)[1] <- "bin_uri"
  colnames(dfAllSeq)[18] <- "nucleotides"
  #Adding an index column to reference later with Match overall dataframe
  dfAllSeq$ind <- row.names(dfAllSeq)
  #Deleting any possible duplicate entries
  dfAllSeq <- (by(dfAllSeq, dfAllSeq["bin_uri"], head, n=1))
  dfAllSeq <- Reduce(rbind, dfAllSeq)
  
} else {
  #Else if there are no bins with more than one member than we would simply merge latlon with initial 
  #to get dfAllSeq and thus all of the sequences we want
  dfAllSeq <- merge(dfLatLon, dfInitial, by.x = "bin_uri", by.y = "bin_uri")
  dfAllSeq <- (dfAllSeq[,c("bin_uri.x","binSize","record_id","phylum_taxID","phylum_name",
                           "class_taxID","class_name","order_taxID","order_name","family_taxID",
                           "family_name","subfamily_taxID","subfamily_name","genus_taxID","genus_name",
                           "species_taxID","species_name","nucleotides.x","medianLatAbs","medianLatMap","latMin","latMax",
                           "medianLon"
  )])
  dfAllSeq$ind <- row.names(dfAllSeq)
}

#To simplify the number of dataframes, can remove dfLatLon and dfBinList since these are redundant
rm(dfBinList) 
rm(dfLatLon)

#############
#Multiple Sequence Alignment of Sequences with Reference Sequence and Beginning of Automation by Class

#Manual input of reference sequences into dfRefSeq dataframe
dfRefSeq <- data.frame(taxa = c("Clitellata","Polychaeta"),
                       nucleotides = c("AACTCTATACTTTATTTTAGGCGTCTGAGCAGGAATAGTGGGGGCTGGTATAAGACTCCTAATTCGAATTGAGCTAAGACAGCCGGGAGCATTTCTAGGAAGGGATCAACTCTATAACACTATTGTAACTGCTCACGCATTTGTAATAATTTTCTTTCTAGTAATACCTGTATTTATTGGGGGGTTCGGTAATTGACTTCTGCCTTTAATACTTGGAGCCCCTGACATGGCATTCCCACGTCTTAACAACATAAGATTTTGACTCCTTCCCCCATCACTAATCCTTCTAGTATCCTCTGCTGCAGTAGAAAAGGGTGCGGGAACTGGATGAACTGTTTATCCACCCCTAGCAAGAAACATTGCTCATGCCGGCCCATCTGTAGACTTAGCTATTTTTTCTCTTCATTTAGCAGGTGCTTCATCAATCTTGGGTGCCATTAATTTTATTACTACTGTTATTAATATACGATGAAGAGGCTTACGACTTGAACGAATCCCATTATTCGTTTGAGCCGTACTAATTACAGTGGTCCTTCTACTCTTATCCTTACCAGTATTAGCCGGTGCAATTACTATACTACTTACCGATCGAAATCTAAATACCTCCTTCTTTGACCCTGCTGGAGGTGGGGATCCAATCCTATATCAACACTTATTC",
                                       "CACCCTCTATTTTATCTTCGGAGTCTGATCTGGTCTTCTAGGTACATCAATAAGGCTACTAATCCGAATTGAGCTTGGGCAACCCGGCTCTTTTTTAGGTAGAGACCAACTGTACAATACGGTTGTAACTGCACATGCTTTTCTTATAATTTTCTTTCTTGTCATACCCGTCTTTATTGGTGGGTTTGGAAACTGACTTGTCCCATTAATACTTGCTGCTCCAGACATGGCATTTCCCCGAATAAACAACATAAGTTTCTGACTTCTCCCTCCGGCCCTGATTCTTCTTCTGAGCTCCGCTGCAGTTGAAAAGGGCGTTGGTACAGGTTGAACAGTATACCCCCCCTTATCAAGAAACCTCGCACATGCGGGTCCATCTGTAGATTTGGCCATCTTTTCTCTCCACTTAGCCGGGATCTCATCTATTCTCGGAGCTATTAACTTTATTACTACCGTAATCAACATGCGGTCTAAAGGCCTTCGATTAGAACGAGTCCCTCTATTCGTGTGGGCAGTAAAGATTACTGCTATCCTTCTTCTTTTATCACTTCCTGTTTTAGCAGGAGCCATTACTATACTCCTAACTGACCGTAATTTAAACACCTCATTCTTCGATCCAGCAGGGGGAGGAGACCCCATTCTTTACCAACACCTCTTT"))
colnames(dfRefSeq)[2] <- "nucleotides"
dfRefSeq$nucleotides <- as.character(dfRefSeq$nucleotides)

#Symmetrical trimming of the references to a standard 620 bp from 658 bp
dfRefSeq$nucleotides <- substr(dfRefSeq$nucleotides, 20, nchar(dfRefSeq$nucleotides)-19)

#Check of sequence length
dfRefSeq$seqLength <- nchar(dfRefSeq$nucleotides)

#Subset dfAllSeq by entries in the reference sequence dataframe
dfAllSeq <- subset(dfAllSeq, dfAllSeq$class_name %in% dfRefSeq$taxa)
    
#Break dfAllSeq down into the various classes, could be edited for orders instead
taxaListComplete <- lapply(unique(dfAllSeq$class_taxID), function(x) dfAllSeq[dfAllSeq$class_taxID == x,])
    
#Revise dfRefSeq dataframe to reflect this
classList <- sapply( taxaListComplete , function(x) unique( x$class_name ) )
#This command will ensure the reference sequence dataframe is in the same order as the allseq dataframe
dfRefSeq <- dfRefSeq[match(classList, dfRefSeq$taxa),]
    
#Extract sequences and bin_uri from each class
classBin <- sapply( taxaListComplete , function(x) ( x$bin_uri ) )
classSequences <- sapply( taxaListComplete, function(x) ( x$nucleotides ) )
classSequencesNames <- classBin
    
#Take our reference sequences
alignmentRef <- as.character(dfRefSeq$nucleotides)
dfRefSeq$reference <- "reference"
#Name our reference as reference for each class so it can be identified as such in the alignment
alignmentRefNames <- dfRefSeq$reference
    
#Merge our reference sequences with each of our class sequences
alignmentSequencesPlusRef <- mapply(c, classSequences, alignmentRef)
#Merge the names together
alignmentNames <- mapply(c, classSequencesNames, alignmentRefNames)

#Converting all sequences in dfAllSeq plus reference to DNAStringSet format, this is 
#the format required for the alignment
dnaStringSet2 <- sapply( alignmentSequencesPlusRef, function(x) DNAStringSet( x ) )

#Name the DNAStringSet List with the appropriate bin uri's
for (i in seq(from=1, to=nrow(dfRefSeq), by = 1)){
  names(dnaStringSet2[[i]]) <- alignmentNames[[i]]
}

#Run a multiple sequence alignment of all sequences including the reference
#Using default settings of ClustalOmega of the msa package to speed up the alignment, 
#again default scoring matrix for Clustal Omega is known as "Gonnet"
#Sometimes you may see this message during this alignment: 
#Transfer:hhalign/hhalignment-C.h:2968: profile has no leading and/or trailing residues (h=-1:t=0:#=1)
#Those messages simply indicate a particular sequence being aligned that is flanked by gaps
#This could take several minutes depending on the taxa
alignment2 <- foreach(i=1:nrow(dfRefSeq)) %do% msaClustalOmega(dnaStringSet2[[i]])

##############
#Sequence Trimming according to the Reference Sequence

#For trimming of the sequences we have to determine where in the alignment the 
#reference sequence is and determine its start and stop positions relative to the
#other sequences we can then use these positions to trim the rest of the sequences 
#in the alignment
refSeqPos <- foreach(i=1:nrow(dfRefSeq)) %do% which(alignment2[[i]]@unmasked@ranges@NAMES == "reference")
refSeqPos <- foreach(i=1:nrow(dfRefSeq)) %do% alignment2[[i]]@unmasked[refSeqPos[[i]]]

#Finding start position by searching for the first nucleotide position of the 
#reference sequence
refSeqPosStart <- foreach(i=1:nrow(dfRefSeq)) %do% regexpr("[ACTG]", refSeqPos[[i]])
refSeqPosStart <- as.numeric(refSeqPosStart)

#Finding last nucleotide position of the reference sequence
refSeqPosEnd <- foreach(i=1:nrow(dfRefSeq)) %do% (nchar(dfRefSeq$nucleotides[i]) + refSeqPosStart[i])
refSeqPosEnd <- as.numeric(refSeqPosEnd)

#Then we can substr the alignment by these positions to effectively trim the alignment
alignment2Trim <- foreach(i=1:nrow(dfRefSeq)) %do% substr(alignment2[[i]], refSeqPosStart[i]+1, refSeqPosEnd[i])

#Again convert to dnaStringSet format
dnaStringSet3 <- sapply( alignment2Trim, function(x) DNAStringSet( x ) )

#Remove our reference sequence from this as we dont want this to be included in further analysis
refSeqRemove <- foreach(i=1:nrow(dfRefSeq)) %do% which(dnaStringSet3[[i]]@ranges@NAMES == "reference")
dnaStringSet3 <- foreach(i=1:nrow(dfRefSeq)) %do% subset(dnaStringSet3[[i]][-refSeqRemove[[i]]])
alignment2Trim <- foreach(i=1:nrow(dfRefSeq)) %do% alignment2Trim[[i]][-refSeqRemove[[i]]]

###############
#Pairwise Distance Determination with TN93

#Conversion to DNAbin format before using genetic distance matrix
dnaBIN <- foreach(i=1:length(taxaListComplete)) %do% as.DNAbin(dnaStringSet3[[i]])

#Now for the computation of genetic distance, several models can be used - "raw", 
#"N", "TS", "TV", "JC69", "K80" (the default), 
#"F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet"
#Details on each model can be found here: 
#http://www.inside-r.org/packages/cran/ape/docs/dist.dna
#Using the TN93 model for our data
matrixGeneticDistance <- foreach(i=1:length(taxaListComplete)) %do% 
  dist.dna(dnaBIN[[i]], model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#convert to dataframe
geneticDistanceMatrixList <- foreach(i=1:length(taxaListComplete)) %do%  as.data.frame(matrixGeneticDistance[i])

#Putting it into a stack (each column concatenated into one long column of indexes and values) 
#so it can be easily subsetted
geneticDistanceStackList <- foreach(i=1:length(taxaListComplete)) %do% stack(geneticDistanceMatrixList[[i]])

################
#Finding appropriate pairings according to genetic distance criteria (latitude done further down in script)
#and dividing based on lineage

#These values can easily be edited to add more or less stringency to the matches
#Will produce lists with indexes (bins, row and column names) of each match according to genetic distance criteria of 0.15
pairingResultCandidates <- foreach(i=1:length(taxaListComplete)) %do% which(geneticDistanceMatrixList[[i]]<=0.15)

#Next we can delete the orders that dont have any pairings
pairingResultCheck <- foreach(i=1:length(taxaListComplete)) %do% which(length(pairingResultCandidates[[i]]) == 0)
pairingResultCheck <- which(pairingResultCheck>0)
#Subset taxalistcomplete, geneticDistance, pairingResultCandidate, dnaStringSet lists according to pairingResultCheck
#so that we arent keeping classes without pairings in these lists
if(length(pairingResultCheck>0)){
  taxaListComplete <- taxaListComplete[-pairingResultCheck]
  geneticDistanceStackList <- geneticDistanceStackList[-pairingResultCheck]
  pairingResultCandidates <- pairingResultCandidates[-pairingResultCheck]
  dnaStringSet3 <- dnaStringSet3[-pairingResultCheck]
}

#Use these pairing results and reference against the genetic distance stack list to subset it
pairingResultCandidates2 <- foreach(i=1:length(taxaListComplete)) %do% 
    geneticDistanceStackList[[i]][c(pairingResultCandidates[[i]]),]
#this command will grab the precise row and column numbers of each pairing so they can be used to determine a unique identifier
#for each pairing
pairingResultCandidates3 <- foreach(i=1:length(taxaListComplete)) %do%
    data.table(which(geneticDistanceMatrixList[[i]]<=0.15, arr.ind = TRUE))

#Merge pairingResultCandidate list 2 and 3 together and combine into the pairing results dataframe
dfPairingResultsL1L2 <- do.call(rbind, Map(data.frame, pairingResultCandidates2, pairingResultCandidates3))
#Creating a key that can be used to uniquely identify each pairing based on its ingroup distance and its precise position in its respective pairwise distance matrix
#from the row and column number
#This is to avoid an issue where bins meeting the the 0.15 criteria were being incorrectly paired together because of duplicated ingroupdistance values
#Also settings number of significant digits to longer so values are more likely to be unique
options(digits = 10)
dfPairingResultsL1L2$pairingKey <- dfPairingResultsL1L2$values * (dfPairingResultsL1L2$row + dfPairingResultsL1L2$col)

#combine AllSeq with pairingResultCandidates to get all relevant taxonomic and latitudinal data
dfPairingResultsL1L2 <- suppressWarnings(merge(dfPairingResultsL1L2, dfAllSeq, by.x = "ind", by.y = "bin_uri"))

#order by pairingKey so all pairings are ordered correctly
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$pairingKey),]

#Changing significant digits back
options(digits =5)
#get rid of zero ingroupdistance entries
zeroDistance <- which(dfPairingResultsL1L2$values == 0)
if(length(zeroDistance)>0){
  dfPairingResultsL1L2 <- dfPairingResultsL1L2[-zeroDistance,]
}

#Then we can multiply these ingroupdistance values by 1.3 to determine the minimum outgroup 
#distance from the pairings, 
#we put this in another column in the pairing results dataframe
#This minimum outgroup distance could also be a user adjustable parameter and can be easily
#modified
dfPairingResultsL1L2$inGroupDistx1.3 <- dfPairingResultsL1L2$values * 1.3

#Also grouping dataframe every 2 rows to reflect each unique pairing/matching, pairing 
#column will give a number value for each pairing ordered 
dfPairingResultsL1L2$inGroupPairing <- rep(1:(nrow(dfPairingResultsL1L2)/2), each = 2)

#Reorganizing and renaming some columns in pairing result candidate dataframe to make more easily 
#readable 
colnames(dfPairingResultsL1L2)[1] <- "bin_uri"
dfPairingResultsL1L2 <- (dfPairingResultsL1L2[,c("inGroupPairing","record_id","bin_uri","values",
                                                 "inGroupDistx1.3","medianLatAbs","medianLatMap","latMin","latMax",
                                                 "binSize","phylum_taxID","phylum_name","class_taxID",
                                                 "class_name","order_taxID","order_name","family_taxID",
                                                 "family_name","subfamily_taxID","subfamily_name",
                                                 "genus_taxID","genus_name","species_taxID","species_name",
                                                 "nucleotides","ind","medianLon")])
colnames(dfPairingResultsL1L2)[4] <- "inGroupDist"
colnames(dfPairingResultsL1L2)[26] <- "indexNo"
dfPairingResultsL1L2$index <- seq.int(nrow(dfPairingResultsL1L2))

#Merge aligned and trimmed nucleotide sequences with pairing candidates, replacing raw sequences from intial tsv
alignment2TrimUnlist <- unlist(alignment2Trim)
#New dataframe with trimmed and aligned sequences
dfTrimmedSeq <- data.frame(alignment2TrimUnlist)
dfTrimmedSeq$bin_uri <- names(alignment2TrimUnlist)
colnames(dfTrimmedSeq)[1] <- "trimmedNucleotides"

#merge sequences to dfPairingResultsL1L2
dfPairingResultsL1L2 <- merge(dfPairingResultsL1L2, dfTrimmedSeq, 
                              by.x = "bin_uri", by.y = "bin_uri")
#Again reorganization of pairingresults to have trimmedNucleotides added instead of nucleotides
dfPairingResultsL1L2 <- (dfPairingResultsL1L2[,c("inGroupPairing","record_id","bin_uri","inGroupDist",
                                                 "inGroupDistx1.3","medianLatAbs","medianLatMap","latMin","latMax",
                                                 "binSize","phylum_taxID","phylum_name","class_taxID",
                                                 "class_name","order_taxID","order_name","family_taxID",
                                                 "family_name","subfamily_taxID","subfamily_name",
                                                 "genus_taxID","genus_name","species_taxID","species_name",
                                                 "trimmedNucleotides","indexNo","medianLon","index")])
#order by ingrouppairing between pairings
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupPairing),]

#Dividing Pairing Candidates into two lineages:

#First lineage, this is done by picking even numbered rows of pairing results dataframe: 2,4,6...
dfPairingResultsL1 <- dfPairingResultsL1L2[dfPairingResultsL1L2$index%%2==0,]
#Second lineage, this is done with odd numbered rows
dfPairingResultsL2 <- dfPairingResultsL1L2[dfPairingResultsL1L2$index%%2>0,]

##############
#Eliminating Pairings based on Latitude Difference and Overlapping Latitudinal Range
#This section may take a few mins to process

#Determing latitude difference between pairings
latitudeDiffCheck <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  abs(dfPairingResultsL1$medianLatAbs[i] - dfPairingResultsL2$medianLatAbs[i])

#Setting 20 degrees criteria
#Ensuring all pairings meet the latitude difference of 20 degrees (can be modified by user) 
latitudeDiffCheck <- which(latitudeDiffCheck<20)
if(length(latitudeDiffCheck)>0){
  #subsetting Lineage 1 by this latitude check
  dfPairingResultsL1 <- dfPairingResultsL1[-latitudeDiffCheck,]
  #Now subsetting Lineage 2
  dfPairingResultsL2 <- dfPairingResultsL2[-latitudeDiffCheck,]
  #Now subsetting for dataframe with both lineages
  dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, dfPairingResultsL1L2$inGroupPairing %in% 
                                   dfPairingResultsL1$inGroupPairing)
}

#Next we can work on establishing latitudinal ranges for each pairing
#If the two lineages of a pairing have overlapping latitude regions of greater than 25% then 
#we would not consider that pairing as a viable pairing
#This is because we want each lineage of a pairing to meet an appropriate difference in latitude

#we can define the overlap range threshold as 25% of the latitude range of L1
#If an overlap is greater than this value we would discard with this pairing
#Of course this value could be easily modified to add more or less stringency to the script
rangeThresholdL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  (abs(dfPairingResultsL1$latMax[i] - dfPairingResultsL1$latMin[i]) * 0.25) 

#Same process for L2
rangeThresholdL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  (abs(dfPairingResultsL2$latMax[i] - dfPairingResultsL2$latMin[i]) * 0.25) 

#Define our latitude ranges for each lineage of a pairing
rangeL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  range(dfPairingResultsL1$latMax[i], dfPairingResultsL1$latMin[i])
rangeL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  range(dfPairingResultsL2$latMax[i], dfPairingResultsL2$latMin[i])

#Then we can determine the overlap region between them using the Overlap function from 
#the Desctools package
#Overlap will return an absolute value so we dont have to worry about negatives for 
#overlap values
overlapValueL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% Overlap(rangeL1[[i]], 
                                                                     rangeL2[[i]])
overlapValueL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% Overlap(rangeL2[[i]], 
                                                                     rangeL1[[i]])

#Then if there is a range overlap between two lineages in a pairing, we can determine
#if this overlap is actually larger than the 25% value of rangeThreshold for each 
#individual pairing
rangeOverlapCheckL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(rangeThresholdL1[[i]]<overlapValueL1[[i]])

#same process for L2
rangeOverlapCheckL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(rangeThresholdL2[[i]]<overlapValueL2[[i]])

#Then overlap values exceeding that 25% value we set should be returned as an 
#integer of 1, if an overlap does not exceed this value, then it will return an integer of 0
#If there is a value that meets these criteria, its associated pairing will be removed 
#from the pairing results dataframe
#This will name each element of rangeOverlapCheck with the inGroupPairing number to 
#identify the pairing we need to eliminate
names(rangeOverlapCheckL1) <- paste0(dfPairingResultsL1$inGroupPairing)
names(rangeOverlapCheckL2) <- paste0(dfPairingResultsL2$inGroupPairing)
#Identify which pairing in rangeOverlapCheck is greater than 0, do this with respect to both lineages
overlapIndL1 <- which(rangeOverlapCheckL1>0)
overlapIndL2 <- which(rangeOverlapCheckL2>0)
overlapIndTotal <- append(overlapIndL1,overlapIndL2)
overlapIndTotal <- unique(overlapIndTotal)
#if overlapIndTotal is not empty:
if(length(overlapIndTotal)>0){
  #Eliminate based on that pairing number(s) for each pairing result dataframe
  dfPairingResultsL1 <- dfPairingResultsL1[-overlapIndTotal,]
  dfPairingResultsL2 <- dfPairingResultsL2[-overlapIndTotal,]
  dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, dfPairingResultsL1L2$inGroupPairing %in% 
                                   dfPairingResultsL1$inGroupPairing)
}

################
#Ensuring Every Pairing is a unique set of bins

#Once again order by ingroup distance
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupDist),]
#Take the first instance of each bin, the first instance of it being the lowest ingroup distance it has in a pairing
dfPairingResultsL1L2 <- by(dfPairingResultsL1L2, dfPairingResultsL1L2["bin_uri"], head, n=1)
dfPairingResultsL1L2 <- Reduce(rbind, dfPairingResultsL1L2)
dfPairingResultsL1 <- subset(dfPairingResultsL1L2,duplicated(dfPairingResultsL1L2$inGroupPairing))
#Subsetting against L1L2 to get finalized pairings in L1L2
dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, dfPairingResultsL1L2$inGroupPairing %in% dfPairingResultsL1$inGroupPairing)
#now rebuilding L1 and L2 based on a subsetted L1L2
dfPairingResultsL1 <- dfPairingResultsL1L2[dfPairingResultsL1L2$index%%2==0,]
dfPairingResultsL2 <- dfPairingResultsL1L2[dfPairingResultsL1L2$index%%2>0,]

#Order and renumber pairings starting at 1 
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupDist),]
dfPairingResultsL1L2$inGroupPairing <- rep(1:(nrow(dfPairingResultsL1L2)/2), each = 2)
dfPairingResultsL1 <- dfPairingResultsL1[order(dfPairingResultsL1$inGroupDist),]
dfPairingResultsL1$inGroupPairing <- 1:nrow(dfPairingResultsL1)
dfPairingResultsL2 <- dfPairingResultsL2[order(dfPairingResultsL2$inGroupDist),]
dfPairingResultsL2$inGroupPairing <- 1:nrow(dfPairingResultsL2)

################
#Outgroup determination for each Pairing

#Now we can search for the best possible outgroupings for each pairing
#This will involve searching for outgroups relative to each pairing and finding one that
#is far enough away from both lineages (1.3x criteria)

#Fist we can take all genetic distance values from each distance matrix and concatenate them into one
#stack of values to be easily subsetted from
dfGeneticDistanceStack <- do.call(rbind, geneticDistanceStackList)

#Then we can search for outgroupings relative to lineage 1
#We can use the bin_uri's to subset the dfGeneticDistanceStackList according to 
#bin_uris represented in lineage 1, this will be called dfOutGroupL1
#This will essentially limit our outgroup distances to those associated with lineage1
dfOutGroupL1 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% 
                         dfPairingResultsL1$bin_uri)

#Order outgroups by Lineage1
dfOutGroupL1 <- dfOutGroupL1[order(match(dfOutGroupL1[,2],
                                         dfPairingResultsL1[,3])),]

#Find which bins in the outgroup genetic distances match the bins in the pairing results
outGroupCandidatesL1a <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfOutGroupL1$ind == dfPairingResultsL1$bin_uri[i])

#Find which outgroups meet the 1.3x criteria
outGroupCandidatesL1b <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfOutGroupL1$values >= dfPairingResultsL1$inGroupDistx1.3[i])

#Intersection of the two using mapply to find the correct outgroupings for each pairing
outGroupCandidatesL1c <- mapply(intersect,outGroupCandidatesL1a,outGroupCandidatesL1b)
#Unlist to make into one vector 
outGroupCandidatesL1c <- unlist(outGroupCandidatesL1c)
#Adding an rownum column to dfBestOutGroup, this represents the second index of 
#each outgroup candidate
dfOutGroupL1$rownum<-seq.int(nrow(dfOutGroupL1))
#Then we can subset to rownum column based on the outgroup candidates
dfOutGroupL1 <- dfOutGroupL1[dfOutGroupL1$rownum %in% 
                               outGroupCandidatesL1c, ]
#Same process for lineage 2
dfOutGroupL2 <- subset(dfGeneticDistanceStack, dfGeneticDistanceStack$ind %in% 
                         dfPairingResultsL2$bin_uri)
dfOutGroupL2 <- dfOutGroupL2[order(match(dfOutGroupL2[,2],
                                         dfPairingResultsL2[,3])),]
outGroupCandidatesL2a <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfOutGroupL2$ind == dfPairingResultsL2$bin_uri[i])
outGroupCandidatesL2b <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfOutGroupL2$values >= dfPairingResultsL2$inGroupDistx1.3[i])
outGroupCandidatesL2c <- mapply(intersect,outGroupCandidatesL2a,outGroupCandidatesL2b)
outGroupCandidatesL2c <- unlist(outGroupCandidatesL2c)
dfOutGroupL2$rownum <- seq.int(nrow(dfOutGroupL2))
dfOutGroupL2 <- dfOutGroupL2[dfOutGroupL2$rownum %in% 
                               outGroupCandidatesL2c, ]

#this next command will determine which outgroup candidates are shared between both lineages
dfOutGroupL1 <- dfOutGroupL1[dfOutGroupL1$rownum %in% dfOutGroupL2$rownum, ]
dfOutGroupL2 <- dfOutGroupL2[dfOutGroupL2$rownum %in% dfOutGroupL1$rownum, ]

#Then we can determine averages between outgroup distances and find the outgroup with the min average distance 
#to ensure the outgroup is still relatively close in distance to each lineage but at least at the 1.3x divergence value determined
#for each pairing
dfOutGroupL1$distAverage <- (dfOutGroupL1$values + dfOutGroupL2$values) / 2
dfOutGroupL2$distAverage <- (dfOutGroupL1$values + dfOutGroupL2$values) / 2
dfOutGroupL1 <- dfOutGroupL1[!duplicated(dfOutGroupL1$distAverage), ]
dfOutGroupL2 <- dfOutGroupL2[!duplicated(dfOutGroupL2$distAverage), ]
outGroupList <- lapply(unique(dfOutGroupL1$ind), function(x) dfOutGroupL1[dfOutGroupL1$ind == x,])
minOutGroupAverage <- sapply( outGroupList , function(x) min( x$distAverage ) )

#Subsetting each outgroup dataframe by our min outgroup averages
dfOutGroupL1 <- subset(dfOutGroupL1, distAverage %in% minOutGroupAverage)
dfOutGroupL2 <- subset(dfOutGroupL2, distAverage %in% minOutGroupAverage)

#Next find the bin_uri's for each outgroup candidate by subsetting the dfGeneticDistanceStack dataframe with values from the L1 and L2 dataframes
#L1
dfOutGroupL1 <- merge(dfOutGroupL1, dfGeneticDistanceStack, by.x = "values", by.y = "values", all.x = TRUE)
dfOutGroupL1 <- dfOutGroupL1[order(match(dfOutGroupL1[,2], dfPairingResultsL1[,3])),]
outGroupL1Check <- which(dfOutGroupL1$ind.x != dfOutGroupL1$ind.y)
dfOutGroupL1 <- dfOutGroupL1[outGroupL1Check,]

#L2 (same process)
dfOutGroupL2 <- merge(dfOutGroupL2, dfGeneticDistanceStack, by.x = "values", by.y = "values", all.x = TRUE)
dfOutGroupL2 <- dfOutGroupL2[order(match(dfOutGroupL2[,2], dfPairingResultsL2[,3])),]
outGroupL2Check <- which(dfOutGroupL2$ind.x != dfOutGroupL2$ind.y)
dfOutGroupL2 <- dfOutGroupL2[outGroupL2Check,]

#merging of outgroup dataframes so final filtering of outgroups can be performed
dfOutGroupMerge <- merge(dfOutGroupL1, dfOutGroupL2, by.x = "ind.y", by.y = "ind.y", all.x = TRUE)
outGroupMergeCheck <- which(dfOutGroupMerge$rownum.x == dfOutGroupMerge$rownum.y)
dfOutGroupMerge <- dfOutGroupMerge[outGroupMergeCheck,]
dfOutGroupMerge <- dfOutGroupMerge[!duplicated(dfOutGroupMerge$distAverage.x), ]

#Reorganization and merging of outgroup dataframes with dfAllSeq to get all taxonomic information
#L1
dfOutGroupL1 <- (dfOutGroupMerge[,c("ind.y","values.x","ind.x.x")])
dfOutGroupL1 <- merge(dfOutGroupL1, dfAllSeq, by.x = "ind.y", by.y = "bin_uri")
#Then merge with dfTrimmedSeq to replace raw sequences with trimmed and aligned ones
dfOutGroupL1 <- merge(dfOutGroupL1, dfTrimmedSeq, by.x = "ind.y", by.y = "bin_uri")
dfOutGroupL1 <- subset(dfOutGroupL1, select = -c(nucleotides) )
#renaming a few columns
colnames(dfOutGroupL1)[1] <- "outGroupBin"
colnames(dfOutGroupL1)[2] <- "outGroupDistance"
colnames(dfOutGroupL1)[3] <- "associatedInGroupBin"

#Same process for L2
dfOutGroupL2 <- (dfOutGroupMerge[,c("ind.y","values.y","ind.x.y")])
dfOutGroupL2 <- merge(dfOutGroupL2, dfAllSeq, by.x = "ind.y", by.y = "bin_uri")
dfOutGroupL2 <- merge(dfOutGroupL2, dfTrimmedSeq, by.x = "ind.y", by.y = "bin_uri")
dfOutGroupL2 <- subset(dfOutGroupL2, select = -c(nucleotides) )
#renaming a few columns
colnames(dfOutGroupL2)[1] <- "outGroupBin"
colnames(dfOutGroupL2)[2] <- "outGroupDistance"
colnames(dfOutGroupL2)[3] <- "associatedInGroupBin"

#remove outGroupMerge and dfTrimmedSeq since these are redundant and no longer needed
rm(dfOutGroupMerge)
rm(dfTrimmedSeq)

#Some classes may not have any viable outgroupings so we can filter those out
noOutGroupCheck <- setdiff(dfPairingResultsL1$bin_uri, dfOutGroupL1$associatedInGroupBin)
noOutGroupCheck <- foreach(i=1:length(noOutGroupCheck)) %do% which(dfPairingResultsL1$bin_uri == noOutGroupCheck[[i]])
noOutGroupCheck <- unlist(noOutGroupCheck)

#If there is at least one pairing without an outgroup then subset pairing dataframes by those pairings
if(length(noOutGroupCheck)>0){
  #subsetting Lineage 1 by this outgroup check
  dfPairingResultsL1 <- dfPairingResultsL1[-noOutGroupCheck,]
  #Now subsetting Lineage 2 by this outgroup check
  dfPairingResultsL2 <- dfPairingResultsL2[-noOutGroupCheck,]
  #Now subsetting for dataframe with both lineages
  dfPairingResultsL1L2 <- subset(dfPairingResultsL1L2, dfPairingResultsL1L2$inGroupPairing %in% 
                                   dfPairingResultsL1$inGroupPairing)
}

#Renumber pairings again
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupDist),]
dfPairingResultsL1L2$inGroupPairing <- rep(1:(nrow(dfPairingResultsL1L2)/2), each = 2)
dfPairingResultsL1 <- dfPairingResultsL1[order(dfPairingResultsL1$inGroupDist),]
dfPairingResultsL1$inGroupPairing <- 1:nrow(dfPairingResultsL1)
dfPairingResultsL2 <- dfPairingResultsL2[order(dfPairingResultsL2$inGroupDist),]
dfPairingResultsL2$inGroupPairing <- 1:nrow(dfPairingResultsL2)

###############
#Finalizing the Pairing Results Dataframes 

#First order outgroup lineage dataframes again by bin
dfOutGroupL1 <- dfOutGroupL1[order(match(dfOutGroupL1[,3],
                                         dfPairingResultsL1[,3])),]
dfOutGroupL2 <- dfOutGroupL2[order(match(dfOutGroupL2[,3],
                                         dfPairingResultsL2[,3])),]

#Adding an additional column for the delta lat difference of each pairing
latDelta <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  abs(dfPairingResultsL1$medianLatAbs[i] - dfPairingResultsL2$medianLatAbs[i])

#Finally merge outgroups to both L1, L2 and L1L2 dataframes
#Using the data tables package for this for reliable merging
#Making these dataframes data tables
dfPairingResultsL1 <- data.table(dfPairingResultsL1)
dfOutGroupL1 <- data.table(dfOutGroupL1)
dfOutGroupL1$inGroupPairing <- seq.int(nrow(dfOutGroupL1))
#Setting keys for merge
setkey(dfOutGroupL1, inGroupPairing)
setkey(dfPairingResultsL1, inGroupPairing)
dfPairingResultsL1 <- merge(dfPairingResultsL1, dfOutGroupL1, by.x = "inGroupPairing", 
                            by.y = "inGroupPairing")
dfPairingResultsL1 <- data.frame(dfPairingResultsL1)
dfPairingResultsL1$latDelta <- latDelta

#L2, same process
dfPairingResultsL2 <- data.table(dfPairingResultsL2)
dfOutGroupL2 <- data.table(dfOutGroupL2)
dfOutGroupL2$inGroupPairing <- seq.int(nrow(dfOutGroupL2))
setkey(dfOutGroupL2, inGroupPairing)
setkey(dfPairingResultsL2, inGroupPairing)
dfPairingResultsL2 <- merge(dfPairingResultsL2, dfOutGroupL2, by.x = "inGroupPairing", 
                            by.y = "inGroupPairing")
dfPairingResultsL2 <- data.frame(dfPairingResultsL2)
dfPairingResultsL2$latDelta <- latDelta

#Find pairings that where each lineage is equal in distance to the outgroup and remove these
equalOutGroupDistCheck <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  (dfPairingResultsL1$outGroupDistance[i] - dfPairingResultsL2$outGroupDistance[i])
equalOutGroupDistCheck <- which(equalOutGroupDistCheck == 0)
#If there is at least one pairing with equal outgroup distance from both lineages
if(length(equalOutGroupDistCheck)>0){
  #subsetting Lineage 1 by this outgroup check
  dfPairingResultsL1 <- dfPairingResultsL1[-equalOutGroupDistCheck,]
  #Now subsetting Lineage 2 by this outgroup check
  dfPairingResultsL2 <- dfPairingResultsL2[-equalOutGroupDistCheck,]
}
#Renumber pairing again
dfPairingResultsL1$inGroupPairing <- seq.int(nrow(dfPairingResultsL1))
dfPairingResultsL2$inGroupPairing <- seq.int(nrow(dfPairingResultsL2))

#Now some dataframe reorganization and ordering
dfPairingResultsL1 <- (dfPairingResultsL1[,c("inGroupPairing","associatedInGroupBin","inGroupDist","inGroupDistx1.3",
                                             "medianLatAbs.x","latDelta","latMin.x","latMax.x","record_id.x",
                                             "binSize.x","phylum_taxID.x","phylum_name.x","class_taxID.x",
                                             "class_name.x","order_taxID.x","order_name.x","family_taxID.x",
                                             "family_name.x","subfamily_taxID.x","subfamily_name.x",
                                             "genus_taxID.x","genus_name.x","species_taxID.x","species_name.x",
                                             "trimmedNucleotides.x","medianLatMap.x","medianLon.x",
                                             "outGroupBin","outGroupDistance","record_id.y","binSize.y",
                                             "phylum_taxID.y","phylum_name.y","class_taxID.y",
                                             "class_name.y","order_taxID.y","order_name.y","family_taxID.y",
                                             "family_name.y","subfamily_taxID.y","subfamily_name.y",
                                             "genus_taxID.y","genus_name.y","species_taxID.y","species_name.y",
                                             "trimmedNucleotides.y")])
colnames(dfPairingResultsL1)[2] <- "inGroupBin"
dfPairingResultsL1 <- dfPairingResultsL1[order(dfPairingResultsL1$inGroupPairing),]

#L2
dfPairingResultsL2 <- (dfPairingResultsL2[,c("inGroupPairing","associatedInGroupBin","inGroupDist","inGroupDistx1.3",
                                             "medianLatAbs.x","latDelta","latMin.x","latMax.x","record_id.x",
                                             "binSize.x","phylum_taxID.x","phylum_name.x","class_taxID.x",
                                             "class_name.x","order_taxID.x","order_name.x","family_taxID.x",
                                             "family_name.x","subfamily_taxID.x","subfamily_name.x",
                                             "genus_taxID.x","genus_name.x","species_taxID.x","species_name.x",
                                             "trimmedNucleotides.x","medianLatMap.x","medianLon.x",
                                             "outGroupBin","outGroupDistance","record_id.y","binSize.y",
                                             "phylum_taxID.y","phylum_name.y","class_taxID.y",
                                             "class_name.y","order_taxID.y","order_name.y","family_taxID.y",
                                             "family_name.y","subfamily_taxID.y","subfamily_name.y",
                                             "genus_taxID.y","genus_name.y","species_taxID.y","species_name.y",
                                             "trimmedNucleotides.y")])
colnames(dfPairingResultsL2)[2] <- "inGroupBin"
dfPairingResultsL2 <- dfPairingResultsL2[order(dfPairingResultsL2$inGroupPairing),]

#Rebuilding of L1L2 dataframe with outgroups
dfPairingResultsL1L2 <- rbind(dfPairingResultsL1, dfPairingResultsL2)
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupPairing),]

###############
#Identifying PseudoReplicates 

#To check for the phylogenetics problem of pseudoreplication we can generate another 
#smaller distance matrix with our selected pairings only for each class
#If a bin from one pairing is actually closer to a bin from another pairing 
#(as opposed to its paired lineage) 
#than we would have to average the results of those two pairings in the statistics section

#To do this we can subset the original pairwise distance matrices with the pairing 
#bins of each lineage only for each order

#First we need the bins associated with each genetic/pairwise distance matrix
geneticDistanceMatrixNames <- foreach(i=1:length(geneticDistanceMatrixList)) %do%  
  as.character(names(geneticDistanceMatrixList[[i]]))
#Then bins associated with each pairing lineage
pairingNamesL1 <- as.character(dfPairingResultsL1$inGroupBin)
pairingNamesL2 <- as.character(dfPairingResultsL2$inGroupBin)

#Then subset this list of bins based on bins present in each of the lineages of each pairing
geneticDistanceMatrixNamesL1 <- foreach(i=1:length(geneticDistanceMatrixNames)) %do% 
  subset(geneticDistanceMatrixNames[[i]], geneticDistanceMatrixNames[[i]] %in% pairingNamesL1)
geneticDistanceMatrixNamesL2 <- foreach(i=1:length(geneticDistanceMatrixNames)) %do% 
  subset(geneticDistanceMatrixNames[[i]], geneticDistanceMatrixNames[[i]] %in% pairingNamesL2)

#Then subset the list of distance matrices using the lists above
#Can call this new list of subsetted matrices pseudoRepMatrixList
pseudoRepMatrixList <- foreach(i=1:length(geneticDistanceMatrixList)) %do%  
  geneticDistanceMatrixList[[i]][geneticDistanceMatrixNamesL1[[i]], geneticDistanceMatrixNamesL2[[i]]]

#Then we have to remove matrices in this list which either only have one pairing (2 bins) or 
#no pairings present since these cannot have pseudoreplicates present
pseudoRepMatrixLengthCheck <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  which(length(pseudoRepMatrixList[[i]])<2)
pseudoRepMatrixLengthCheck <- sapply( pseudoRepMatrixLengthCheck , function(x) length( x ) )
pseudoRepMatrixLengthCheck <- which(pseudoRepMatrixLengthCheck>0)
if(length(pseudoRepMatrixLengthCheck)>0){
  pseudoRepMatrixList <- pseudoRepMatrixList[-pseudoRepMatrixLengthCheck]
}

#For each column and row of each distance matrix (each matrix belonging to a specific class), if a distance value is lower than the
#pairwise distances of each pairing than we know there is another bin which is actually closer
#We check for this by using our pairwise distances from the dfPairingResultsLineage
#dataframes

#First we check the minimum values for each row of each matrix
pseudoRepMatrixMinRow <- foreach(i=1:length(pseudoRepMatrixList)) %do% apply(pseudoRepMatrixList[[i]],1,min)

#Then we check the minimum values for each column of each matrix
pseudoRepMatrixMinColumn <- foreach(i=1:length(pseudoRepMatrixList)) %do% apply(pseudoRepMatrixList[[i]],2,min)

#if bin distance values intersect between these two minimums (row and column) then 
#we know there are no pseudoreplicates present for these bins
#if bin distance values do not intersect, then we know there is a pseudoreplicate bin that is actually closer in distance than one of the pairing bins

#these next few commands identify bins which do not intersect thus identifying pseudoreplicates
pseudoRepMatrixCandidates <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  which(pseudoRepMatrixMinRow[[i]] != pseudoRepMatrixMinColumn[[i]])
pseudoRepMatrixCandidates2 <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  pseudoRepMatrixMinColumn[[i]][pseudoRepMatrixCandidates[[i]]]
pseudoRepMatrixCandidates <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  names(pseudoRepMatrixCandidates[[i]])
pseudoRepMatrixCandidates <- unlist(pseudoRepMatrixCandidates)
pseudoRepMatrixCandidates2 <- foreach(i=1:length(pseudoRepMatrixList)) %do% 
  names(pseudoRepMatrixCandidates2[[i]])
pseudoRepMatrixCandidates2 <- unlist(pseudoRepMatrixCandidates2)

#Now we can get the inGroupPairing numbers associated with the pseudoreplicates
pseudoRepL1 <- subset(dfPairingResultsL1, dfPairingResultsL1$inGroupBin %in% pseudoRepMatrixCandidates)
pseudoRepL1 <- pseudoRepL1[order(match(pseudoRepL1[,2],pseudoRepMatrixCandidates)),]
pseudoRepL1 <- (pseudoRepL1[,c("inGroupPairing")])  
pseudoRepL2 <- subset(dfPairingResultsL2, dfPairingResultsL2$inGroupBin %in% pseudoRepMatrixCandidates2)
pseudoRepL2 <- pseudoRepL2[order(match(pseudoRepL2[,2],pseudoRepMatrixCandidates2)),]
pseudoRepL2 <- (pseudoRepL2[,c("inGroupPairing")])  

#If there is at least 1 pseudoreplicate found in any of the distance matrices
if(length(pseudoRepMatrixCandidates)>0){
  #Make the pseudoreplicates into a dataframe
  dfPseudoRep <- as.data.frame(pseudoRepL1)
  dfPseudoRep$inGroupPairing2 <- pseudoRepL2
  colnames(dfPseudoRep)[1] <- "inGroupPairing"
  #delete duplicate inGroupPairing across columns (this is a rare occurence but does happen sometimes)
  duplicatePseudoRep <- which(dfPseudoRep$inGroupPairing == dfPseudoRep$inGroupPairing2)
  dfPseudoRep <- dfPseudoRep[-duplicatePseudoRep,]
  
  #also to remove duplicated pairings across the inGroupPairing and inGroupPairing2 columns:
  #Ex: 1,2 and 2,1 in separate rows of dfPseudoRep since these are the same
  dfPseudoRep <- unique(t(apply(dfPseudoRep, 1, sort)))
  dfPseudoRep <- data.frame(dfPseudoRep)
  colnames(dfPseudoRep)[1] <- "inGroupPairing"
  colnames(dfPseudoRep)[2] <- "inGroupPairing2"
  dfPseudoRep <- dfPseudoRep[order(dfPseudoRep$inGroupPairing),]
  dfPseudoRep <- suppressWarnings(by(dfPseudoRep, dfPseudoRep["inGroupPairing2"], head, n=1))
  dfPseudoRep <- Reduce(rbind, dfPseudoRep)
}

#These identified pseudoreplicates will now have their signed 
#relative outgroup distances averaged in the statistics section

################
#Relative OutGroup Distance Determination for Pairings

#In this section, relative outgroup distance is determined for each pairing
#and also the sign of this distance is determined

#Positive relative distance indicates a pairing where the lower latitude lineage also has the 
#further outgroup distance thus supporting the latitudinal diversity gradient

#Negative would be where the lower latitude lineage has the smaller outgroup distance

#create variables that hold outgroup distances
testOutGroupDist<-as.numeric(dfPairingResultsL1L2$outGroupDistance) #stored as a vector
#set trueTestOutGroupDist to null to ensure empty variable
trueTestOutGroupDist<-NULL
#counter for trueTestOutGroupDist
x=1
#for loop that increments by 2 to divide adjacent outgroup distances found in
#outGroupDist
#This will give our relative genetic distances in the vector testOutGroupDist
for (i in seq(from=1, to=length(testOutGroupDist), by = 2)){
  if (testOutGroupDist[i] > testOutGroupDist[i+1]){
    trueTestOutGroupDist[x]<-testOutGroupDist[i]/testOutGroupDist[i+1]
  }
  else if (testOutGroupDist[i+1] > testOutGroupDist[i]){
    trueTestOutGroupDist[x]= testOutGroupDist[i+1]/testOutGroupDist[i]
  }
  else if (testOutGroupDist[i+1] == testOutGroupDist[i]){
    trueTestOutGroupDist[x]= testOutGroupDist[i+1]/testOutGroupDist[i]
  }
  x<-x+1
}

#Checking to see which latitudes are lower for lineage 1 compared to lineage2
lowerLatL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfPairingResultsL1$medianLatAbs.x[i]<dfPairingResultsL2$medianLatAbs.x[i]) 
#Checking to see which genetic distances are greater for lineage 1 compared to 
#lineage2
greaterDistanceL1 <- foreach(i=1:nrow(dfPairingResultsL1)) %do% 
  which(dfPairingResultsL1$outGroupDistance[i]>dfPairingResultsL2$outGroupDistance[i])
#Then an integer of 1 will be assigned to those in Lineage1 that meet those 
#criteria, meeting these criteria is defined as a success
successValL1 <- mapply(intersect,lowerLatL1,greaterDistanceL1)
#This will give the index of each pairing relative to Lineage 1 that met those criteria
successValL1 <- which(successValL1>0) 

#Do the same process for lineage 2 compared to lineage1
lowerLatL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do%
  which(dfPairingResultsL2$medianLatAbs.x[i]<dfPairingResultsL1$medianLatAbs.x[i]) 
greaterDistanceL2 <- foreach(i=1:nrow(dfPairingResultsL2)) %do% 
  which(dfPairingResultsL2$outGroupDistance[i]>dfPairingResultsL1$outGroupDistance[i])
successValL2 <- mapply(intersect,lowerLatL2,greaterDistanceL2)
successValL2 <- which(successValL2>0) 

#Append both success lists together representing our total number of successes
successOverall <- append(successValL1, successValL2)
#Sort the list
successOverall <- sort(successOverall)

#Also making a vector for pairings that were not successes
failureOverall <- dfPairingResultsL1$inGroupPairing[-successOverall]

#Can then subset our trueTestOutGroupDist vector containing relative distances 
#with our successess
relativeDistPos <- trueTestOutGroupDist[successOverall]
#Then subset based on which relative distances were not successes
relativeDistNeg <- trueTestOutGroupDist[failureOverall]
#These values will also be assigned a negative value
relativeDistNeg <- relativeDistNeg * -1
#Also identifying which pairing was positive and which was negative
names(relativeDistPos) <- paste0(successOverall)
names(relativeDistNeg) <- paste0(failureOverall)

#Creating a new named vector with pos and negative values appended together 
#but with pairing numbers as names 
relativeDistOverall <- append(relativeDistNeg, relativeDistPos)
#Make realtiveDistOverall into a dataframe
dfRelativeDist <- data.frame(relativeDistOverall)
#Adding rownames variable column of dfRelativeDistOverall
dfRelativeDist$variable <- rownames(dfRelativeDist)
#Change to numeric
dfRelativeDist$variable <- as.numeric(dfRelativeDist$variable)
#Order by this column
dfRelativeDist <- dfRelativeDist[order(dfRelativeDist$variable),]
#Adding relative distances to a new column called value
colnames(dfRelativeDist)[1] <- "value"
#Creating another value called sign to determine which is positive and which is negative
dfRelativeDist[["sign"]] = ifelse(dfRelativeDist[["value"]] >= 0, 
                                  "positive", "negative")

#Appending relative outgroup distance results to ResultsSummary dataframe and ResultsL1L2 dataframe
dfPairingResultsL1$relativeOutGroupDist <- dfRelativeDist$value
dfPairingResultsL2$relativeOutGroupDist <- dfRelativeDist$value
dfPairingResultsL1L2 <- rbind(dfPairingResultsL1, dfPairingResultsL2)
dfPairingResultsL1L2 <- dfPairingResultsL1L2[order(dfPairingResultsL1L2$inGroupPairing),]

#Also creating another dataframe called pairing results summary which will give a quick idea of pairing results
#and make the results easier to read, this will include relative outgroup distance
dfPairingResultsSummary <- (dfPairingResultsL1L2[,c("inGroupPairing","relativeOutGroupDist","inGroupBin","inGroupDist","inGroupDistx1.3",
                                                    "outGroupDistance","outGroupBin","medianLatAbs.x","latDelta","latMin.x",
                                                    "latMax.x","class_name.x","order_name.x","family_name.x","genus_name.x",
                                                    "species_name.x","trimmedNucleotides.x","class_name.x","order_name.y","family_name.y",
                                                    "genus_name.y","species_name.y","trimmedNucleotides.y")])
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
#Reordering
dfPairingResultsSummary <- dfPairingResultsSummary[order(dfPairingResultsSummary$inGroupPairing),]

#Merge Relative Dist dataframe to a pairing results dataframe to grab the order names so results can
#be broken down by class
dfRelativeDist <- merge(dfRelativeDist, dfPairingResultsL1, by.x = "variable", by.y = "inGroupPairing")
dfRelativeDist <- (dfRelativeDist[,c("variable","value","sign","class_name.x")])

##################
#PseudoReplicate Relative Outgroup Distance Averaging 

#Now to average pseudoreplicates so they can be added to dfRelativeDist
#Again if at least 1 pseudoreplicate is present in the dfPseudoRep dataframe
if(nrow(dfPseudoRep)>0){
  #First we have to add the signed relative outgroup distances to dfPseudoRep, merging dfRelativeDist
  #to do this
  dfPseudoRep <- merge(dfPseudoRep, dfRelativeDist, by.x = "inGroupPairing", by.y = "variable")
  dfPseudoRep <- merge(dfPseudoRep, dfRelativeDist, by.x = "inGroupPairing2", by.y = "variable")
  dfPseudoRep <- (dfPseudoRep[,c("inGroupPairing","value.x","inGroupPairing2","value.y","class_name.x.x")])
  colnames(dfPseudoRep)[2] <- "relativeDist1"
  colnames(dfPseudoRep)[4] <- "relativeDist2"
  #Also ordering dfPseudoRep
  dfPseudoRep <- dfPseudoRep[order(dfPseudoRep$inGroupPairing),]
  #rounding
  dfPseudoRep$relativeDist1 <- round(dfPseudoRep$relativeDist1, 4)
  dfPseudoRep$relativeDist2 <- round(dfPseudoRep$relativeDist2, 4)
  
  #Now subtracting 1 from positive values and adding one to negative values for averaging 
  #This is to address an issue where averages of pairings differing in sign were producing decimal results
  #RelativeDist1:
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if(dfPseudoRep$relativeDist1[i]>0){
      dfPseudoRep$relativeDist1[i] <- dfPseudoRep$relativeDist1[i] - 1
    }
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if(dfPseudoRep$relativeDist1[i]<0){
      dfPseudoRep$relativeDist1[i] <- dfPseudoRep$relativeDist1[i] + 1
    }
  #RelativeDist2
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if(dfPseudoRep$relativeDist2[i]>0){
      dfPseudoRep$relativeDist2[i] <- dfPseudoRep$relativeDist2[i] - 1
    }
  foreach(i=1:nrow(dfPseudoRep)) %do% 
    if(dfPseudoRep$relativeDist2[i]<0){
      dfPseudoRep$relativeDist2[i] <- dfPseudoRep$relativeDist2[i] + 1
    }
  
  #breaking pseudoreplicates down to a list
  pseudoRepList <- lapply(unique(dfPseudoRep$inGroupPairing), 
                          function(x) dfPseudoRep[dfPseudoRep$inGroupPairing == x,])
  
  #Number of pseudoreplicates per pairing
  pseudoRepLength <- sapply( pseudoRepList , function (x) length( x$relativeDist2 ) )
  #distances and names associated with pseudoreplicates only
  pseudoRepRelativeDist2 <- sapply( pseudoRepList , function (x) ( x$relativeDist2 ) )
  pseudoRepRelativeDist2Names <- sapply( pseudoRepList , function (x) ( x$inGroupPairing2 ) )
  pseudoRepRelativeDist2Names <- pseudoRepRelativeDist2Names
  #distances and names associated with ingrouppairings only
  pseudoRepRelativeDist1 <- unique(dfPseudoRep$relativeDist1)
  pseudoRepRelativeDist1Names <- sapply( pseudoRepList , function (x) ( x$inGroupPairing ) )
  pseudoRepRelativeDist1Names <- sapply( pseudoRepRelativeDist1Names , function (x) unique( x ) )
  #append these distances together for averaging using map
  #map will append to a list format
  pseudoRepAllRelativeDist = suppressWarnings(Map(c, pseudoRepRelativeDist2, pseudoRepRelativeDist1))
  pseudoRepAllRelativeDistNames = Map(c, pseudoRepRelativeDist2Names, pseudoRepRelativeDist1Names)
  
  #Now we can finally average the relative distances based on the values in the pseudoRepAllRelativeDist list
  pseudoRepAverage <- sapply( pseudoRepAllRelativeDist , function(x) mean( x ) )
  
  #Making another dataframe with the averages for the pseudoreplicates called dfPseudoRepAverage
  dfPseudoRepAverage <- data.frame(pseudoRepAverage)
  
  #Now subtracting 1 if negative or adding 1 if positive to averages
  #Once again this will ensure that all pseudoreplicate average values are above 1 or below -1
  foreach(i=1:nrow(dfPseudoRepAverage)) %do% 
    if(dfPseudoRepAverage$pseudoRepAverage[i]<0){
      dfPseudoRepAverage$pseudoRepAverage[i] <- dfPseudoRepAverage$pseudoRepAverage[i] - 1
    }
  foreach(i=1:nrow(dfPseudoRepAverage)) %do% 
    if(dfPseudoRepAverage$pseudoRepAverage[i]>0){
      dfPseudoRepAverage$pseudoRepAverage[i] <- dfPseudoRepAverage$pseudoRepAverage[i] + 1
    }
  
  #Adding another column for the pairings associated with each average
  dfPseudoRepAverage$variable <- pseudoRepAllRelativeDistNames
  trim <- function (x) sub("[c(]","", x)
  dfPseudoRepAverage$variable <- trim(dfPseudoRepAverage$variable)
  colnames(dfPseudoRepAverage)[1] <- "value"
  #Also adding the signs of each average for plotting later on
  dfPseudoRepAverage[["sign"]] = ifelse(dfPseudoRepAverage[["value"]] >= 0, 
                                        "positive", "negative")
  #Now lets subset dfRelativeDist to remove pairings that were averaged
  #First making a vector with all pairings averaged
  pseudoRepPairings <- append(dfPseudoRep$inGroupPairing, dfPseudoRep$inGroupPairing2)
  pseudoRepPairings <- pseudoRepPairings[!duplicated(pseudoRepPairings)]
  #Then subsetting dfRelativeDist based on this
  dfRelativeDist <- dfRelativeDist[setdiff(dfRelativeDist$variable, pseudoRepPairings),]  
  #Adding class names to dfPseudoRepAverage
  dfPseudoRepAverage$className <- sapply( pseudoRepList, function(x) unique( x$class_name.x.x ) )
  #Column name change for dfRelativeDist
  colnames(dfRelativeDist)[4] <- "className"
  #Now dfPseudoRepAverage will be appended to the relativeDist dataframe
  dfRelativeDist <- rbind(dfRelativeDist,dfPseudoRepAverage)
  #Can remove dfPseudoRep since it is redundant with dfPseudoRepAverage
  rm(dfPseudoRep)
}

##############
#Statistical Analysis of Pairings

#Binomial Test

#Performing a binomial test on our pairings to test to see if pairings with lower 
#latitude AND greater outgroup distance are more prevalent than the null expectation of 50% 
#Null expectation being that pairings with high lat/high outgroup dist are equally 
#prevalent compared with low lat/high outgroup dist pairings

#Binomial test on relative branch lengths for sister pairs
#Number of successes defined as number of positive values in the dfRelativeDist dataframe:
successOverall <- length(which(dfRelativeDist$sign == "positive"))
#total number of trials is number of rows of dfRelativeDist dataframe
#Null expectation set to 50%

#Binomial test for all orders
binomialTestOutGroup <- binom.test(successOverall, 
                                   nrow(dfRelativeDist), 
                                   p= 0.5)

#Binomial test for each order separately

#Break relative outgroup distances down by class
classBinomList <- lapply(unique(dfRelativeDist$className), 
                         function(x) dfRelativeDist[dfRelativeDist$className == x,])

#Total number of observations for each order
classNumObservations <- sapply( classBinomList , function (x) length( x$variable ) )

#Number of successes for each order
classNumSuccesses <- sapply( classBinomList , function (x) which( x$sign == "positive" ) )
classNumSuccesses <- sapply( classNumSuccesses , function (x) length( x ) )

#Binom test for each order
binomClass <- foreach(i=1:length(classNumObservations)) %do% 
  binom.test(classNumSuccesses[i],
             classNumObservations[i], 
             p=0.5)

#p-value extraction for binomial test from each order
pvalBinomialTotal <- binomialTestOutGroup$p.value
pvalBinomialClass <- sapply( binomClass , function (x) ( x$p.value ) )

#Creation of a dataframe for p-values
allClassVariable <- "AllClasses"
classNames1 <- sapply( classBinomList , function (x) unique( x$className ) )
classNames2 <- append(allClassVariable, classNames1)
pValBinomial <- append(pvalBinomialTotal, pvalBinomialClass)
dfPVal <- data.frame(classNames2)
dfPVal$pValueBinomial <- pValBinomial

#Wilcoxon Test

#Next we do a Wilcoxon test on all of the signed relative distances (pos or neg) from the
#value column of the relativeDist dataframe
#to compare the median for a significant difference from the null expectation of zero.
#This test will consider both the magnitude and direction but is also non-parametric

#wilcoxon test for all classes
wilcoxTestOutGroup <- wilcox.test(dfRelativeDist$value, mu=0)

#Wilcoxon test for each class separately:

#Relative outgroup values for each class to be fed into test
classValue <- sapply( classBinomList , function (x) ( x$value ) )

#Wilcoxon test for each class
wilcoxonClass <- foreach(i=1:length(classValue)) %do% wilcox.test(classValue[[i]], mu=0) 

#p-value extraction for wilcoxon test
pvalWilcoxon <- wilcoxTestOutGroup$p.value
pvalWilcoxonClass <- sapply( wilcoxonClass , function (x) ( x$p.value ) )
pValWilcoxonTotal <- append(pvalWilcoxon, pvalWilcoxonClass)

#Addition of wilcoxon p-values to dfPVal dataframe
dfPVal$pValueWilcoxon <- pValWilcoxonTotal 

###############
#Plotting of Relative Outgroup Distance Results

#Plotting of our relative distances based on pairing number
#Will plot red if the value is below 0 (meaning not a success) and 
#blue if above 0 (success!)

#Make the variable column a factor for dfRealtiveDistOverall so ggplot classes pairings correctly
dfRelativeDist$variable <- factor(dfRelativeDist$variable, 
                                  levels = dfRelativeDist$variable)

#breaking dfRelativeDist down to a list by order so plots can be generated according to class
relativeDistClass <- lapply(unique(dfRelativeDist$className), 
                            function(x) dfRelativeDist[dfRelativeDist$className == x,])

#Variables for the title of each plot, will include name of the class and associated p-values
classTitle <- foreach(i=1:length(classNames1)) %do% 
  paste("Relative Outgroup Distances of Each Pairing of", classNames1[i])
pValBinomialTitle <- foreach(i=1:length(pvalBinomialClass)) %do% 
  paste("Binomial Test p-value:", round(pvalBinomialClass[i], digits = 4))
pvalWilcoxonTitle <- foreach(i=1:length(pvalBinomialClass)) %do% 
  paste("Wilcoxon Test p-value:", round(pvalWilcoxonClass[i], digits = 4))

#Plot of signed relative distance per pairing using ggplot:

#Separate plots will be generated for each class 

foreach(i=1:length(relativeDistClass)) %do%  
  (print(ggplot(relativeDistClass[[i]], aes(x = variable, y = value, color = sign))
         + geom_point(stat="identity", size = 2.5) 
         + theme(text = element_text(size=13), axis.text.x = element_text(face="bold",angle=90, vjust=1))
         + ggtitle(paste0(classTitle[i],"\n", pValBinomialTitle[i], "\n", pvalWilcoxonTitle[i])) 
         + labs(x="Pairing Number", y="Signed Relative OutGroup Distance", color = "sign")))


#Click on zoom icon in the bottom right hand corner to see plot more clearly

#Plots can be scrolled through using the arrow icons in the bottom right hand corner

#Plots that are important to keep can then be exported to any directory using the Export function in the viewer
#(bottom right hand of RStudio)

#If "Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : invalid graphics state" message appears,
#use this command:
#dev.off()

################
#Map Plotting of latitude separated pairings using Plotly
#Info can be found here on this: https://plot.ly/r/scatter-plots-on-maps/

#Code adapted from Plotly website
#note there are more options for map customization that can be found on plotly
#this just represents a basic world map representation

#new variable for basic map layout characteristics:
#using eckert4 projection which will reduce area near poles and expand area near equator

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

#New dataframe column with data for hovering over points on the map
#can add more columns to hover if you want more detail on the map
dfPairingResultsL1L2$hover <- paste("PairNum:",dfPairingResultsL1L2$inGroupPairing,
                                    "Ingroup:",dfPairingResultsL1L2$inGroupBin, 
                                    dfPairingResultsL1L2$species_name.x, 
                                    round(dfPairingResultsL1L2$inGroupDist, 3), 
                                    "Outgroup:", dfPairingResultsL1L2$outGroupBin,
                                    dfPairingResultsL1L2$species_name.y,
                                    sep = "<br>")

#Title for the map:

#First taking all unique class names to name the map
uniqueClass <- unique(dfPairingResultsL1L2$class_name.x)
#Pasting them together
classConcatenate <- paste(uniqueClass , collapse = ',')

#Map title command
mapTitle <- paste("Map of Latitude Separated Sister Pairings for", classConcatenate)

#this command will show a map organized by the pairing number with each pairing number having its
#own distinctive color on a spectrum of purple (shortest distance) to red (longest distance)
#map can become crowded with large numbers of pairings so individual maps may need to be manually created
#for each class
plot_ly(dfPairingResultsL1L2, lat = medianLatMap.x, lon = medianLon.x, text = hover,
        color = as.ordered(inGroupPairing), colors = "Spectral", mode = "markers+lines", type = 'scattergeo') %>%
  layout(title = paste0(mapTitle) , geo = mapLayout)

#A map can be generated by class with each order having its own distinctive color as well
#but due to limitations with plotly, lines will not connect pairings:
plot_ly(dfPairingResultsL1L2, lat = medianLatMap.x, lon = medianLon.x, text = hover,
        color = class_name.x, mode = "markers", type = 'scattergeo') %>%
  layout(title = paste0(mapTitle) , geo = mapLayout)

#***you will need to click on the icon in the viewer (bottom right corner) that says "show in new window" 
#(little box with arrow beside the refresh icon), unfortunately does not show actual map directly in Rstudio
#the map will appear in a web browser window though you dont have to be online to do this***

###############
#***Only important if you want to post maps online***
#Posting of Map to Plotly Server for Online Viewing

#For uploading to plotly server for online viewing of map

#You will first have to create a plotly account to do this:
#https://plot.ly/

#Note there is a limit of one plot for the free version of the Plotly account
#and the plot is public meaning other people on plotly can view the plot though
#it is not easily found on the website without the direct link

#to obtain additional plots on the server you have to pay for a package

#Run these commands for uploading user details, enter username and API key in the empty quotations to run commands:
#(obtained from making an account and in settings of account details) 

#Sys.setenv("plotly_username"="") 
#Sys.setenv("plotly_api_key"="")

#run these commands to make plot as a variable:

#plot <- plot_ly(dfPairingResultsL1L2, lat = medianLatMap.x, lon = medianLon.x,
#text = hover,
#color = as.ordered(inGroupPairing), colors = "Spectral", mode = "markers+lines", type = 'scattergeo') %>%
#layout(title = 'Latitude Separated Sister Pairings', geo = mapLayout)

#run this command for posting of map to plotly server (can rename in quotations to a name you prefer):

#plotly_POST(plot, filename = "LatitudeSeparatedPairingsMap")

###############
#Output of pairing results to TSV or CSV
#***Only important if you want to export results outside of R***

#uncomment the appropriate commands for either TSV or CSV output

#First create this dataframe for dfPairingResultsL1L2, this will allow output 
#of this file:

#dfPairingResultsOut <- data.frame(lapply(dfPairingResultsL1L2, as.character),
#stringsAsFactors=FALSE)

#or for results summary:

#dfPairingResultsOut <- data.frame(lapply(dfPairingResultsSummary, as.character),
#stringsAsFactors=FALSE)

#Defining another variable to give a unique name to the CSV
#Name would be what you specify it to be, R will prompt you to insert a name for
#the file in the console:
#filename <- readline(prompt="")

#Then you can uncomment one of these write.table commands to output
#**will output file to current working directory of R**

#CSV

#write.table(dfPairingResultsOut, file=paste(filename, ".csv", sep=""), quote=FALSE,
#sep=',', col.names = NA)

#TSV

#write.table(dfPairingResultsOut, file=paste(filename, ".tsv", sep=""), quote=FALSE, 
#sep='\t', col.names = NA)
