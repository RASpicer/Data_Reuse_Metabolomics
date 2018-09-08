---
title: "Data Reuse in Metabolomics"
author: "Analysis by Rachel Spicer, github:RASpicer"
date: "03/06/2018"
output: 
  html_document:
    code_folding: hide
    number_sections: yes
    theme: cerulean
    keep_md: true
---


# Data Processing 
This RMarkdown contains the code used for analysis for the subsection <b>Data Reuse in Metabolomics</b> of Chapter 3 of the thesis <b>Fit for purpose? A metascientific analysis of metabolomics data in public repositories</b>. 

## Functions
This sections contains functions that are used for processing the reuse data.

```r
# Functions Extract which studies have been reused
ReusedStudies <- function(x) {
    Studies <- as.character(x$`Study(-ies)`)
    Studies <- unlist(strsplit(Studies, ","))
    # remove instances where a single study is not referenced
    Studies <- Studies[!grepl("-", Studies)]
    Studies <- Studies[!grepl("All", Studies)]
    Studies <- Studies[!grepl("GNPS", Studies)]
    return(Studies)
}

# Produce an ordered frequency table of the number of reuses per study
UniqueStudiesCount <- function(x) {
    Studies <- ReusedStudies(x)
    # Produce a frequency table of the number of reuses per study
    UniqueStudiesCount <- as.data.frame(table(Studies))
    # Order by Frequency of reuse
    UniqueStudiesCountOrdered <- UniqueStudiesCount[order(-UniqueStudiesCount$Freq), 
        ]
    return(UniqueStudiesCountOrdered)
}

# Find the count of MetaboLights studies being reused
MetaboLightsReuse <- function(x) {
    MetaboLightsC <- x[grepl("MTBLS", x)]
    MetaboLightsCount <- as.data.frame(table(MetaboLightsC))
    MetaboLightsCount <- MetaboLightsCount[order(-MetaboLightsCount$Freq), ]
    MetaboLightsCount$Repository <- rep("MetaboLights", nrow(MetaboLightsCount))
    colnames(MetaboLightsCount) <- c("Study", "Frequency", "Repository")
    return(MetaboLightsCount)
}

# Find the count of Metabolomics Workbench studies being reused
MWReuse <- function(x) {
    MetabolomicsWorkbenchC <- x[grepl("ST000", x)]
    MetabolomicsWorkbenchCount <- as.data.frame(table(MetabolomicsWorkbenchC))
    MetabolomicsWorkbenchCount <- MetabolomicsWorkbenchCount[order(-MetabolomicsWorkbenchCount$Freq), 
        ]
    MetabolomicsWorkbenchCount$Repository <- rep("Metabolomics Workbench", nrow(MetabolomicsWorkbenchCount))
    colnames(MetabolomicsWorkbenchCount) <- c("Study", "Frequency", "Repository")
    return(MetabolomicsWorkbenchCount)
}

# Find the count of GNPS studies being reused
GNPSReuse <- function(x) {
    GNPSC <- x[grepl("MSV", x)]
    GNPSCount <- as.data.frame(table(GNPSC))
    GNPSCount <- GNPSCount[order(-GNPSCount$Freq), ]
    GNPSCount$Repository <- rep("GNPS", nrow(GNPSCount))
    colnames(GNPSCount) <- c("Study", "Frequency", "Repository")
    return(GNPSCount)
}

# Combine the amount of reuse per repository into a single data frame
# Requires input of the number of GNPS, MetaboLights and Metabolomics
# Workbench studies which have not been reused
CombStudies <- function(Studylist, GNPS, ML, MW) {
    GNPSs <- GNPSReuse(Studylist)
    MLs <- MetaboLightsReuse(Studylist)
    MWs <- MWReuse(Studylist)
    AllStudies <- rbind(GNPSs, MLs, MWs)
    MLnr <- data.frame(Study = rep("Study", ML), Frequency = rep(0, ML), Repository = rep("MetaboLights"))
    # MW 609 - nrow(MetabolomicsWorkbenchCount) = 364
    MWnr <- data.frame(Study = rep("Study", MW), Frequency = rep(0, MW), Repository = rep("Metabolomics Workbench"))
    # GNPS 688 - nrow(GNPSCount) = 604
    GNPSnr <- data.frame(Study = rep("Study", GNPS), Frequency = rep(0, GNPS), 
        Repository = rep("GNPS"))
    # Combine into single dataframe
    AllStudies <- rbind(AllStudies, MLnr, MWnr, GNPSnr)
    AllStudies$Repository <- as.factor(AllStudies$Repository)
    levels(AllStudies$Repository) <- gsub(" ", "\n", levels(AllStudies$Repository))
    return(AllStudies)
}

# Create data frame of the percentage of reuse at each frequency
PerReuse <- function(x) {
    # The number of studies with each frequency of reuse
    Per0 <- nrow(x[x$Frequency == 0, ])/nrow(x) * 100
    Per1 <- nrow(x[x$Frequency == 1, ])/nrow(x) * 100
    Per2 <- nrow(x[x$Frequency == 2, ])/nrow(x) * 100
    Per3 <- nrow(x[x$Frequency == 3, ])/nrow(x) * 100
    Per4 <- nrow(x[x$Frequency == 4, ])/nrow(x) * 100
    df <- data.frame(FrequencyofReuse = c(0, 1, 2, 3, 4), Percentage = c(Per0, 
        Per1, Per2, Per3, Per4))
    return(df)
}

# Create a data frame of the reuse per repository and the reuse across
# repositories
PerReuseRepository <- function(x) {
    GNPS <- x[x$Repository == "GNPS", ]
    ML <- x[x$Repository == "MetaboLights", ]
    MW <- x[x$Repository == "Metabolomics\nWorkbench", ]
    GNPSPer <- PerReuse(GNPS)
    MLPer <- PerReuse(ML)
    MWPer <- PerReuse(MW)
    AllPer <- PerReuse(x)
    df = cbind(GNPSPer, MLPer[, 2], MWPer[, 2], AllPer[, 2])
    colnames(df) <- c("Frequency", "GNPS", "MetaboLights", "Metabolomics Workbench", 
        "All Studies")
    return(df)
}

# Percentage of studies published in each year that have been reused
# Requires input of list of studies that have been reused and all studies in
# repository
AgeReuse <- function(Repos, ReStud) {
    # Reused studies
    Reused <- Repos[Repos$StudyID %in% ReStud$Studies, ]
    # Not reused studies
    Not <- Repos[!Repos$StudyID %in% ReStud$Studies, ]
    # Publication Year of studies that have been reused
    ReusedAge <- as.data.frame(table(Reused$Year))
    # Publication Year of studies that have not been reused
    NotAge <- as.data.frame(table(Not$Year))
    df <- cbind(ReusedAge, NotAge[1:nrow(ReusedAge), 2])
    colnames(df) <- c("Year", "Reused", "Not")
    # Calculate percentage of studies reused per year
    df$Percentage <- round(df$Reused/(df$Reused + df$Not) * 100, 2)
    # Remove frequency columns
    df <- df[, -(2:3)]
    return(df)
}
```

## Reuse per Year

Calculate the number of publications that reuse data from each repository per year.

```r
# Load Reuse data
DataReuse <- read.csv("../data/DataReuseMetabolomics.csv", check.names = FALSE)
DataReuse$Year <- as.factor(DataReuse$Year)

# Convert data re-use per year to a factor (including 0 for the year 2012)
DataReuse$Year <- factor(DataReuse$Year, levels = c("2012", levels(DataReuse$Year)))
# Find the frequency of data reuse per year
DataReuseYear <- as.data.frame(table(DataReuse$Year))
colnames(DataReuseYear) <- c("Year", "Frequency")
DataReuseYear$Repository <- rep("All")

# Count number of reuses per repository MetaboLights
MLNo <- DataReuse[grepl("MetaboLights", DataReuse$Repository), ]
# Metabolomics Workbench
MWNo <- DataReuse[grepl("Metabolomics Workbench", DataReuse$Repository), ]
# GNPS
GNPSNo <- DataReuse[grepl("GNPS", DataReuse$Repository), ]

# Create contigency tables of the frequency of reuse per year adding a level
# for the year prior to any examples of reuse for plotting MetaboLights
MLYear <- as.data.frame(table(MLNo$Year))
colnames(MLYear) <- c("Year", "Frequency")
MLYear$Repository <- rep("MetaboLights")
MLYear$Year <- as.factor(MLYear$Year)
MLYear$Year <- factor(MLYear$Year, levels = c("2012", levels(MLYear$Year)))
# Metabolomics Workbench
MWYear <- as.data.frame(table(MWNo$Year))
colnames(MWYear) <- c("Year", "Frequency")
MWYear$Repository <- rep("Metabolomics Workbench")
MWYear$Year <- as.factor(MWYear$Year)
MWYear$Year <- factor(MWYear$Year, levels = c("2016", levels(MWYear$Year)))
# GNPS
GNPSYear <- as.data.frame(table(GNPSNo$Year))
colnames(GNPSYear) <- c("Year", "Frequency")
GNPSYear$Repository <- rep("GNPS")
GNPSYear$Year <- as.factor(GNPSYear$Year)
GNPSYear$Year <- factor(GNPSYear$Year, levels = c("2015", levels(GNPSYear$Year)))

# Recombine data into single data frame
Freqperyear <- rbind(DataReuseYear, MLYear, MWYear, GNPSYear)
Freqperyear$Repository <- as.factor(Freqperyear$Repository)

# Remove 2018 for ploting
Y2018 <- grepl("2018", Freqperyear$Year)
No2018 <- Freqperyear[!Y2018, ]
```

## Reuse per Study

Calculate the frequency of reuse of each study.


```r
# Find the total number of studies being reused (including reuse in Spicer
# et al. 2018)
TotalUniqueStudies <- ReusedStudies(DataReuse)

# As of 15/2/18 there are 329 MetaboLghts studies, 614 Metabolomics
# Workbench studies and 691 GNPS studies GNPS 691 - nrow(GNPSCount) = 606 ML
# 332 - nrow(MetaboLightsCount) = 187 MW 614 - nrow(MWCount) = 369
AllStudies <- CombStudies(TotalUniqueStudies, 606, 187, 369)

# Find the percentage reuse at each frequency across all repositories
Reuseper <- PerReuseRepository(AllStudies)
# melt data for ploting
Reuseperg <- melt(Reuseper, "Frequency")
colnames(Reuseperg) <- c("Frequency", "Repository", "Percentage")
Reuseperg$Frequency <- factor(Reuseperg$Frequency, levels = c("4", "3", "2", 
    "1", "0"))
levels(Reuseperg$Repository) <- gsub(" ", "\n", levels(Reuseperg$Repository))

# Find the total number of studies being reused (excluding reuse in Spicer
# et al. 2018)
TotalUniqueStudiesnoSpicer <- ReusedStudies(DataReuse[!grepl("Spicer", DataReuse$`Author(-s)`), 
    ])

# As of 15/2/18 there are 329 MetaboLghts studies, 614 Metabolomics
# Workbench studies and 691 GNPS studies GNPS 691 - nrow(GNPSCountnoSpicer)
# = 606 ML 329 - nrow(MetaboLightsCountnoSpicer) = 275 MW 614 -
# nrow(MWCountnoSpicer) = 592
AllStudiesnospicer <- CombStudies(TotalUniqueStudiesnoSpicer, 606, 275, 592)

# Find the percentage reuse at each frequency across all repositories
Reusepernospicer <- PerReuseRepository(AllStudiesnospicer)
# melt data for ploting
Reusepernospicerg <- melt(Reusepernospicer, "Frequency")
colnames(Reusepernospicerg) <- c("Frequency", "Repository", "Percentage")
Reusepernospicerg$Frequency <- factor(Reusepernospicerg$Frequency, levels = c("4", 
    "3", "2", "1", "0"))
levels(Reusepernospicerg$Repository) <- gsub(" ", "\n", levels(Reusepernospicerg$Repository))
```

## Same Authors

Studies were examined to see whether they reused data produced by the same set of authors, submitters or study owners. This is important as one of the most common reasons researchers cite for not sharing data is the fear that they will be able to generate less publications from their data, and other researchers will scoop them.

Percentage of studies that reuse data that share at least one author/submitter/study owner with the original study:

```r
# Calculate percentage of studies where the original authors have reused
# their own data
SameAuth <- round(sum(DataReuse$`Same Group`)/(nrow(DataReuse)) * 100, digits = 2)
SameAuth
```

```
## [1] 47.06
```

## Time until Reuse
Calculate the time between data being made public and its reuse. The average number of years is:

```r
# Load repository data MetaboLights
ML <- read.csv("../data/MLStudiesTime.csv", check.names = FALSE, stringsAsFactors = FALSE)
ML$Year <- format(as.Date(ML$StudyPublicationDate, "%d/%m/%y"), "%Y")
# Metabolomics Workbench
MW <- read.csv("../data/MWStudiesTime.csv", check.names = FALSE, stringsAsFactors = FALSE)
MW$Year <- format(as.Date(MW$ReleaseDate, "%d/%m/%Y"), "%Y")
# GNPS
GNPS <- read.csv("../data/GNPSStudiesTime.csv", check.names = FALSE, stringsAsFactors = FALSE)
GNPS$Year <- format(as.Date(GNPS$`Upload Date`, "%b. %d,%Y"), "%Y")

# Split DataReuse by study, but maintain details of reuse
UniqueStudiesReuse <- DataReuse %>% mutate(`Study(-ies)` = strsplit(as.character(`Study(-ies)`), 
    ",")) %>% unnest(`Study(-ies)`)
# rename colname `Study(-ies)` to StudyID for joining
colnames(UniqueStudiesReuse)[9] <- "StudyID"
# rename year to ReuseYear
colnames(UniqueStudiesReuse)[3] <- "ReuseYear"

# Create data frames of just StudyID and Year
MLSY <- as.data.frame(cbind(ML[, 1], ML[, 7]))
MWSY <- as.data.frame(cbind(MW[, 1], MW[, 10]))
GNPSSY <- as.data.frame(cbind(GNPS[, 2], GNPS[, 10]))
# join data frames
StIDY <- rbind(MLSY, MWSY, GNPSSY)
colnames(StIDY) <- c("StudyID", "Year")

# add year to studies
UniqueStudiesReuseYear <- left_join(UniqueStudiesReuse, StIDY, by = "StudyID")
# Convert years from factor to numeric
UniqueStudiesReuseYear$Year <- as.numeric(as.character(UniqueStudiesReuseYear$Year))
UniqueStudiesReuseYear$ReuseYear <- as.numeric(as.character(UniqueStudiesReuseYear$ReuseYear))
# Calculate time between data being published and reuse
UniqueStudiesReuseYear$Time2Reuse <- UniqueStudiesReuseYear$ReuseYear - UniqueStudiesReuseYear$Year

# Create contigency table
Time2Reuse <- as.data.frame(table(UniqueStudiesReuseYear$Time2Reuse))
colnames(Time2Reuse) <- c("Years", "Frequency")
Time2Reuse$Percentage <- Time2Reuse$Frequency/sum(Time2Reuse$Frequency) * 100

# Average time until reuse
AveTime <- round(mean(UniqueStudiesReuseYear$Time2Reuse, na.rm = T), 2)
AveTime
```

```
## [1] 1.8
```


## Time Data made Public

Calculate the time until data reuse per study in each repository.


```r
# Unique Studies that have been reused + Frequency
AllStudiesReused <- UniqueStudiesCount(DataReuse)

# Percentage reused Metabolights studies
MLAge <- AgeReuse(ML, AllStudiesReused)

# Percentage reused Metabolomics Workbench studies
MWAge <- AgeReuse(MW, AllStudiesReused)

# Percentage reused Metabolomics Workbench studies
GNPSAge <- AgeReuse(GNPS, AllStudiesReused)

# Create contigency tables of the frequency of reuse per year adding a level
# for the year prior to any examples of reuse for plotting MetaboLights
MLAge$Repository <- rep("MetaboLights")
MLAge$Year <- as.factor(MLAge$Year)
# Metabolomics Workbench
MWAge$Repository <- rep("Metabolomics Workbench")
MWAge$Year <- as.factor(MWAge$Year)
# GNPS
GNPSAge$Repository <- rep("GNPS")
GNPSAge$Year <- as.factor(GNPSAge$Year)

# Recombine data into single data frame
PerAge <- rbind(MLAge, MWAge, GNPSAge)
PerAge$Repository <- as.factor(PerAge$Repository)
# Add row showing 0 reuse for GNPS in 2017
GNPS2017 <- c("2017", 0, "GNPS")
PerAge <- rbind(PerAge, GNPS2017)
PerAge$Percentage <- as.numeric(PerAge$Percentage)

# Reuse excluding Spicer et al, (2017)
NoSpicerReused <- UniqueStudiesCount(DataReuse[!grepl("Spicer", DataReuse$`Author(-s)`), 
    ])

# Percentage reused Metabolights studies
MLNSAge <- AgeReuse(ML, NoSpicerReused)

# Percentage reused Metabolomics Workbench studies
MWNSAge <- AgeReuse(MW, NoSpicerReused)

# Percentage reused Metabolomics Workbench studies
GNPSNSAge <- AgeReuse(GNPS, NoSpicerReused)

# Create contigency tables of the frequency of reuse per year adding a level
# for the year prior to any examples of reuse for plotting MetaboLights
MLNSAge$Repository <- rep("MetaboLights")
MLNSAge$Year <- as.factor(MLNSAge$Year)
# Metabolomics Workbench
MWNSAge$Repository <- rep("Metabolomics Workbench")
MWNSAge$Year <- as.factor(MWNSAge$Year)
# GNPS
GNPSNSAge$Repository <- rep("GNPS")
GNPSNSAge$Year <- as.factor(GNPSNSAge$Year)

# Recombine data into single data frame
PerNSAge <- rbind(MLNSAge, MWNSAge, GNPSNSAge)
PerNSAge$Repository <- as.factor(PerNSAge$Repository)
# Add row showing 0 reuse for Metabolomics Workbench and GNPS in 2017
MW2017 <- c("2017", 0, "Metabolomics Workbench")
PerNSAge <- rbind(PerNSAge, MW2017, GNPS2017)
PerNSAge$Percentage <- as.numeric(PerNSAge$Percentage)
```

Calculate the number of studies released per repository per year.


```r
# Average deposition date of MetaboLights studies
MLSDT <- round(mean(as.numeric(ML$Year), na.rm = T), 2)

# Average deposition date of Metabolomics Workbench studies
MWSDT <- round(mean(as.numeric(MW$Year), na.rm = T), 2)

# Average deposition date of GNPS studies
GNPSDT <- round(mean(as.numeric(GNPS$Year), na.rm = T), 2)

# Average deposition date of all studies
Years <- cbind(ML$Year, MW$Year, GNPS$Year)
YearsDT <- round(mean(as.numeric(Years), na.rm = T), 2)

# Distribution of age of studies MetaboLights
MLSYear <- as.data.frame(table(ML$Year))
colnames(MLSYear) <- c("Year", "Frequency")
MLSYear$Repository <- rep("MetaboLights")
# Metabolomics Workbench
MWSYear <- as.data.frame(table(MW$Year))
colnames(MWSYear) <- c("Year", "Frequency")
MWSYear$Repository <- rep("Metabolomics Workbench")
# GNPS
GNPSSYear <- as.data.frame(table(GNPS$Year))
colnames(GNPSSYear) <- c("Year", "Frequency")
GNPSSYear$Repository <- rep("GNPS")
# All studies
AllSYears <- as.data.frame(table(Years))
colnames(AllSYears) <- c("Year", "Frequency")
AllSYears$Repository <- rep("All")

# Recombine data into single data frame
StudYear <- rbind(AllSYears, MLSYear, MWSYear, GNPSSYear)
StudYear$Repository <- as.factor(StudYear$Repository)

# Remove 2018, 2019 and 2050 for ploting
y2018 <- grepl("2018", StudYear$Year)
Till2017 <- StudYear[!y2018, ]
y2019 <- grepl("2019", Till2017$Year)
Till2017 <- Till2017[!y2019, ]
y2050 <- grepl("2050", Till2017$Year)
Till2017 <- Till2017[!y2050, ]

# Add zero years for plotting
MW0 <- data.frame(Year = 2012, Frequency = 0, Repository = "Metabolomics Workbench")
GNPS0 <- data.frame(Year = 2013, Frequency = 0, Repository = "GNPS")

# recombine data frames
Till2017$Year <- factor(Till2017$Year, levels = c("2011", levels(Till2017$Year)))
Till2017 <- rbind(Till2017, MW0, GNPS0)
```

# Tables

## Table 3.7. Studies that reuse public available metabolomics data



```r
Paper <- DataReuse$Link
kable(cbind(DataReuse[, 1:4], Paper), caption = "Studies that reuse public available metabolomics data, as of 15^th^ February 2018. The table shows articles that reuse publicly available data, the repository(-ies) that were the source of the data, the year the article was published and classification of how the data were re-used.") %>% 
    kable_styling(full_width = F, bootstrap_options = c("hover", "responsive"))
```

<table class="table table-hover table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Studies that reuse public available metabolomics data, as of 15^th^ February 2018. The table shows articles that reuse publicly available data, the repository(-ies) that were the source of the data, the year the article was published and classification of how the data were re-used.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Title </th>
   <th style="text-align:left;"> Repository </th>
   <th style="text-align:left;"> Year </th>
   <th style="text-align:left;"> Classification </th>
   <th style="text-align:left;"> Paper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Predicting Network Activity from High Throughput Metabolomics </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2013 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1371/journal.pcbi.1003123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> The Risa R/Bioconductor package: integrative data analysis from experimental metadata and back again. </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2014 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1186%2F1471-2105-15-S1-S11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PredRet: Prediction of Retention Time by Direct Mapping between Multiple Chromatographic Systems </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2015 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1021/acs.analchem.5b02287 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> The influence of scaling metabolomics data on model classification accuracy </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2015 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1007/s11306-014-0738-7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Joint Analysis of Dependent Features within Compound Spectra Can Improve Detection of Differential Features </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2015 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.3389/fbioe.2015.00129 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiNChE: A web tool and library for chemical enrichment analysis based on the ChEBI ontology </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2015 </td>
   <td style="text-align:left;"> Resource </td>
   <td style="text-align:left;"> https://doi.org/10.1186/s12859-015-0486-3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Approaches to sample size determination for multivariate data: Applications to PCA and PLS-DA of omics data </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1021/acs.jproteome.5b01029 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Galaxy-M: a Galaxy workflow for processing and analyzing direct infusion and liquid chromatography mass spectrometry-based metabolomics data </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1186/s13742-016-0115-8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Performance Evaluation and Online Realization of Data-driven Normalization Methods Used in LC/MS based Untargeted Metabolomics Analysis </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1038/srep38881 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Effect of Insulin Resistance on Monounsaturated Fatty Acid Levels: A Multi-cohort Non-targeted Metabolomics and Mendelian Randomization Study </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Biological studies </td>
   <td style="text-align:left;"> https://doi.org/10.1371/journal.pgen.1006379 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Non-targeted metabolomics combined with genetic analyses identifies bile acid synthesis and phospholipid metabolism as being associated with incident type 2 diabetes </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Biological studies </td>
   <td style="text-align:left;"> https://doi.org/10.1007/s00125-016-4041-1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Partial least squares with structured output for modelling the metabolomics data obtained from complex experimental designs: A study into the Y-block coding </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.3390/metabo6040038 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Dereplication of peptidic natural products through database search of mass spectra </td>
   <td style="text-align:left;"> GNPS </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1038/nchembio.2219 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DES-ncRNA: A knowledgebase for exploring information about human micro and long noncoding RNAs based on literature-mining </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Resource </td>
   <td style="text-align:left;"> https://doi.org/10.1080/15476286.2017.1312243 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DES-TOMATO: A Knowledge Exploration System Focused On Tomato Species </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Resource </td>
   <td style="text-align:left;"> https://doi.org/10.1038/s41598-017-05448-0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MsPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry-Based Fragmentation in Metabolomics </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1021/acs.analchem.6b04358 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NOREVA: normalization and evaluation of MS-based metabolomics data </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1093/nar/gkx449 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Untargeted metabolomics suffers from incomplete data analysis </td>
   <td style="text-align:left;"> MetaboLights, Metabolomics Workbench </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1007/s11306-017-1246-3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Joint Bounding of Peaks Across Samples Improves Differential Analysis in Mass Spectrometry-Based Metabolomics </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1021/acs.analchem.6b04719 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LiverWiki: a wiki-based database for human liver </td>
   <td style="text-align:left;"> MetaboLights, Metabolomics Workbench </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Resource </td>
   <td style="text-align:left;"> https://doi.org/10.1186/s12859-017-1852-0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mzML2ISA &amp; nmrML2ISA: generating enriched ISA-Tab metadata files from metabolomics XML data </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1093/bioinformatics/btx169 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mass Spectral Feature List Optimizer (MS-FLO): A Tool To Minimize False Positive Peak Reports in Untargeted Liquid Chromatography-Mass Spectroscopy (LC-MS) Data Processing </td>
   <td style="text-align:left;"> Metabolomics Workbench </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1021/acs.analchem.6b04372 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xMSannotator: an R package for network-based annotation of high-resolution metabolomics data </td>
   <td style="text-align:left;"> Metabolomics Workbench </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Software </td>
   <td style="text-align:left;"> https://doi.org/10.1021/acs.analchem.6b01214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Distribution based nearest neighbor imputation for truncated high dimensional data with applications to pre-clinical and clinical metabolomics studies </td>
   <td style="text-align:left;"> Metabolomics Workbench </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1186/s12859-017-1547-6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Significance estimation for large scale untargeted metabolomics annotations </td>
   <td style="text-align:left;"> GNPS </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1038/s41467-017-01318-5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proposal for a common nomenclature for fragment ions in mass spectra of lipids. </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Metadata </td>
   <td style="text-align:left;"> https://doi.org/10.1371/journal.pone.0188394 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Metadata analyser: Measuring metadata quality </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Metadata </td>
   <td style="text-align:left;"> https://doi.org/10.1007/978-3-319-60816-7_24 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Assessing Public Metabolomics Metadata, Towards Improving Quality. </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Metadata </td>
   <td style="text-align:left;"> https://doi.org/10.1515/jib-2017-0054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Molecular structures enumeration and virtual screening in the chemical space with RetroPath2.0. </td>
   <td style="text-align:left;"> MetaboLights </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Resource </td>
   <td style="text-align:left;"> https://doi.org/10.1186/s13321-017-0252-9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chemical Similarity Enrichment Analysis (ChemRICH) as alternative to biochemical pathway mapping for metabolomic datasets </td>
   <td style="text-align:left;"> Metabolomics Workbench </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Resource </td>
   <td style="text-align:left;"> https://doi.org/10.1038/s41598-017-15231-w </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Meta-mass shift chemical profiling of metabolomes from coral reefs. </td>
   <td style="text-align:left;"> GNPS </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1073/pnas.1710248114 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Compliance with minimum information guidelines in public metabolomics repositories </td>
   <td style="text-align:left;"> MetaboLights, Metabolomics Workbench, MetaPhen, MeRy-B </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:left;"> Metadata </td>
   <td style="text-align:left;"> https://doi.org/10.1038/sdata.2017.137 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Evaluation and comparison of bioinformatic tools for the enrichment analysis of metabolomics data </td>
   <td style="text-align:left;"> MetaboLights, Metabolomics Workbench </td>
   <td style="text-align:left;"> 2018 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1186/s12859-017-2006-0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Increased diversity of peptidic natural products revealed by modification-tolerant database search of mass spectra </td>
   <td style="text-align:left;"> GNPS </td>
   <td style="text-align:left;"> 2018 </td>
   <td style="text-align:left;"> Methods </td>
   <td style="text-align:left;"> https://doi.org/10.1038/s41564-017-0094-2 </td>
  </tr>
</tbody>
</table>

# Figures
Code that was used to generate raw figures. Figures were further processed in Adobe Illustrator.

## Figure 3.10. The number of articles that reuse metabolomics data over time. 
Data reuse across all repositories is shown in red, GNPS reuse is shown in blue, MetaboLights reuse is shown in green and Metabolomics Workbench is shown in purple. Some articles reuse data from multiple repositories (both MetaboLights and Metabolomics Workbench), so the total number of articles that reuse metabolomics data per year is not the sum of the articles that reuse data from each repository per year. The launch year of each repository is also highlighted.


```r
ggplot(No2018, aes(Year, Frequency, color=Repository, group=Repository, shape=Repository))  + 
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 18), linetype="dashed", color = "black", size=0.2) +
  geom_segment(aes(x = 2, y = 0, xend = 2, yend = 17), linetype="dashed", color = "black", size=0.2) +
  geom_segment(aes(x = 3, y = 0, xend = 3, yend = 18), linetype="dashed", color = "black", size=0.3) +
  theme_bw() +
  theme(
    legend.position="top",
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title.align=0.5) +
  scale_color_manual(values = c("#D55E00","#0072B2", "#009E73", "#CC79A7")) +
  guides(color = guide_legend(title.position = "top")) +
  annotate("text", x = 1, y = 19, label = "MetaboLights \n launched") +
  annotate("text", x = 2, y = 18.5, label = "Metabolomics \n Workbench \n launched") +
  annotate("text", x = 3, y = 19, label = "GNPS \n launched")
```

<img src="figs/reuseperyear-1.png" style="display: block; margin: auto;" />

## Figure 3.11A. The percentage of studies reused at each frequency including all studies

The percentage of studies reused at each frequency: 0, 1, 2, 3, or 4 times, as of 15^th^ February 2018, including reuse in all studies.

```r
ggplot(Reuseperg, aes(x=Repository, y=Percentage, fill=Frequency))  +
  geom_bar(colour="black", stat="identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_fill_manual(values = c('#0868ac','#43a2ca','#7bccc4',"#bae4bc",'#f0f9e8')) +
  ylab("Frequency of Reuse per Study (%)")+
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        legend.text.align = 0,
        #legend.text = element_text(face = "italic"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
```

<img src="figs/reuseperstudyall-1.png" style="display: block; margin: auto;" />

## Figure 3.11B. The percentage of studies reused at each frequency excluding Spicer *et al.* (2017)

The percentage of studies reused at each frequency: 0, 1, 2, 3, or 4 times, as of 15^th^ February 2018 excluding reuse by Spicer \textit{et al.} (2017).


```r
ggplot(Reusepernospicerg, aes(x=Repository, y=Percentage, fill=Frequency))  +
  geom_bar(colour="black", stat="identity", width = 0.7) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  scale_fill_manual(values = c('#0868ac','#43a2ca','#7bccc4',"#bae4bc",'#f0f9e8')) +
  ylab("Frequency of Reuse per Study (%)")+
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        legend.text.align = 0,
        #legend.text = element_text(face = "italic"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
```

<img src="figs/reuseperstudynospicer-1.png" style="display: block; margin: auto;" />

## Figure 3.12A. The frequency of metabolomics studies released per year

The total number of studies released per year is shown in red, GNPS studies are in blue, MetaboLights are in green and Metabolomics Workbench are in purple.


```r
ggplot(Till2017, aes(Year, Frequency, color=Repository, group=Repository, shape=Repository))  + 
  geom_line() +
  geom_point(size = 3) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
  theme_bw() +
  theme(
    legend.position="top",
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title.align=0.5) +
  scale_color_manual(values = c("#D55E00","#0072B2", "#009E73", "#CC79A7")) +
  guides(color = guide_legend(title.position = "top"))
```

<img src="figs/sAge-1.png" style="display: block; margin: auto;" />


## Figure 3.12B. The frequency of publicly available studies and time until data reuse


```r
ggplot(Time2Reuse, aes(Years, Percentage, group = 1))  + 
  geom_line() +
  geom_point() +
  ylab("Percentage of Reused Studies") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  theme_bw() +
  theme(
    legend.position="top",
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title.align=0.5)
```

<img src="figs/time2reuse-1.png" style="display: block; margin: auto;" />


## Supplementary Figure 1. Types of Reuse


```r
# Extract the frequency of types of reuse
TypeReuse <- as.data.frame(table(DataReuse$Classification))
colnames(TypeReuse) <- c("Classification", "Frequency")
TypeReuse$Percentage <- TypeReuse$Frequency/sum(TypeReuse$Frequency) * 100
TypeReuseOrdered <- TypeReuse[order(-TypeReuse$Percentage),]

colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)', 'rgb(114,147,203)')

plot_ly(TypeReuseOrdered, labels = ~Classification, values = ~Percentage, type = 'pie',
        textposition = 'inside',
        textinfo = 'label+percent',
        insidetextfont = list(color = '#FFFFFF'),
        hoverinfo = 'text',
        text = ~paste(Frequency, "Studies"),
        marker = list(colors = colors,
                      line = list(color = '#FFFFFF', width = 1)),
        #The 'pull' attribute can also be used to create space between the sectors
        showlegend = FALSE) %>%
  layout(title = 'Study Classification',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
```

<!--html_preserve--><div id="14e91a291ecd" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="14e91a291ecd">{"x":{"visdat":{"14e97e68eb31":["function () ","plotlyVisDat"]},"cur_data":"14e97e68eb31","attrs":{"14e97e68eb31":{"labels":{},"values":{},"textposition":"inside","textinfo":"label+percent","insidetextfont":{"color":"#FFFFFF"},"hoverinfo":"text","text":{},"marker":{"colors":["rgb(211,94,96)","rgb(128,133,133)","rgb(144,103,167)","rgb(171,104,87)","rgb(114,147,203)"],"line":{"color":"#FFFFFF","width":1}},"showlegend":false,"alpha":1,"sizes":[10,100],"type":"pie"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Study Classification","xaxis":{"showgrid":false,"zeroline":false,"showticklabels":false},"yaxis":{"showgrid":false,"zeroline":false,"showticklabels":false},"hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"data":[{"labels":["Methods","Software","Resource","Metadata","Biological studies"],"values":[35.2941176470588,29.4117647058824,17.6470588235294,11.7647058823529,5.88235294117647],"textposition":["inside","inside","inside","inside","inside"],"textinfo":"label+percent","insidetextfont":{"color":"#FFFFFF"},"hoverinfo":["text","text","text","text","text"],"text":["12 Studies","10 Studies","6 Studies","4 Studies","2 Studies"],"marker":{"fillcolor":"rgba(31,119,180,1)","color":"rgba(31,119,180,1)","colors":["rgb(211,94,96)","rgb(128,133,133)","rgb(144,103,167)","rgb(171,104,87)","rgb(114,147,203)"],"line":{"color":"#FFFFFF","width":1}},"showlegend":false,"type":"pie","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1}},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":{"render":[{"code":"function(el, x) { var ctConfig = crosstalk.var('plotlyCrosstalkOpts').set({\"on\":\"plotly_click\",\"persistent\":false,\"dynamic\":false,\"selectize\":false,\"opacityDim\":0.2,\"selected\":{\"opacity\":1}}); }","data":null}]}}</script><!--/html_preserve-->

## Supplementary Figure 2. The percentage of studies published per year that have been reused including all studies


```r
ggplot(PerAge, aes(Year, Percentage, color=Repository, group=Repository, shape=Repository))  + 
  geom_line() +
  geom_point(size = 3) +
  ylab("Reused Studies (%)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 49), linetype="dashed", color = "black", size=0.3) +
  geom_segment(aes(x = 2, y = 0, xend = 2, yend = 49), linetype="dashed", color = "black", size=0.15) +
  geom_segment(aes(x = 3, y = 0, xend = 3, yend = 22), linetype="dashed", color = "black", size=0.35) +
  theme_bw() +
  theme(
    legend.position="top",
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title.align=0.5) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#CC79A7")) +
  guides(color = guide_legend(title.position = "top")) +
  annotate("text", x = 1.5, y = 45, label = "MetaboLights \n launched") +
  annotate("text", x = 2.55, y = 42.5, label = "Metabolomics \n Workbench \n launched") +
  annotate("text", x = 3.35, y = 18, label = "GNPS \n launched")
```

<img src="figs/timeall-1.png" style="display: block; margin: auto;" />

## Supplementary Figure 3. The percentage of studies published per year that have been reused excluding Spicer *et al.* (2017)


```r
ggplot(PerNSAge, aes(Year, Percentage, color=Repository, group=Repository, shape=Repository))  + 
  geom_line() +
  geom_point(size = 3) +
  ylab("Reused Studies (%)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 49), linetype="dashed", color = "black", size=0.3) +
  geom_segment(aes(x = 2, y = 0, xend = 2, yend = 6), linetype="dashed", color = "black", size=0.1) +
  geom_segment(aes(x = 3, y = 0, xend = 3, yend = 22), linetype="dashed", color = "black", size=0.35) +
  theme_bw() +
  theme(
    legend.position="top",
    axis.text = element_text(colour = "black"),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    # Remove gridlines and borders
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title.align=0.5) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#CC79A7")) +
  guides(color = guide_legend(title.position = "top")) +
  annotate("text", x = 1, y = 56, label = "MetaboLights \n launched") +
  annotate("text", x = 2, y = 16, label = "Metabolomics \n Workbench \n launched") +
  annotate("text", x = 3, y = 30, label = "GNPS \n launched")
```

<img src="figs/timenospicer-1.png" style="display: block; margin: auto;" />

## Supplementary Table 1. The frequency of reuse per study

```r
ReuseStud <- UniqueStudiesCount(DataReuse)
colnames(ReuseStud) <- c("Study", "Frequency")
kable(ReuseStud, row.names = F) %>% kable_styling(full_width = F)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Study </th>
   <th style="text-align:right;"> Frequency </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MTBLS36 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS93 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS1 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS124 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS126 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS146 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS17 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS2 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS20 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS214 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS28 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS87 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS90 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000075 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000354 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078568 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078598 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078839 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078936 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079450 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS103 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS125 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS127 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS140 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS19 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS213 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS229 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS24 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS265 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS266 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS273 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS289 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS315 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS32 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS341 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS364 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS38 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS40 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS424 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS74 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS79 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS88 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000011 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000077 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000091 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000163 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000220 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000284 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000320 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000321 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000326 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000340 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000342 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000382 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000383 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000387 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000403 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A05001 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A06001 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX10 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX101 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX104 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX106 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX11 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX115 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX15 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX16 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX18 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX27 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX28 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX36 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX42 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX44 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX45 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX5 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX53 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX6 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX66 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX68 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX69 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX70 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX71 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX72 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX74 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX75 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX76 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX77 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX82 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX86 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX90 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX93 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX96 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX97 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX98 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MEX99 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078552 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078557 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078567 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078577 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078584 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078586 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078589 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078603 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078604 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078606 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078607 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078611 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078612 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078628 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078635 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078649 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078658 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078670 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078683 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078708 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078710 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078711 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078719 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078726 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078744 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078787 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078803 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078805 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078811 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078812 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078816 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078817 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078832 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078836 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078892 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078903 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078922 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078937 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078960 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000078993 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079029 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079040 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079050 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079069 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079091 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079098 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079104 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079105 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079146 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079243 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079329 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079339 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079341 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079344 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079356 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079398 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079416 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079421 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079447 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079558 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079573 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079581 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079598 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079651 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079652 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079679 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079758 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079760 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079772 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079773 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079777 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079778 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079787 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079808 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079813 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079825 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079838 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079888 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079905 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MSV000079907 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS100 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS102 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS104 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS105 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS11 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS117 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS12 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS123 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS129 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS13 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS137 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS143 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS144 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS150 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS152 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS154 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS155 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS156 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS157 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS160 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS161 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS162 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS163 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS164 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS166 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS169 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS172 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS173 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS174 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS176 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS177 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS178 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS188 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS189 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS191 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS197 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS198 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS200 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS216 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS217 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS218 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS219 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS22 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS225 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS226 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS227 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS228 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS23 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS233 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS234 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS237 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS240 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS242 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS243 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS247 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS253 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS26 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS263 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS264 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS267 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS27 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS277 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS279 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS280 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS282 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS298 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS30 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS307 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS31 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS313 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS320 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS327 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS33 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS337 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS338 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS34 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS345 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS35 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS354 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS358 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS368 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS374 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS376 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS378 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS385 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS39 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS394 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS4 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS404 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS41 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS414 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS419 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS42 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS422 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS427 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS45 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS46 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS47 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS52 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS57 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS59 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS61 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS67 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS7 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS71 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS72 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS81 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS92 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS95 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MTBLS96 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000001 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000002 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000004 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000009 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000010 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000013 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000016 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000026 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000027 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000028 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000029 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000030 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000031 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000032 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000033 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000034 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000035 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000036 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000037 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000038 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000039 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000041 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000042 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000043 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000044 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000045 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000046 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000047 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000056 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000058 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000061 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000062 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000063 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000065 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000069 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000070 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000071 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000074 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000076 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000081 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000082 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000083 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000087 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000089 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000090 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000092 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000093 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000095 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000096 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000099 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000104 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000105 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000106 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000110 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000111 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000113 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000114 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000115 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000121 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000122 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000133 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000134 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000135 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000137 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000138 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000140 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000142 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000144 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000145 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000146 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000147 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000149 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000150 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000153 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000154 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000158 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000159 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000160 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000161 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000164 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000166 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000168 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000169 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000171 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000172 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000176 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000182 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000188 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000189 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000192 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000193 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000194 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000195 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000196 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000198 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000199 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000201 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000202 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000203 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000206 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000207 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000212 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000213 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000215 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000218 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000221 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000222 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000223 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000224 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000225 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000226 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000230 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000231 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000232 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000233 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000236 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000242 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000245 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000246 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000248 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000249 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000250 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000254 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000255 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000257 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000259 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000261 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000270 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000272 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000273 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000274 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000276 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000278 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000279 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000282 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000283 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000285 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000287 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000288 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000291 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000292 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000293 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000295 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000296 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000298 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000299 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000301 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000302 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000303 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000304 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000306 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000310 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000311 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000313 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000314 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000315 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000316 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000317 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000318 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000319 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000322 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000324 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000325 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000327 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000329 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000330 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000331 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000332 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000336 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000337 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000338 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000341 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000344 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000346 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000352 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000355 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000356 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000367 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000368 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000369 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000370 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000371 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000374 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000375 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000376 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000379 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000380 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000381 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000385 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000386 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000388 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000389 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000390 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000391 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000392 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000396 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000397 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000398 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000404 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000412 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000413 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000419 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000421 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000422 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000425 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000426 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000427 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000428 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000432 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000433 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000434 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000435 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000438 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000439 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000440 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000442 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000443 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000445 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000450 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000451 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000452 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000465 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000477 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000483 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000502 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000510 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000539 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000542 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ST000543 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>
