#Generating the a priori graph for the windows from known chromosomes
#Halie Rando, November 2017, rando2@illinois.edu

#load and clean cnv data (just need the original output with the # reads mapping in each window)
#I used rstudio to load this file as a dataframe cnv10000
cnv_scaff <- gsub("scaffold", "", cnv10000$chromosome)
#I'm using a cutoff of the largest 12625 scaffs because we have 676878 and that's where they get smaller than 1 kb
cnv_isLong <- as.numeric(cnv_scaff)<=12625 
cnv_crop <- data.frame(cnv10000, isLong=cnv_isLong)
#now adjust to allow only large scaffs & a min of 100 reads for data quality
cnv_crop <- cnv_crop[which(cnv_crop$isLong == TRUE & 
                             (cnv_crop$test + cnv_crop$ref)>100),1:6]
#now calculate the % of reads found in the male (test)
cnv_restore <- data.frame(cnv_crop, perc.male=cnv_crop$test/(cnv_crop$test + cnv_crop$ref))
cnv_restore <- data.frame(cnv_restore, perc.male.norm=(cnv_restore$perc.male-
                                                         mean(cnv_restore$perc.male))/sd(cnv_restore$perc.male))
hist(cnv_restore$perc.male) #visualize to check ~50% (or other baserate if known)

#load & clean ygs data
#reads ygs output in (I cleaned it up a bit first in vi to add headers, etc) as df ygs_forR
names(ygs_forR)[1]<- "chromosome"
ygs_scaff <- gsub("scaffold", "", ygs_forR$chromosome)
ygs_isLong <- as.numeric(ygs_scaff)<=12625
ygs_valid <- ygs_forR[which(ygs_isLong==TRUE),] #make sure we only have top 12625 scaff

ygs_valid$P_VSC_UK <- as.numeric(ygs_valid$P_VSC_UK) #change data type
ygs_valid <- ygs_valid[which(is.na(ygs_valid$P_VSC_UK)==FALSE),] #drop NA
ygs_norm <- (ygs_valid$P_VSC_UK - mean(ygs_valid$P_VSC_UK))/sd(ygs_valid$P_VSC_UK)
ygs_valid <- data.frame(ygs_valid, P_VSC_UK.norm=(ygs_valid$P_VSC_UK - mean(ygs_valid$P_VSC_UK))/sd(ygs_valid$P_VSC_UK))

hist(ygs_valid$P_VSC_UK.norm) #check normalized data (mine is bi or tri modal)

#merge two datasets
combined <- (merge(cnv_restore, ygs_valid, by = 'chromosome')) 
#this is a little slow

#pull data from your predictions. I made this in Excel & it's a table of the scaffolds
#with their predicted location (X, Y, unknown, autosome) and then the color I want that
#type plotted as (red, blue, grey, black)
#apriori_colors is a dataframe made from this Excel table reads in as a csv
defined_apriori <- apriori_colors[which(apriori_colors$type != "unknown"),] #drop unknown

#this loop grabs the windows falling onto already-studied scaffolds and assigns them
#the correct color. Loops in R are slow, so it takes about 1-2 minutes for my data.
combined_apriori <- data.frame()
color <- c()
for(r in 1:nrow(defined_apriori)){ 
  scaff <- as.character(defined_apriori[r,1])
  temp <- combined[which(combined$chromosome == scaff),]
  combined_apriori <- rbind(combined_apriori, temp)
  color <- c(color, rep(as.character(defined_apriori[r,3]),nrow(temp)))
}

#Two basically identical plots, one with normalized and one with raw data
plot(combined_apriori$perc.male.norm, combined_apriori$P_VSC_UK.norm,col=color,
     main="Putative Chromosomal Trends", ylab = "% unmatched valid single-copy k-mers (normalized)", 
     xlab="% reads mapped to males (normalized)")
#Above plot shows the distribution of the scaffolds as predicted in genome paper
plot(combined_apriori$perc.male, combined_apriori$P_VSC_UK,col=color,
     main="Putative Chromosomal Trends", ylab = "% unmatched valid single-copy k-mers", 
     xlab="% reads mapped to males")
legend("bottomright", #this works because it doesn't cover any points
       c("Autosome","X", "Y"), # legend labels
       pch=1, # open dots
       col=c("black","red","blue")) # colors of dots