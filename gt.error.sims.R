#!/usr/bin/env Rscript

#created by Nick Sard
#created on 7/14/2015

#This script was created to simulate a data set and assign parents to offspring with an estimated genotyping error rate.
#assignments were used to calculate the expected proportion of offspring with one parent assigned to them.

#loading libraries
library(dplyr)

#Note: This script was originally run on a Linux based computer cluster.

#Loading in source scripts that are used simulate genotypic data, and to assign parents to offspring
#These scripts are modified from the original functions written in R library SOLOMON.
source("gp.analysis.R")

#using the 2010 parents as an example
my.parents <- 374

#first simulate some allele frequencies and makes a fake dataset
alfs <- read.table("alfs.txt", header=T,stringsAsFactors=F)
all.gts <- sim.gt.data(afreq = alfs,Nparents = my.parents,Noffs_perpair = 5,error = 0.02,Nunrelated = 0,Nsibs = 0,write2file = F)

#seperating out parents and kids
moms <- all.gts[grepl(pattern = "Mom",x = all.gts$IDs),-1]
dads <- all.gts[grepl(pattern = "Dad",x = all.gts$IDs),-1]
off <- all.gts[grepl(pattern = "Offspring",x = all.gts$IDs),-1]

#Assigning kids to potential moms
mom.filename <- paste0(paste("bayes.output.moms",Sys.time(),Sys.getpid(), sep = "_"),".txt")
no.par.bayes.ns(adults = moms,offspring = off,reps = 1000,Ngts = 5e7)
file.rename(from = "Output_SOLOMON_Posterior_Probabilities.txt",to = mom.filename)

#Assigning kids to potential dads
dad.filename <- paste0(paste("bayes.output.dads",Sys.time(),Sys.getpid(), sep = "_"),".txt")
no.par.bayes.ns(adults = dads,offspring = off,reps = 1000,Ngts = 5e7)
file.rename(from = "Output_SOLOMON_Posterior_Probabilities.txt",to = dad.filename)

#post processing for each
output1 <- read.table(mom.filename,header = T,sep = "\t")
colnames(output1) <- c("parent","offspring","mm","prob")
output1 <- filter(output1, mm <=1 & prob < 0.05)

output2 <- read.table(dad.filename, header = T,sep = "\t")
colnames(output2) <- c("parent","offspring","mm","prob")
output2 <- filter(output2, mm <=1 & prob < 0.05)

#Assembling the genetic pedigree
ped <- data.frame(off = off[,1],
                  mom = "UK",
                  dad = "UK", stringsAsFactors=F)

#Putting pedigree together

i <- NULL
for(i in 1:nrow(ped)){
  
  #for mom side of things
  tmp <- output1[output1$offspring == ped$off[i],]
  if(nrow(tmp) == 1){
    ped$mom[i] = as.character(tmp$parent)
  } else if(nrow(tmp) > 1) {
    ped$mom[i] = "duplicated"
  } else {
    ped$mom[i] = "UK"
  } #end of mom if statements

  #for dad side of things
  tmp <- output2[output2$offspring == ped$off[i],]
  if(nrow(tmp) == 1){
    ped$dad[i] = as.character(tmp$parent)
  } else if(nrow(tmp) > 1) {
   ped$dad[i] = "duplicated"
  } else {
   ped$dad[i] = "UK"
  } # end of dad if statements
} #end of ped for loop

#making a column describing the assignment type
ped$type <- "none"
ped$type[ped$mom != "UK" & ped$dad != "UK"] <- "both"
ped$type[ped$mom != "UK" & ped$dad == "UK"] <- "mom.only"
ped$type[ped$mom == "UK" & ped$dad != "UK"] <- "dad.only"
ped$type[ped$mom == "UK" & ped$dad == "UK"] <- "none"

#Summarizing data in a table
df <- data.frame(t(table(ped$type)))[,-1]
df$prop <- df$Freq/nrow(off)

df1 <- data.frame(type = c("both","mom.only","dad.only","none"),
		  freq = 0,
		  prop = 0,
		  run = 0, stringsAsFactors=F)

for(a in 1:nrow(df1)){
	if(nrow(df[df$Var2 == df1$type[a],])==1){ 
	   df1$freq[a] <- df[df$Var2 == df1$type[a],2]
	   df1$prop[a] <- df[df$Var2 == df1$type[a],3]
	}
}

#removing some files I don't want
#file.remove(mom.filename)
#file.remove(dad.filename)

#Writing summary to standard out
write.table(df1, "", sep="\t", quote=F, row.names=F)

#fin!