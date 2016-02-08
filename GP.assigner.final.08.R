#originally created on: September 30, 2015
#by: Nick Sard

#ABOUT: This script was written to assign parent pairs to grandoffspring
#Assignments must match all all loci
#Parents and offspring should match at all loci, if not a warning is invoked

#loading in libraries
library(dplyr)
library(reshape2)

#loading my own functions
# - none -

#setting working directory and reading in data
setwd("C:/Users/Nick.Sard/Documents/Research/R/Github/Grandparentage/")
list.files()

#reading in grandparents
gp1 <- read.table(file = "grandmothers.2008.txt",header = T,stringsAsFactors = F)
gp2 <- read.table(file = "grandfathers.2008.txt",header = T,stringsAsFactors = F)

#reading in parent and offpsring - assumes parent offspring assignments are on the same rows in the two data frames
par <- read.table(file = "grandkids.mothers.2012.txt",header = T,stringsAsFactors = F)
off <- read.table(file = "grandkids.missing.dads.2012.txt",header = T,stringsAsFactors = F)


head(par)
head(off)
head(gp1)
head(gp2)

#fixing the first name
names(par)[1] <- "Sample.Name"
names(off)[1] <- "Sample.Name"

### Step one - remove all the alleles the offspring matches to the parent

#recording the number of parents
num.par <- nrow(par)

#melting the dataframe
par1 <- par[duplicated(par$Sample.Name)==F,]
par1 <- melt(par1,"Sample.Name")
par1$variable <- as.character(par1$variable)
par1 <- arrange(par1, Sample.Name)

#getting the loci all the came name because loci are have two names E.g (Locus_1.1, Locus_1.2 - one for each allele)
loci <- unique(par1$variable)
locus.locs <- seq(from=1,to = length(loci)-1,by = 2)
loci <- loci[locus.locs]
my.loci <- loci
loci.2 <- rep(rep(x = loci,each=2),times=length(unique(par1$Sample.Name)))
par1$variable <- loci.2

#fixing the colnames to make the more understandable
colnames(par1) <- c("Sample.Name","loci","allele")
head(par1)

#recording the number of offspring
num.off <- nrow(off)

#melting the dataframe
off1 <- off[duplicated(off$Sample.Name)==F,]
off1 <- melt(off1,"Sample.Name")
off1$variable <- as.character(off1$variable)
off1 <- arrange(off1, Sample.Name)

#getting the loci all the came name because loci are have two names E.g (Locus_1.1, Locus_1.2 - one for each allele)
loci <- unique(off1$variable)
locus.locs <- seq(from=1,to = length(loci)-1,by = 2)
loci <- loci[locus.locs]
num.loci <- length(loci)
loci.2 <- rep(rep(x = loci,each=2),times=length(unique(off1$Sample.Name)))
off1$variable <- loci.2

#fixing the colnames to make the more understandable
colnames(off1) <- c("Sample.Name","loci","allele")
head(off1)

if(num.off != num.par) {stop("Number of parents and offspring do not match")}

#now for the removing of the alleles
head(off)
head(par)


#this for loop takes the each parent offspring pair (POP) and compares their alleles at each locus
#depending on how many alleles match it will either leave the one allele that doesn match or
#it will leave both alelles if they both match
i <- NULL
ep.offspring2 <- NULL
pops2check_same.gt <- NULL
pops2check_mm <- NULL
for(i in 1:num.off){
  
  #this assumes that each parent and offspring are the same row of the two datasets (par and off)
  #it grabs their genotypes in long form in par1 and off1
  each.parent <- par1[par1$Sample.Name == par[i,1],]
  each.offspring <- off1[off1$Sample.Name == off[i,1],]
  
  #now for each locus its compares the alleles
  j <- NULL
  keep.alleles <- NULL
  for(j in 1:length(loci)){

    #getting the alleles for the given locus for each POP
    ep.locus <- each.parent[each.parent$loci == loci[j],]
    eo.locus <-   each.offspring[each.offspring$loci == loci[j],]
    
    #figuring out which alleles they match at using offspring alleles first
    allele.logic <- eo.locus[order(eo.locus$allele),3] %in% ep.locus[order(ep.locus$allele),3]
    
    #counting how many alleles they match from this prospective
    checker <- length(allele.logic[allele.logic == T])
    
    #making suring that match all loci - if now the program dies
    if(checker == 0){
      warning(paste(par[i,1],"and",off[i,1],"do not match at all loci.","Check locus",loci[j], "See pops2check_mm"))
      pops2check_mm <- c(pops2check_mm,paste(par[i,1],off[i,1],loci[j], sep="++"))
    }
    
    #figuring out which alleles they match at using the parents alleles first
    allele.logic2 <-  ep.locus[order(ep.locus$allele),3] %in% eo.locus[order(eo.locus$allele),3] 
    
    #making sure they match at least one allele
    checker2 <- length(allele.logic2[allele.logic2 == T])
  
    #if they have the same genotype then I keep both alleles
    #this builds in flexibility down the line when I assign the offspring
    #grandparent pairs because they can match either allele
    if(checker == 2 & checker2 == 2){
      warning(paste(par[i,1],"and",off[i,1],"match at both alleles, which may be an issue in the grandparentage assignments.","See pops2check_same.gt"))
      keep.alleles <- c(keep.alleles,allele.logic)
      pops2check_same.gt <- c(pops2check_same.gt,paste(par[i,1],off[i,1],loci[j], sep="++"))
    }
    
    #if the offsprings is homozygous and the parent is heterozygous and the match
    #at one allele then I just remove one of the homozygous alleles
    if(checker == 2 & checker2 == 1){
      rand.loc <- sample(c(1,2),size = 1)
      allele.logic[rand.loc] <- F
      allele.logic <- ifelse(allele.logic == T,F,T)
      keep.alleles <- c(keep.alleles,allele.logic)
    }

    #if offspring is heterozygous and the parent is homosygous and they match
    #at one allele then I switch the logic so that I select the allele they
    #dont match at
    if(checker ==  1 & checker2 == 2){
      #switching that logic beacuse I only want one of them
      allele.logic <- ifelse(allele.logic == T,F,T)
      keep.alleles <- c(keep.alleles,allele.logic)
    }
    
    #if offspring is heterozygous and the parent is heteroygous and they match
    #at one allele then I switch the logic so that I select the allele they
    #dont match at
    if(checker ==  1 & checker2 == 1){
      #switching that logic beacuse I only want one of them
      allele.logic <- ifelse(allele.logic == T,F,T)
      keep.alleles <- c(keep.alleles,allele.logic)
    }

  }
  
  #adding the offspring with their alleles that could match a grandparent pair to the outfile
  ep.offspring2 <- rbind(ep.offspring2,each.offspring[keep.alleles,])
}
warnings()

#making that output into something easier to read/interpret/work with downstream
pops2check_mm <- data.frame(pops2check_mm)
pops2check_mm$parent <- lapply(strsplit(as.character(pops2check_mm$pops2check_mm), "\\++"), "[", 1)
pops2check_mm$offspring <- lapply(strsplit(as.character(pops2check_mm$pops2check_mm), "\\++"), "[", 2)
pops2check_mm$locus <- lapply(strsplit(as.character(pops2check_mm$pops2check_mm), "\\++"), "[", 3)
pops2check_mm$pops2check_mm <- NULL
head(pops2check_mm)


pops2check_same.gt <- data.frame(pops2check_same.gt)
pops2check_same.gt$parent <- lapply(strsplit(as.character(pops2check_same.gt$pops2check_same.gt), "\\++"), "[", 1)
pops2check_same.gt$offspring <- lapply(strsplit(as.character(pops2check_same.gt$pops2check_same.gt), "\\++"), "[", 2)
pops2check_same.gt$locus <- lapply(strsplit(as.character(pops2check_same.gt$pops2check_same.gt), "\\++"), "[", 3)
pops2check_same.gt$pops2check_same.gt <- NULL
head(pops2check_same.gt)


#now these are the alleles at each locus that I will compare to each grandparent parent pair
off.names <- unique(ep.offspring2$Sample.Name)

# now getting grandparent alleles in the same format

#making sure that each grandparent pair is only found once in the datasets
keepers <- duplicated(paste0(gp1$Sample.Name,gp2$Sample.Name))
gp1 <- gp1[keepers==F,]
gp2 <- gp2[keepers==F,]

#recording the number of grandparent1
num.gp1 <- nrow(gp1)
gp1$order <- 1:num.gp1

#melting the dataframe
gp1.1 <- melt(gp1,c("Sample.Name","order"))
gp1.1$variable <- as.character(gp1.1$variable)
gp1.1 <- arrange(gp1.1, Sample.Name,variable)
gp1.1$variable <- gsub('\\.1','',gp1.1$variable)
head(gp1.1)

#getting the loci all the came name because loci are have two names E.g (Locus_1.1, Locus_1.2 - one for each allele)
loci <- unique(gp1.1$variable)
locus.locs <- seq(from=1,to = length(loci)-1,by = 2)
loci <- loci[locus.locs]

#fixing the colnames to make the more understandable
colnames(gp1.1) <- c("Sample.Name","order","loci","allele")

#recording the number of grandparent1
num.gp2 <- nrow(gp2)
gp2$order <- 1:num.gp2

#melting the dataframe
gp2.1 <- melt(gp2,c("Sample.Name","order"))
gp2.1$variable <- as.character(gp2.1$variable)
gp2.1 <- arrange(gp2.1, Sample.Name)
gp2.1$variable <- gsub('\\.1','',gp2.1$variable)
head(gp2.1)

#fixing the colnames to make the more understandable
colnames(gp2.1) <- c("Sample.Name","order","loci","allele")

#adding a new columns to each
gp1.1$type <- "GP1"
gp2.1$type <- "GP2"

#now combining to grandparent pairs (gps)
gps <- rbind(gp1.1,gp2.1)
gps <- gps[order(gps$order),]

#getting the number of grandparent pairs
num.gp.pairs <- unique(gps$order)

#making in a fake dataset to be filled
output <- data.frame(GP1 = NA, GP2 = NA,Off = NA,gps = 0,epoff = 0,locus = 0,counter = 0)

#I know this is ugly because its three for loops, but I just had to get this done quickly
#The idea is that it compares each grandparent pair to each possible grandoffspring at each locus
#for loop "i" grabs the putative grandparent pair
#for loop "j" grabs the putative grandoffspring
#for loop "k" compares the grandparent pair and the grandoffspring at each locus to see if they 
#   match at least one allele. If so it adds their names to a outfile
i <-1
i <- NULL
for(i in num.gp.pairs){
  
  #getting the putative grandparent genotypes
  gps1 <- gps[gps$order == i,]

  #making sure there are only 44 alleles in there, if not the program is killed
  if(nrow(gps1)>(4*num.loci)){stop("Check GPS",i,"- there may be too many loci or GPs in there")}

  j <- NULL
  for(j in 1:length(off.names)){
    
    #grabbing the putative grandoffspring
    ep.off3 <- ep.offspring2[ep.offspring2$Sample.Name == off.names[j],]

    #now comparing the grandparent pair and the grandoffspring at all loci
    k <- NULL
    counter <- 0
    for(k in 1:num.loci){
      
      #grabing the given locus for the grandparent pair and the grandoffspring
      egp.locus <- gps1[gps1$loci == my.loci[k],]
      eo.locus <-   ep.off3[ep.off3$loci == my.loci[k],]

      #figuring out if they match at any alleles
      allele.logic <- egp.locus[,4] %in% eo.locus[,3]

      #making sure they only match at one allele by 
      #counint how many TRUEs there are
      checker <- length(allele.logic[allele.logic == T])
      
      #this break will break out of the loop early if they mismatch at a single locus
      if(checker == 0){break}
      
      #if the grandparent pair and the grandoffspring match at atleast one allele then the counter adds one
      if(checker >  0){
      counter = counter+1
      }
      
    }
    
    #checking to make sure they match at all loci
    #if they do then their names get added to an outfile
    #along with some other information, I used for troubleshooting
    if(counter >= num.loci){
      
      #filling there information into a dataframe
      output1 <- data.frame(GP1 = unique(gps1[gps1$type == "GP1",1]),
                           GP2 = unique(gps1[gps1$type == "GP2",1]),
                           Off = unique(ep.off3$Sample.Name),
                           gps = i,
                           epoff = j,
                           locus = k,
                           counter = counter)
      
      #adding it to the outfile
      output <- rbind(output,output1)
    }
  }
}

#removing that first column
output <- output[-1,]
head(output)

out1 <- output[,1:3]
head(out1)

#writing to file
write.table(x = out1,file = "Output/Gtrios.08.txt",append = F,sep = "\t",quote = F,row.names = F,col.names = T)

#fin!