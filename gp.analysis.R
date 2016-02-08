####################
### sim.afreqs() ###
####################

#Description
#This fuction is a modified version of Mark Christie's script to create simulate genotype allele frequencies.
#Note that there are high and low frequencies for each locus, and the alleles and their associated frequencies
#are exactly the same. 

#Input Parameters:
#nlcoi - define the number of loci you want
#n.alleles - define the number of alleles you want


sim.afreqs <- function(nloci=10,n.alleles=10, simple=F){
  
  #making sure nloci and n.alleles are numeric
  nloci <- as.numeric(nloci)
  n.alleles <- as.numeric(n.alleles)
  
  #getting the alleles frequencies for each allele
  a <- c(1:n.alleles)
  b <- rev(a)
  c <- b/a
  d <- sum(c)
  
  #now making actual allele frequencies
  freqs <- c/d
  
  #setting the lowest allele to an arbitrary number
  lowest.allele <- 140
  
  #making afreqs
  
  #getting loci names
  loci <- paste0("Locus.",1:nloci)
  loci <- rep(loci, each=n.alleles)
  
  #for alternative output
  loci.1 <- rep(1:nloci, each=n.alleles)
  
  #getting the allele names
  alleles2 <- cbind(seq(lowest.allele,lowest.allele+n.alleles-1,1),freqs)
  alleles <- rep(alleles2[,1],nloci)
  
  
  #getting their frequencies
  freqs  <-  rep(alleles2[,2],nloci)
  
  #returning finished product
  if(simple == F){
    
    afreqs <- data.frame(loci,alleles,freqs, stringsAsFactors=F)
    afreqs
    return(afreqs)  
    
  } else {
    
    afreqs <- data.frame(loci.1,freqs, stringsAsFactors=F)
    afreqs
    return(afreqs)  
  }
  
  
} # end of sim.afreqs

#####################
### sim.gt.data() ###
#####################

#Description
#This fuction is a modified version of Mark Christie's script to create parents and offspring.
#Note that one can create random genotyping error, add unrelated indviduals, and add full siblings in the mix

#Input Parameters:
#afreq - this is a data frame with two rows in it, the first is the locus, and the second is the allele frequencies
#Nparents - defining thenumber of parents one wants
#Noffs_perpair - defining the number of offspring per parent one wants
#error - defining the amount of random error in the dataset
#Nunrelated - defining the number of unrelated individuals wanted in the dataset
#Nsibs - defining the number of full siblings one wants

#Create Simulated Data Sets==========================================#

sim.gt.data <-  function(afreq, Nparents= 2, Noffs_perpair= 1, error= 0, Nunrelated= 0, Nsibs=0, write2file = F){
  
  #defining SimOffs.txt name
  off.filename <- paste("Off",Sys.time(), Sys.getpid(), sep = "_")
  
  #making sure everything is numeric
  Nparents <- as.numeric(Nparents)
  Noffs_perpair <- as.numeric(Noffs_perpair)
  error <- as.numeric(error)
  Nunrelated <- as.numeric(Nunrelated)
  
  #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
  Nadults <- Nparents*2             
  
  #getting afreqs input
  afreqs <- afreq
  afreqs
  
  #making a function to simulate genotypes for parents
  
  OUT <- NULL
  sims <- function(sims){
    
    #table allele frequencies
    #getting just the frequencies 
    alleles2 <- afreqs[afreqs[,1] == loci[i],]
    alleles2
    
    #changing generic name to three character genotypes
    alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                  
    alleles3
    
    #create homozygote allele frequencies
    homos <- (alleles3[,2])^2                                                          
    homos
    homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
    homos2
    
    #create heterozygote allele frequencies
    
    #first create all possible combinations of freqs
    hets <- t(combn(alleles3[,2],2))                                                   
    hets
    
    #now make expected freqs
    hetfreq <- 2*(hets[,1]*hets[,2])
    hetfreq
    
    #create heterozygote allele names
    #now getting the het. genotype combinations
    hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
    hetvals
    
    #combing them
    hets2 <- cbind(hetvals,hetfreq)
    hets2
    
    #combine hets and homos and create genotypes
    gfreqs <- rbind(hets2,homos2)                                                      
    gfreqs
    
    #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
    n <- 1000000                                    
    
    #create genotypes(by coloumn, 1 for each allele)
    #for the first allele
    
    gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
    gfreqs1
    
    #now the second
    gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
    gfreqs2
    
    #combining them
    gtypes <- cbind(gfreqs1,gfreqs2)
    head(gtypes)
    
    #mixing them up
    gtypes <- gtypes[sample(1:length(gtypes[,1]),replace= F),]
    head(gtypes)
    
    #now getting a random sample of the gentoypes for the parents
    sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
    sg1
    
    OUT<<-cbind(OUT,sg1)
    
  } # end of sim function
  
  #defining the number of loci there are
  loci.length <- length(unique(afreqs[,1]))
  loci <- unique(afreqs[,1])
  loci
  
  #doing the simulation process for each locus
  for(i in 1:length(loci)){
    lapply(1,sims)
  } # end for loop
  OUT
  warnings()
  OUT
  #saving OUT into object called parents
  parents <- OUT
  head(parents)
  
  #next making a matrix to add to parents so that gts in every two columns don't overlap
  #first making a matrix the same size of parents and filling it ever two columns with varying
  #numbers in the 1000s
  c <- c(1:(ncol(OUT)))
  odd <- 2*(unique(round(((c-2))/2)))+1
  l <- length(odd) * 1000
  codes <- seq(from= 1,to= l,by= 1000)
  cols <- sort(rep(codes,2))-1
  Anumbs <- matrix(data= cols,nrow= Nadults,ncol= ncol(OUT),byrow= T)
  Anumbs
  
  #now adding it to parents
  parents <- as.numeric(parents)+Anumbs
  parents
  
  #create full sib families (go down the list in pair) ======================#
  OUT2 <- NULL
  sims <- function(sims){
    
    #first grabbing the genotypes of the parents
    p1 <- parents[i,]
    p1
    p2 <- parents[i+1,]
    p2
    #defining the order of teh alleles
    als <- rep(1:2,length(p1)/2)
    als
    
    #defining the number of offspring that need to be full sibs
    Noffs <- Noffs_perpair                                                          
    Noffs
    
    #using a for loop to making full sib genotypes
    OUT2 <- NULL
    for (b in 1:Noffs){
      
      #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
      #sampling alleles from parent one
      pos1 <- sample(als,length(p1)/2,replace <- TRUE)                                      
      pos1
      
      #sampling alleles from parent two
      pos2 <- sample(als,length(p1)/2,replace <- TRUE)
      pos2
      
      #getting the position for those alleles from parent one
      pos11 <- pos1+(seq(0,(length(p1)-1),2))
      pos11
      
      #getting the position for those alleles from parent two
      pos22 <- pos2+(seq(0,(length(p2)-1),2))
      pos22
      
      #getting those alleles from parent one
      o1 <- p1[pos11]
      o1
      #getting those alleles from parent two
      o2 <- p2[pos22]
      o2
      
      #putting them together to make the offsprings genotype
      o3 <- sort(c(o1,o2))
      o3
      
      #identifying what parents the offspring came fromm
      o3 <- t(c(i,o3))
      o3
      
      #writting to file and appending
      write.table(o3,file = off.filename ,row.names=FALSE,col.names=F,sep="\t",append=T)
    } #end of full sib for loop
  } # end of new sims function
  
  #defining which parents are going to have full sibs
  my.picks <- seq(from=1,to=Nadults-1, by=2)
  my.picks
  
  #making full sib gts for each parent pair
  for(i in my.picks){ 
    lapply(i,sims)
  } # end of full sib genotyping loop
  
  #removing the 1000s from the gts
  Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)                                   
  parents <- as.numeric(parents)-Anumbs
  parents
  
  #getting parent genotypes
  #first the moms
  Moms <- parents[seq(from=2,to= Nadults,by= 2),]
  Moms
  #now the dads
  Dads <- parents[seq(from=1,to= Nadults,by= 2),]
  Dads
  
  #reading SimOffs.txt back in and removing the id column
  Offs2 <- read.table(off.filename, header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Offs  <-  as.matrix(Offs2[,-1])
  
  #removing issue with 1000s
  Anumbs <- matrix(cols,length(Offs[,1]),ncol(OUT),byrow = T)
  Offs <- as.numeric(Offs)-Anumbs
  Offs
  
  #now making offspring names
  Noffs_perpair
  if (Noffs_perpair>1){
    
    Offnames <- ceiling(Offs2[,1]/2)                
    Offnames2 <- paste0(Offnames,".",1:length(Offnames))
    Offnames2
    Offs <- cbind(paste("Offspring",Offnames2),Offs)
    Offs
    
  } else {
    
    Offnames <- ceiling(Offs2[,1]/2)
    Offs <- cbind(paste("Offspring",Offnames),Offs)
  } # end of if statement about offspring
  
  Moms <- cbind(paste("Mom",1:length(Moms[,1])),Moms)
  Moms
  Dads <- cbind(paste("Dad",1:length(Dads[,1])),Dads)
  Dads
  
  #add error======================================================================#
  error
  #first for the dads
  
  #getting the gts
  Dadsg <- Dads[,-1]
  
  #figuring out how many gts to make errors in
  ldad <- length(Dadsg)*error
  ldad
  #randomly selecting where to put the errors
  pdad1 <- sample(1:length(Dadsg),ldad,replace= FALSE)
  pdad1
  
  #randomly selecting some gets to replace in the those locations
  pdad2 <- Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
  pdad2
  
  #actually putting in the error
  Dadsg2 <- replace(Dadsg,pdad1,pdad2)
  Dadsg2    
  #replacing the old Dads gts
  Dads <- cbind(Dads[,1],Dadsg2)
  Dads
  
  #doing the same thing for the kids
  Offsg=Offs[,-1]
  loff=length(Offsg)*error
  poff1=sample(1:length(Offsg),loff,replace=FALSE)
  poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
  Offsg2=replace(Offsg,poff1,poff2)
  Offs=cbind(Offs[,1],Offsg2)
  
  #no error rate for moms
  #Momsg <- Moms[,-1]
  #ldad <- length(Momsg)*error
  #pdad1 <- sample(1:length(Momsg),ldad,replace=FALSE)
  #pdad2 <- Momsg[sample(1:length(Momsg),ldad,replace=FALSE)]
  #Momsg2 <- replace(Momsg,pdad1,pdad2)
  #Moms <- cbind(Moms[,1],Momsg2)
  
  #===============================================================================#
  #making generic locus names
  colsnmz <- rep("Locus",ncol(Offs)-1)
  colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
  colsnmz3 <- paste(colsnmz,colsnmz2)
  colsnmz3 <- c("IDs",colsnmz3)
  colsnmz3
  
  #usig them to make the colnames for the Offs, Dads, and Moms
  colnames(Offs)<- colsnmz3
  colnames(Dads)<- colsnmz3
  colnames(Moms)<- colsnmz3
  
  #Add unrealated individuals <- ====================================================#
  Nunrelated
  
  if (Nunrelated>0){
    
    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nunrelated*3                                                            
    afreqs <- afreq
    OUT <- NULL
    
    #making the allele frequencies again, this way they are completely random and unrelated
    #doing the same thing as lines 34-113
    sims <- function(sims) {
      alleles2 <- afreqs[which(afreqs[,1] == z),]
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])
      
      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2                                                         
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
      
      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))                                                   
      hetfreq <- 2*(hets[,1]*hets[,2])
      
      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
      hets2 <- cbind(hetvals,hetfreq)
      
      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)                                                      
      
      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000         
      
      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
      
      OUT <<-cbind(OUT,sg1)
      
    } #end of sims function
    
    #now applying function
    z1 <- length(unique(afreqs[,1]))
    for(z in 1:z1) {lapply(z,sims)}
    
    #repeating lines 122-204
    parents <- OUT
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    parents <- as.numeric(parents)-Anumbs
    
    #to called this thing unrelated
    unrelated <- cbind(paste("Individual",1:length(parents[,1])),parents)
    Lid <- length(unrelated[,1])/3
    
    #splitting the unrelated individuals among dads, moms, and kids
    m1 <- unrelated[1:Lid,]
    d1 <- unrelated[(Lid+1):(Lid*2),]
    j1 <- unrelated[((Lid*2)+1):(Lid*3),]
    
    #now for the column names again
    colsnmz <- rep("Locus",ncol(Offs)-1)
    colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
    colsnmz3 <- paste(colsnmz,colsnmz2)
    colsnmz3 <- c("IDs",colsnmz3)
    colnames(Offs)<- colsnmz3
    colnames(Dads)<- colsnmz3
    colnames(Moms)<- colsnmz3
    
    #adding the unrelated individuals in the three files
    Moms <- rbind(Moms,m1)
    Dads <- rbind(Dads,d1)
    Offs <- rbind(Offs,j1)
    
  } #end of unrelated if statement
  
  #removing the file SimOff.txt from working directory
  unlink(off.filename)
  
  #Begin add siblings=============================================================#
  #This scripts creates a desired number of siblings and splits them between adults and offspring files
  Nsibs <- as.numeric(Nsibs)
  if (Nsibs > 0) {
    #could change this, but as stands only evaluating 1 pair of siblings per parent-pair
    Noffs_per_pair <- 2                                                                
    
    #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
    Nadults <- Nsibs*2                                                                 
    afreqs <- afreq
    
    #repeating allele frequency simulation from above
    OUT <- NULL
    sims <- function(sims){
      alleles2 <- afreqs[which(afreqs[,1]==z),]
      
      #table allele frequencies
      alleles3 <- cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   
      
      #create homozygote allele frequencies
      homos <- (alleles3[,2])^2                                                          
      homos2 <- cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
      
      #create heterozygote allele frequencies
      hets <- t(combn(alleles3[,2],2))                                                   
      hetfreq <- 2*(hets[,1]*hets[,2])
      
      #create heterozygote allele names
      hetvals <- t(combn(as.character(alleles3[,1]),2))                                  
      hets2 <- cbind(hetvals,hetfreq)
      
      #combine hets and homos and create genotypes
      gfreqs <- rbind(hets2,homos2)                                                      
      
      #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
      n <- 1000000                                                                       
      
      #create genotypes(by coloumn, 1 for each allele)
      gfreqs1 <- rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            
      gfreqs2 <- rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
      gtypes <- cbind(gfreqs1,gfreqs2)
      gtypes <- gtypes[sample(1:length(gtypes[,1]),replace=F),]
      sg1 <- gtypes[sample(1:length(gtypes[,1]),Nadults),]
      
      OUT <<-cbind(OUT,sg1)
      
    } #end of allele frequency simulation
    
    #applyingn to making genotypes of parengs
    z <- length(unique(afreqs[,1]))
    for(z in 1:z) {lapply(z,sims)}
    parents <- OUT
    
    #repeating the part where he addes 1000 to each line
    c <- c(1:(ncol(OUT)))
    odd <- 2*(unique(round(((c-2))/2)))+1
    l <- length(odd) * 1000
    codes <- seq(from=1,to=l,by=1000)
    cols <- sort(rep(codes,2))-1
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)
    parents <- as.numeric(parents)+Anumbs
    
    
    #create full sib families (go down the list in pair)============================#
    OUT2 <- NULL
    sims <- function(sims)  {
      N <- 1:length(parents[,1])
      u <- sample(N,1)
      u2 <- sample(N,1)
      p1 <- parents[u,]
      p2 <- parents[u2,]
      als <- rep(1:2,length(p1)/2)
      
      #number of offspring per pair
      Noffs <- Noffs_per_pair                                                            
      OUT2 <- NULL
      for (b in 1:Noffs){
        
        #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
        pos1 <- sample(als,length(p1)/2,replace=TRUE)                                      
        pos2 <- sample(als,length(p1)/2,replace=TRUE)
        pos11 <- pos1+(seq(0,(length(p1)-1),2))
        pos22 <- pos2+(seq(0,(length(p2)-1),2))
        o1 <- p1[pos11]
        o2 <- p2[pos22]
        o3 <- sort(c(o1,o2))
        o3 <- t(c(z,o3))
        write.table(o3,file="SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
      } #end of for loop
    } # end of simulation for sibs
    
    #applying it to the parents
    z <- length(parents[,1])
    
    #used to move down list, now sample randomly so is just used to produce wanted number of offspring
    C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,sims)                 
    Dads <- parents[seq(from=1,to=length(parents[,1]),by=2),]
    Moms <- parents[seq(from=2,to=length(parents[,1]),by=2),]
    
    #see code before functions for adding 1000s (here am removing 1000s)
    Anumbs <- matrix(cols,Nadults,ncol(OUT),byrow=T)                                   
    parents <- as.numeric(parents)-Anumbs
    Dads <- parents[seq(from=1,to=length(parents[,1]),by=2),]
    Moms <- parents[seq(from=2,to=length(parents[,1]),by=2),]
    
    #reading in the kids again
    Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    Offs  <-  as.matrix(Offs2[,-1])
    Anumbs <- matrix(cols,length(Offs[,1]),ncol(OUT),byrow <- T)
    Offs <- as.numeric(Offs)-Anumbs
    
    #naming offspring
    if (Noffs_per_pair>1){
      #naming of offpsring
      Offnames <- ceiling(Offs2[,1]/2)                                                   
      Offnames2 <- paste(Offnames,".",1:length(Offnames))
      Offs <- cbind(paste("Sibling",Offnames2),Offs)
    } else {
      Offnames <- ceiling(Offs2[,1]/2)
      Offs <- cbind(paste("Sibling",Offnames),Offs)
    } # end of offspring naming if statement
    
    Moms <- cbind(paste("Mom",1:length(Moms[,1])),Moms)
    Dads <- cbind(paste("Dad",1:length(Dads[,1])),Dads)
    
    
    #add error======================================================================#
    Offsg <- Offs[,-1]
    loff <- length(Offsg)*error
    poff1 <- sample(1:length(Offsg),loff,replace=FALSE)
    poff2 <- Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
    Offsg2 <- replace(Offsg,poff1,poff2)
    Offs <- cbind(Offs[,1],Offsg2)
    
    #===============================================================================#
    colsnmz <- rep("Locus",ncol(Offs)-1)
    colsnmz2 <- sort(rep(1:(length(colsnmz)/2),2))
    colsnmz3 <- paste(colsnmz,colsnmz2)
    colsnmz3 <- c("IDs",colsnmz3)
    colnames(Offs)<- colsnmz3
    
    #calculate shared alleles among pairs of siblings===============================#
    sib1 <- Offs[seq(from=1,to=length(Offs[,1]),by=2),]
    sib2 <- Offs[seq(from=2,to=length(Offs[,1]),by=2),]
    write.table(sib1,file="Sib1.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(sib2,file="Sib2.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    
    #appending siblings to file 
    #getting dad siblings
    dadsibs <- read.table("Dads_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    sib1 <- read.table("Sib1.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    dadsibs <- rbind(dadsibs,sib1)
    #write.table(dadsibs,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)
    
    #putting them in the kid file
    juvsibs <- read.table("Juveniles_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    sib2 <- read.table("Sib2.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    juvsibs <- rbind(juvsibs,sib2)
    #write.table(juvsibs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=TRUE)
    
    #adding them to the end of the parent files
    Dads <- rbind(Dads,dadsibs)
    Offs <- rbind(Offs,juvsibs)
    
    #removing files
    unlink("Sib1.txt")
    unlink("Sib2.txt")
    unlink("SimOffs.txt")
  } # end of sibling function
  
  
  if(write2file == T){
    print("Mom, Dad, and Offspring genotypes have been saved to current working directory")  
    write.table(Dads,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(Moms,file="Moms_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
    write.table(Offs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
  } else {
    
    Parents <- rbind(Moms,Dads)
    Parents <- data.frame("Parents",Parents,stringsAsFactors=F)
    colnames(Parents) <- c("Type", colnames(Parents)[-1])
    
    Offspring <- data.frame("Offspring",Offs,stringsAsFactors=F)
    colnames(Offspring) <- c("Type",colnames(Offspring)[-1])
    
    output <- data.frame(rbind(Parents,Offspring),stringsAsFactors=F)
    #making the genotypes numeric
    output <- data.frame(output[,1:2],sapply(output[,3:ncol(output)], as.numeric),stringsAsFactors=F)
    return(output)
  }
} # end of sim.gt.data function

#########################
### no.par.bayes.ns() ###
#########################

#function the same as SOLOMON's no.par.bayes(), but with the slight change to the script to make 
#computation time faster


no.par.bayes.ns <- function(adults,offspring,reps=1000,Ngts=50000000)     {
  
  #defining *.txt name
  out.sims.filename <- paste("out.sims",Sys.time(), Sys.getpid(), sep = "_")
  Sort.txt.filename <- paste("Sort.txt",Sys.time(), Sys.getpid(), sep = "_")
  IdnamesA.txt.filename <- paste("IdnamesA.txt",Sys.time(), Sys.getpid(), sep = "_")
  IdnamesO.txt.filename <- paste("IdnamesO.txt",Sys.time(), Sys.getpid(), sep = "_")
  False_allele_freqs.txt.filename <- paste("False_allele_freqs.txt",Sys.time(), Sys.getpid(), sep = "_")
  True_allele_freqs.txt.filename <- paste("True_allele_freqs.txt",Sys.time(), Sys.getpid(), sep = "_")
  Shared_allele_freqs.txt.filename <- paste("Shared_allele_freqs.txt",Sys.time(), Sys.getpid(), sep = "_")
  lamdaphic.txt.filename <- paste("lamdaphic.txt",Sys.time(), Sys.getpid(), sep = "_")
  Putative.txt.filename <- paste("Putative.txt",Sys.time(), Sys.getpid(), sep = "_")
  True_shared_freqs.txt.filename <- paste("True_shared_freqs.txt",Sys.time(), Sys.getpid(), sep = "_")
  Output_genotypes.txt.filename <- paste("Output_genotypes.txt",Sys.time(), Sys.getpid(), sep = "_")
  #loading in adults and offspring genotypes
  Adults1 <- adults
  Offspring1 <- offspring
  
  #making sure the number of reps is numeric and counting the number of loci
  reps <- as.numeric(reps)
  loci <- ncol(Adults1)
  
  #removing ID column
  Adults <- Adults1[,c(2:loci)]                                                    
  Offspring <- Offspring1[,c(2:loci)]
  
  #getting the number loci, again for Progress bar
  total <- ncol(Adults)/2                                                       
  
  #putting in a progress bar
  
  
  #Begin Master simulation function
  afreqs <- function(afreqs){
    
    #making locus names? *** need to figure out ** looks like used in progress bar
    locus_name <- L
    
    #setting up progress bar again
    
    
    #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
    vect <- c(Adults[,L],Adults[,L+1])                                                 
    alleles <- data.frame(table(vect))
    alleles <- alleles[order(alleles[,1]),]
    if (as.numeric(as.character(alleles[1,1]))==0) {alleles <- alleles[-1,]}           #removes missing data
    if(length(alleles[,1])==1) {                                                    #deals with monomorphic loci by adding 1 very strange allele (799)
      alleles <- cbind(vect[1],alleles[2])
      alleles[2,]<-c(799,1)}
    alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
    homos <- (alleles2[,3])^2                                                          #create homozygote allele frequencies
    homos2 <- cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
    hets <- t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
    hetfreq <- 2*(hets[,1]*hets[,2])
    hetvals <- t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
    hets2 <- cbind(hetvals,hetfreq)
    gfreqs <- rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
    csum <- cumsum(as.numeric(gfreqs[,3]))
    gfreqs1 <- cbind(gfreqs,csum)
    Nadults <- length(Adults[,1])
    Noffs <- length(Offspring[,1])
    
    
    #===============================================================================#end locus-specific HWE genotype frequency calculations
    alength <- length(alleles2[,1])
    for (Y in 1:reps) {
      positions <- 1:length(gfreqs[,1])
      sg3 <- sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                      #sample the repeated genotype positions, by the number of adults
      sadults <- gfreqs[sg3,1:2]                                                         #index gfreqs to create genotypes
      og3 <- sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                        #create juvenile genotyes
      soffs <- gfreqs[og3,1:2]
      soffs <- cbind(as.numeric(soffs[,1]),as.numeric(soffs[,2]))
      asp <- cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
      asp <- rbind(cbind(locus_name,0,0),asp)
      for (i in 1:Nadults) {
        parent1 <- as.numeric(sadults[i,1])                                                #first allele in parent
        parent2 <- as.numeric(sadults[i,2])                                                #second allele in parent
        p1 <- soffs-parent1
        p2 <- soffs-parent2
        pp1 <- which(p1[,1]==0)
        pp2 <- which(p1[,2]==0)
        allele1 <- unique(c(pp1,pp2))
        p21 <- which(p2[,1]==0)
        p22 <- which(p2[,2]==0)
        allele2 <- unique(c(p21,p22))
        Out51 <- cbind(parent1,length(allele1))
        Out52 <- cbind(parent2,length(allele2))
        Out53 <- cbind(0,Noffs-length(unique(c(allele1,allele2))))
        Out5 <- rbind(Out51,Out52,Out53)
        if(parent2==parent1) {Out5 <- Out5[-1,]}                                           #remove 1 of alleles count if homozygous
        if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
          diffs <- sum(Out5[,2])-Noffs                                                     #comment out to be more conservative!
          maxa <- max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
          pos <- which(Out5[,2]==maxa)
          Out5[pos,2]<-Out5[pos,2]-diffs}
        m1 <- match(Out5[,1],asp[,2])
        m2 <- asp[m1,3]+as.numeric(Out5[,2])
        asp[m1,3]<-m2
        asp<-asp
      }
      write.table(asp,file=out.sims.filename,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  L <- ncol(Adults)
  C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
  
  #Bayes P-value==================================================================#
  OUT<- read.table(out.sims.filename, header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  locname <- unique(OUT[,1])                                                         #compile calculations for each locus
  OUTALL <- NULL
  for (z in locname) {
    Loc1 <- OUT[which(OUT[,1]==z),]
    allfreqs <- unique(Loc1[,2])
    OUT2 <- NULL
    for (x in allfreqs) {
      a1<-Loc1[which(Loc1[,2]==x),]
      a2 <- sum(a1[,3])
      a3 <- cbind(x,a2)
      OUT2 <- rbind(OUT2, a3)
    }
    OUT3 <- cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
    OUTALL <- rbind(OUTALL, cbind(z,OUT3))
  }
  #Create multilocus genotypes, 50000 at a time, calculate number of loci that mismatch, and calculate freqs of shared alleles
  NL <- length(unique(OUTALL[,1]))
  ngtypes <- 1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
  Ngts <- as.numeric(Ngts)
  inreps <- 50000     #was tested as 10000 for SNPS
  repnumber <- round(Ngts/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  asp <- cbind(0:NL,rep(0,length(0:NL)))
  
  for (n in 1:repnumber) {
    
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    distm <- apply(OUT, 1, function(x)sum(x == 1))
    distm2 <- data.frame(table(distm))
    m1 <- match(distm2[,1],asp[,1])
    m2 <- asp[m1,2]+distm2[,2]
    asp[m1,2]<-m2
    asp<-asp
  }
  
  #tabulate number of multilocus genotypes with x mismatches
  d2 <- asp
  d3 <- cbind(d2,d2[,2]/sum(d2[,2]))
  #===============================================================================# Create plot of exclusionary power.  Also, some necessary data formatting
  
  
  Adults<- adults
  Offs<- offspring
  Nads <- length(Adults[,1])
  Noffs <- length(Offs[,1])
  asp <- d3
  asp <- cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
  asp <- cbind(asp,cumsum(asp[,2]))
  asp <- cbind(asp,asp[,5]/max(asp[,5]))
  #find minimum Nloci to mismatch (could modify this) and perform exclusion with decided-upon mismatches==#
  distm <- cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
  #mismatch=min(which(round(distm[,6],1)==.9))                                    #deprecated
  mismatch <- which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
  Adults1 <- Adults                                                                   #begin exclusion
  Offspring1 <- Offs
  a <- ncol(Adults1)
  Adults <- Adults1[,c(2:a)]
  Offspring <- Offspring1[,c(2:a)]
  Anames <- Adults1[,1]
  Onames <- Offspring1[,1]
  Adults <- Adults1[,-1]
  Offspring <- Offspring1[,-1]
  categories <- ncol(Adults)
  Aindivids <- length(Adults[,1])
  Oindivids <- length(Offspring[,1])
  A <- 1:Aindivids
  O <- 1:Oindivids
  G <- expand.grid(A,O)
  AG <- G[,1]
  AO <- G[,2]
  Ads <- Adults[AG,]
  Offs <- Offspring[AO,]
  IdnamesA <- Anames[AG]
  IdnamesO <- Onames[AO]
  write.table(IdnamesA,file=IdnamesA.txt.filename,row.names=FALSE,col.names=F,sep="\t",append=F)
  write.table(IdnamesO,file=IdnamesO.txt.filename,row.names=FALSE,col.names=F,sep="\t",append=F)
  IdnamesA<- read.table(IdnamesA.txt.filename, header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  IdnamesO<- read.table(IdnamesO.txt.filename, header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  Names <- cbind(IdnamesA,IdnamesO)
  matches <- function(matches)
  {
    A <- Ads[,z]-Offs[,z]
    B <- Ads[,(z+1)]-Offs[,(z+1)]
    C <- Ads[,z]-Offs[,(z+1)]
    D <- Ads[,(z+1)]-Offs[,z]
    f <- A*B*C*D
    f <- (f^2)*10
    ss <- which(is.na(f)==TRUE)
    f <- replace(f,ss,0)
    identify <- which(f>0)
    f <- replace(f,identify,1)
    f <- cbind(z,f)
    write.table(f,file=Sort.txt.filename,row.names=FALSE,col.names=F,sep="\t",append=T)
  }
  z <- ncol(Ads)
  C1 <-  for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
  Observed<- read.table(Sort.txt.filename, header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  a <- length(U[,1])
  stuff <- rowSums(U)
  Sorted <- cbind(Names,stuff)
  matches <- which(Sorted[,3]<(mismatch+1))
  Actual <- sort(stuff)
  IDS <- which(stuff<(mismatch+1))
  PAdults <- Ads[IDS,]
  POffspring <- Offs[IDS,]
  nput <- length(matches)
  Putativepairs <- Sorted[matches,]
  PAdults <- cbind(Putativepairs[,c(1,3)],PAdults)
  POffspring <- cbind(Putativepairs[,c(2,3)],POffspring)
  names(PAdults)[1]<-"ID"
  names(POffspring)[1]<-"ID"
  head(PAdults)
  head(POffspring)
  PAdults$sorter <- as.numeric(row.names(PAdults))
  POffspring$sorter <- as.numeric(row.names(POffspring))
  Putative <- rbind(PAdults,POffspring)
  Putative <- Putative[order(Putative$sorter),]
  Putative$sorter <- NULL
  #head(Putative)
  #sorts <- function(sorts){
  #  tell <- rbind(PAdults[f,],POffspring[f,])
  #  write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
  #}
  #f <- length(PAdults[,1])
  #C1 <-  for(f in 1:f) lapply(f,sorts)
  
  unlink(Sort.txt.filename)
  unlink(IdnamesA.txt.filename)
  unlink(IdnamesO.txt.filename)
  #Calculate phi for each number mismatching loci=================================#
  
  #Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  
  observed <- data.frame(table(Putative[,2]))
  observed <- cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
  zerom <- 0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
  zerom2 <- which(is.na(match(zerom,observed[,1])))
  if (length(zerom2>0)) {observed <- observed[,3]
                         for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}  #not really 0, to prevent divide by 0
  }   else {observed=observed[,3]}
  expected <- distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
  #expected=distm[1:(mismatch+1),3]*Nads*Noffs                                     #not using cumulative sum
  phi <- expected/observed
  phi <- replace(phi,which(phi>=1),1)
  Offs<- offspring
  actualTrue <- length(grep("Off",Offs[,1]))
  info <- cbind(actualTrue,expected,observed,phi)
  #calculate phi and index values ================================================#
  
  phibase <- phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
  observed <- observed[min(which(phi[]==1))-1]                                       #do the same for observed
  testob <- which(observed==0.000000001)
  phi2 <- cbind(1:length(phi),phi)
  if (length(testob>0)) {phi4 <- phi2[-testob,]} else {phi4 <- phi2}
  nmismatch <- min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
  index <- phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
  index <- index[which(index>-1)]
  if (length(index)>1) {
    if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
  phi <- phi[index]
  index <- index-1
  #Create Plot ===================================================================#
  #pdf(file <- "Output_Dataset_Power.pdf")
  #x <- 0:(length(info[,1])-1)
  #y1 <- info[,3]
  #y2 <- info[,2]
  #p1 <- which(y2==0)
  #y2 <- replace(y2,p1,.000000001)
  #p1 <- which(y1<y2)
  #y3 <- replace(y1,p1,y2[p1])
  #y2 <- y2+1
  #y3 <- y3+1
  #par(mar = c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  #plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
  #ats <- c(0,1,2,3,4,5,6,7,8,9)
  #labs <- c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
  #axis(side=2,ats,labs)
  #lines(x,log10(y3),lwd=2)
  #lines(x,log10(y2),lwd=2)
  #points(x,log10(y3),pch=21,bg="green",cex=2)
  #points(x,log10(y2),pch=21,bg="blue",cex=2)
  #legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
  #yphi <- y2/y3
  #par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
  #plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
  #lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
  #points(x,yphi,cex=2,pch=21,bg="gray")
  #mtext("Number of Mismatching Loci",side=1,line=1.94)
  #dev.off()
  info2 <- cbind(x,info[,4])
  colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
  #write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  #===============================================================================#
  ngtypes <- 20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
  inreps <- 10000
  repnumber <- round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
  #writes values at all loci, to be analyzed further below
  for (n in 1:repnumber) {
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {
      Loco <- OUTALL[which(OUTALL[,1]==r),]
      alleles3 <- cbind(Loco,ngtypes*Loco[,4])
      findo <- which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
      findo2 <- replace(alleles3[,4],findo,1)
      alleles3 <- cbind(alleles3,findo2)
      alleles3 <- alleles3[-which(alleles3[,2]==0),]
      gtrue <- sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
    }
    for (i in index) {                                                              #loop over numbers of mismatched loci
      if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                    #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distp <- cbind(i,DIST)
      write.table(distp,file=False_allele_freqs.txt.filename,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  
  #Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
  OUT9 <- NULL
  for (n in index){
    Putative2 <- Putative[which(Putative[,2]==n),]
    write.table(Putative2, file=Putative.txt.filename,sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                                #need at least 1000 values , else assigned 0 (may be redundant now)
    OBS <- NULL
    afreqs <- function(afreqs) {
      PutL <- Putative2[,c(L+2,L+3)]
      PutLadults <- PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
      PutLoffs <- PutL[seq(from=2,to=length(PutL[,1]),by=2),]
      Puts <- cbind(PutLadults,PutLoffs)
      c1 <- Puts[,1]-Puts[,3]                                                    #find the matching alleles
      c2 <- Puts[,1]-Puts[,4]
      c3 <- Puts[,2]-Puts[,3]
      c4 <- Puts[,2]-Puts[,4]
      Puts2 <- cbind(Puts,c1,c2,c3,c4)
      P5 <- replace(Puts2[,5],which(Puts2[,5]!=0),-10)
      P6 <- replace(Puts2[,6],which(Puts2[,6]!=0),-10)
      P7 <- replace(Puts2[,7],which(Puts2[,7]!=0),-10)
      P8 <- replace(Puts2[,8],which(Puts2[,8]!=0),-10)
      P5 <- replace(P5,which(P5==0),1)
      P6 <- replace(P6,which(P6==0),1)
      P7 <- replace(P7,which(P7==0),1)
      P8 <- replace(P8,which(P8==0),1)
      Puts3 <- cbind(Puts,P5,P6,P7,P8)
      Puts4 <- cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
      alleles2 <- OUTALL[which(OUTALL[,1]==L),]
      alfreq1 <- alleles2[match(Puts4[,1],alleles2[,2]),4]
      alfreq2 <- alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
      alfreq3 <- alleles2[match(Puts4[,3],alleles2[,2]),4]
      alfreq4 <- alleles2[match(Puts4[,4],alleles2[,2]),4]
      Puts5 <- cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
      R1 <- replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
      R2 <- replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
      R3 <- replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
      R4 <- replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
      Puts6 <- cbind(R1,R2,R3,R4)
      Put_share <- apply(Puts6, 1, min)                                          #find row minimum
      Put_share2 <- apply(Puts4, 1, max)                                         #find shared allele name
      Put_share3 <- c(Put_share,Put_share2)
      OBS <<- cbind(OBS,Put_share3)
    }
    L <- ncol(Putative2)-2
    C1 <- for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths <- length(OBS[,1])/2
    if (lengths==1) {OBA3 <- t(OBS[(lengths+1):(2*lengths),])} else {OBA3 <- OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file=True_shared_freqs.txt.filename,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS <- t(OBS[1:lengths,])} else {OBS <- OBS[1:lengths,]}   #formatting for if there is only a single pair
    obsp <- apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
    OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
  }
  #calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#
  
  OBA3<- read.table(True_shared_freqs.txt.filename, header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
    OUT <- NULL
    for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
      vect <- c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
      alleles <- data.frame(table(vect))
      alleles <- alleles[order(alleles[,1]),]
      if (as.numeric(as.character(alleles[1,1]))<=0) {alleles <- alleles[-1,]}
      alleles2 <- cbind(alleles,alleles[,2]/sum(alleles[,2]))
      gtrue <- sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
      OUT <- cbind(OUT,gtrue)
      
      if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
        mm <- alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
        write.table(cbind(r,mm),file=Shared_allele_freqs.txt.filename,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
      }
      }
    }
    for (i in index) {
      if (i==0) {DIST <- as.data.frame(apply(OUT, 1, prod))} else {
        DIST <- NULL                                                                   #sample distribution by appropriate number of loci
        distp2 <- as.matrix(OUT)
        a1 <- NULL
        a2 <- NULL
        for (z in 1:length(distp2[,1])) {a1 <- rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
        for (p in 1:i) {a2 <- rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
        distp2[a2]<-1
        a3 <- apply(distp2, 1, prod)
        DIST<-as.data.frame(a3)
      }
      distt<-cbind(i,DIST)
      write.table(distt,file=True_allele_freqs.txt.filename,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
  }
  #Calculate lamdaphi=============================================================#
  
  Putative3<- read.table(Putative.txt.filename, header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Putadults <- Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
  Putoffs <- Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
  Names <- cbind(as.character(Putadults),as.character(Putoffs))
  empirical <- cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
  distp <- read.table(False_allele_freqs.txt.filename, header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  P2 <- NULL
  for (i in index) {                                                              #loop over numbers of mismatched loci
    empirical2 <- empirical[which(as.numeric(as.character(empirical[,3]))==i),]
    if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
    if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
      a3 <- distp[which(distp[,1]==i),2]
      DIST<-as.data.frame(a3)
      P <- NULL
      for (b in 1:length(empirical2[,1])) {
        p1 <- length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
        if (p1==0) {p1=0.00001}
        p2 <- cbind(empirical2[b,1],empirical2[b,2],p1)
        p3 <- cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
        P <- rbind(P,p3)
      }
    }
    P2<-rbind(P2,cbind(i,P))
  }
  lamdaphi <- as.numeric(as.character(P2[,5]))
  #Calculate lamda|phic ==========================================================#
  
  lamdaphic_dist<- read.table(True_allele_freqs.txt.filename, header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  Observed<- read.table(Shared_allele_freqs.txt.filename, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  mm1 <- which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
  mm1
  Observed[mm1,2] <- 1                                                            #replace NAs with 1
  a <- unique(Observed[,1])
  U <- NULL
  for (i in a) {
    u <- Observed[Observed[,1]==i,2]
    U <- cbind(U,u)
  }
  lamdaphic <- apply(U, 1, prod)
  l1 <- length(which(OUT9[,2]==0))
  if (l1>0) lamdaphic <- c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
  P3 <- cbind(P2,lamdaphic)
  for (i in index) {                                                              #loop over numbers of mismatched loci
    e2 <- P3[which(as.numeric(as.character(P3[,1]))==i),]
    if(length(e2)==6) {e2 <- t(e2)}                                                 #deals with one indiviudal formatting
    if (e2[1,5]==0) {write.table(e2[,6], file=lamdaphic.txt.filename ,sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
      a3 <- lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
      DIST<-as.data.frame(a3)
      for (b in 1:length(e2[,1])) {                                                 #calculate p values
        p1 <- length(which(DIST[,1] <= e2[b,6]))
        p2 <- p1/length(DIST[,1])
        write.table(p2, file=lamdaphic.txt.filename ,sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      }
    }
  }
  lamdaphic<- read.table(lamdaphic.txt.filename , header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
  #Put it all together with Bayes theorem!========================================#
  
  vals <- cbind(P2,lamdaphic[,1])
  philength <- cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
  phis <- rep(philength[,2],philength[,3])
  vals <- cbind(vals,phis)
  colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
  phi <- as.numeric(as.character(vals[,7]))
  lamdaphi <- as.numeric(as.character(vals[,5]))
  lamdaphic <- as.numeric(as.character(vals[,6]))
  lamdaphi <- replace(lamdaphi,which(lamdaphi==0),1)
  lamdaphic <- replace(lamdaphic,which(lamdaphic==0),1)
  pval <- (lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
  pval <- cbind(vals[,2],vals[,3],vals[,1],pval)
  pval <- pval[order(as.numeric(pval[,4])),]
  colnames(pval) <- c("Adult","Offspring","NL_mismatch","Probability of pair being false given frequencies of shared alleles")
  #write.table(vals, file="Posterior_Components_Bayes.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
  
  unlink(False_allele_freqs.txt.filename)                                                #clear all sims files
  unlink(True_allele_freqs.txt.filename)
  unlink(Shared_allele_freqs.txt.filename)
  unlink(lamdaphic.txt.filename )
  unlink(Putative.txt.filename)
  unlink(True_shared_freqs.txt.filename)
  unlink(Output_genotypes.txt.filename)
  unlink(out.sims.filename)
  rm(list=ls())
}

#fin!
