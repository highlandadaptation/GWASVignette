#read in the phenotype
#I bulk highlighted every phenotype I could find on the SEEDs website and pulled down the phenotype, the names in the phenotype don't correspond directly to the genotypes
#for me the most difficult part is connecting the germinate base name to the naming scheme used in the SNPs (which is in the column Sample.ID.of.DNA.from.single.plants.used.in.GWAS) once you get the right files though it's pretty much just a simple db join
fl<-read.table('DAPhenos.txt')

#now with the phenos read in, the next task is to get a single phenotype value that controls for the complicated design, multiple sites, landraces siring seed from elite dams, and 5-0+ replicate phenotypes per landrace genotype
tester<-fl$Tester.GID
plant<-fl$Sample.ID.of.DNA.from.single.plants.used.in.GWAS
test<-as.character(fl$dataset_description)
pheno<-fl$phenotypedata_value
#this just makes a prediction for each phenotype based upon the model
lmod<-lm(pheno ~ plant+test+tester)
pred<-predict(lmod)
#this finds the average prediction for each plant
plants<-tapply(pred,plant,mean)
#note we've substantially reduced the number of phenotypes here...
length(pheno);length(plants)
#get a weight for how many estimates went into the mean (it's unclear if this actually helps but it's possible that individuals with fewer estimates will have less power and more variability than those with more predictions)
weights<-tapply(pred,plant,length)

#read in the SNPs
#a quick example, this is a set of ~6,000 SNPs that are unimputed, relatively evenly polymorphic, and contain very little missing genotypes
#to make the SNPs in this manner I read the SeeDs vcfs into TASSEL and numericalized the sites. Thus, this is how R reads the format of the TASSEL output of numerical genotypes
SNPs<-read.table('HeavyFilterRefProb.txt',header=TRUE,skip=1,row.names = 1) 
dim(SNPs)
#I've only been interested in plants in Mexico/north CA because that is where sampling is the best so I need to remove individuals' SNPs that don't have phenotypes, and phenotypes that don't have SNPs
#split the SNP names and find the ones that are in the phenotypes
keeps<-sapply(names(plants),function(x) grep(paste(x,':',sep=""),rownames(SNPs)))
SNPs2<-SNPs[unlist(keeps),]
#next two lines are just memory management but when I run this on a million SNPs acros 2000 individuals it doesn't make sense to keep around duplicated SNP tables
rm(SNPs)
gc()
#now split the SNP names to be the same as the phenotype names
nms<-sapply(rownames(SNPs2),function(x) strsplit(x,':')[[1]][1])
#there are a few weird duplicated geno names, remove them and rename the SNP table
SNPs2<-SNPs2[!duplicated(nms),]
nms<-nms[!duplicated(nms)]
rownames(SNPs2)<-nms
#now remove anything with a phenotype that we don't have a genotype for
plants<-plants[names(plants) %in% rownames(SNPs2)]
weights<-weights[names(weights) %in% rownames(SNPs2)]
#this gets everything all orderly
SNPs2<-SNPs2[names(plants),]
dim(SNPs2);length(plants);length(weights)
#I do the SNP by SNPs in parallel, if you don't want to you can skip the next three lines and change the parApply command on line 51 and 54 a simple apply command
library(parallel)
cl <- makeCluster(3)
clusterExport(cl=cl,list('plants','weights')) #this passes extra objects to the cluster that you make

#if you I removed some individuals from the analyses I made some loci monomorphic, do the next two lines, otherwise skip
#slopes<-parApply(cl=cl,SNPs2,MARGIN=2, function(x) lm(plants ~ x )$coefficients[2]) #here I'm just pulling out the p-values
#SNPs2<-SNPs2[,-which(is.na(slopes))]
#Now for the real analysis the more you add to the SNP-by-SNP step (e.g. population structure) the slower the whole thing becomes
#also the way that I'm doing it with lm doesn't account for random variables (which is why people tend to prefer things like EMMA
vals<-parApply(cl=cl,SNPs2,MARGIN=2, function(x) -log10(summary(lm(plants ~ x ,weights=weights))$coefficients[2,4]))
stopCluster(cl=cl)

#make the manhattan table
library(qqman)
bp<-as.numeric(sapply(colnames(SNPs2),function(x) strsplit(x,'_')[[1]][2]))
chr<-as.numeric(sapply(sapply(colnames(SNPs2),function(x) strsplit(x,'_')[[1]][1]),substr,2,4))
qqp<-data.frame(SNP=paste('s',1:length(chr),sep=""),CHR=chr,BP=bp,P=vals)
manhattan(qqp,logp = FALSE)

