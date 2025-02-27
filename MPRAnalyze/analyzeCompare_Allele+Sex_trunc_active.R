library(MPRAnalyze)
library(BiocParallel)
library(batchtools)

##### Capture the number of cores and set #####
numCores <- as.numeric(commandArgs(trailingOnly=TRUE)) - 1
print(numCores)
options(mc.cores=numCores)
#register(MulticoreParam(numCores)) #Note:  register used detectCores()
bpparam <- MulticoreParam(numCores)
print(multicoreWorkers())
show(bpparam)

##### IMPORT COUNT AND ANNOTATION FILES #####
print("Importing Annotation Files...")
dnaAnnot=read.table("/standard/vol192/civeleklab/Nik/MPRAnalyze/Data/5F_ByAllele_dnaAnnot_trunc.csv",sep=",",header=T)
rnaAnnot=read.table("/standard/vol192/civeleklab/Nik/MPRAnalyze/Data/5F_ByAllele_rnaAnnot_trunc.csv",sep=",",header=T)

# Make sure annotation files are formatted correctly
rnaAnnot$condition <- as.factor(rnaAnnot$Condition)
rnaAnnot$barcode <- as.factor(rnaAnnot$barcode)
rnaAnnot$Sex <- as.factor(rnaAnnot$Sex)
rnaAnnot$Allele <- as.factor(rnaAnnot$Allele)
dnaAnnot$condition <- as.factor(dnaAnnot$Condition)
dnaAnnot$barcode <- as.factor(dnaAnnot$barcode)
dnaAnnot$Sex <- as.factor(dnaAnnot$Sex)
dnaAnnot$Allele <- as.factor(dnaAnnot$Allele)

# Import DNA/RNA Count Tables
print("Importing Count Files...")
dnaCounts=read.table("{/Your/Directory/Here}/5F_ByAllele_dnaCounts_trunc_sex_active.csv",sep=",",fill=TRUE,row.names=1,header=TRUE)
rnaCounts=read.table("{/Your/Directory/Here}/5F_ByAllele_rnaCounts_trunc_sex_active.csv",sep=",",fill=TRUE,row.names=1,header=TRUE)

# Re-define any NA's as 0's
dnaCounts[is.na(dnaCounts)]<-0
rnaCounts[is.na(rnaCounts)]<-0

# Re-define matrices as count matrices
dnaCounts<-data.matrix(dnaCounts)
storage.mode(dnaCounts)<-"integer"
rnaCounts<-data.matrix(rnaCounts)
storage.mode(rnaCounts)<-"integer"


##### Create MPRA Object #####
print("Creating MpraObject")

# create a logical array where every control sequence is labelled "TRUE", all else = "FALSE
controls = logical(length=length(dnaCounts[,1]))
controls[grep("NEG",row.names(dnaCounts))]=TRUE

obj <- MpraObject(dnaCounts, 
                  rnaCounts, 
                  dnaAnnot = dnaAnnot, 
                  rnaAnnot = rnaAnnot, 
                  colAnnot = NULL, 
                  controls,
                  BPPARAM = bpparam)  

##### Library Size Normalization #####
print("Estimating Library Depths")
obj <- estimateDepthFactors(obj, lib.factor=c("replicate","Allele","Sex"),which.lib="both",depth.estimator = "totsum")

### Quantification Analysis
# Quantification analysis is addressing the question of what is the transription 
# rate for each enhancer in the dataset. These estimates can then be used to 
# identify and classify active enhancers that induce a higher transcription rate. 
# Regarding model design - this data is from a paired experiment, so DNA factors
# are fully transferable to the RNA model. For the RNA, we will be interested in 
# having a separate estimate of transcription rate for each condition (chromosomal
# and episomal), so this is the only factor included in the RNA model. Finally,
# fitting the model is done by calling the analyzeQuantification function:

##### Analyze Allele-by-Sez #####

# NOTE: This code is looking at whether a given SNP shows sex-bias in its allelic effect.

print("analyzeComparative: Sex")
print("dnaDesign = ~ barcode + replicate + Allele + Sex")
print("rnaDesign = ~ Sex*Allele")
print("reducedDesign = ~ Allele")

obj <- analyzeComparative(obj = obj, 
                          dnaDesign = ~ barcode + replicate + Allele + Sex, 
                          rnaDesign = ~ Sex*Allele, 
                          reducedDesign = ~ Allele, 
                          BPPARAM=bpparam)

save.image(file = paste0("{/Your/Directory/Here}/5F_analyzeCompare_Allele+Sex_trunc_active.RData"))
