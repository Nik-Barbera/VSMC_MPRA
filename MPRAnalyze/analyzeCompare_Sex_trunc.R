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
dnaAnnot=read.table("/standard/vol192/civeleklab/Nik/MPRAnalyze/Data/5F_dnaAnnot_trunc.csv",sep=",",header=T)
rnaAnnot=read.table("/standard/vol192/civeleklab/Nik/MPRAnalyze/Data/5F_rnaAnnot_trunc.csv",sep=",",header=T)

rnaAnnot$condition <- as.factor(rnaAnnot$Condition)
rnaAnnot$barcode <- as.factor(rnaAnnot$barcode)
rnaAnnot$Sex <- as.factor(rnaAnnot$Sex)
#rnaAnnot$Allele <- as.factor(rnaAnnot$Allele)
dnaAnnot$condition <- as.factor(dnaAnnot$Condition)
dnaAnnot$barcode <- as.factor(dnaAnnot$barcode)
dnaAnnot$Sex <- as.factor(dnaAnnot$Sex)
#dnaAnnot$Allele <- as.factor(dnaAnnot$Allele)

print("Importing Count Files...")
dnaCounts=read.table("/standard/vol192/civeleklab/Nik/MPRAnalyze/Data/5F_dnaCounts_trunc.csv",sep=",",fill=TRUE,row.names=1,header=TRUE)
rnaCounts=read.table("/standard/vol192/civeleklab/Nik/MPRAnalyze/Data/5F_rnaCounts_trunc.csv",sep=",",fill=TRUE,row.names=1,header=TRUE)
dnaCounts[is.na(dnaCounts)]<-0
rnaCounts[is.na(rnaCounts)]<-0
# re-defining matrices as count matrices
dnaCounts<-data.matrix(dnaCounts)
storage.mode(dnaCounts)<-"integer"
rnaCounts<-data.matrix(rnaCounts)
storage.mode(rnaCounts)<-"integer"


##### Create MPRA Object #####
print("Creating MpraObject")

# create a logical array where every control sequence is labelled "TRUE", all else = "FALSE
controls = logical(length=length(dnaCounts[,1]))
controls[grep("NEG",row.names(dnaCounts))]=TRUE

##*****Added by JMH*****##
obj <- MpraObject(dnaCounts, 
                  rnaCounts, 
                  dnaAnnot = dnaAnnot, 
                  rnaAnnot = rnaAnnot, 
                  colAnnot = NULL, 
                  controls,
                  BPPARAM = bpparam)  
                 ## ^^^Just removed parentheses to test if object works

##### Library Size Normalization #####
print("Estimating Library Depths")
obj <- estimateDepthFactors(obj, lib.factor=c("replicate","Sex"),which.lib="both",depth.estimator = "totsum")

##### Analyze Sex Skew #####

# NOTE: This code is looking at whether a given candidate *sequence* shows sex-bias in enhancer activity
# it is *NOT* considering allelic skew.

print("analyzeComparative: Sex")
print("dnaDesign = ~ barcode + replicate + Sex")
print("rnaDesign = ~ Sex")
print("reducedDesign = ~ 1")

obj <- analyzeComparative(obj = obj, 
                          dnaDesign = ~ barcode + replicate + Sex, 
                          rnaDesign = ~ Sex, 
                          reducedDesign = ~ 1, 
                          BPPARAM=bpparam)

save.image(file = paste0("/standard/vol192/civeleklab/Nik/MPRAnalyze/Outputs/5F_analyzeCompare_Sex_trunc.RData"))
