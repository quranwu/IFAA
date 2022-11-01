set.seed(1)

## create an ID variable for the example data
ID=seq_len(5)

## generate three covariates x1, x2, and x3, with x2 binary
x1<-rnorm(5)
x2<-rbinom(5,1,0.5)
x3<-rnorm(5)
dataC<-data.frame(cbind(ID,x1,x2,x3))

## Coefficients for x1, x2, and x3 among 20 taxa.
beta_1<-c(1,rep(0,19))
beta_2<-c(0,2,rep(0,18))
beta_3<-rnorm(20)
beta_mat<-cbind(beta_1,beta_2,beta_3)

## Generate absolute abundance for 20 taxa in ecosystem.
dataM_eco<-floor(exp(10+as.matrix(dataC[,-1]) %*% t(beta_mat) + rnorm(100,sd=0.05)))

## Generate sequence depth and generate observed abundance
Ci<-runif(5,0.01,0.05)
dataM<-floor(apply(dataM_eco,2,function(x) x*Ci))
colnames(dataM)<-paste0("rawCount",1:20)

## Randomly introduce 0 to make 25% sparsity level.
# dataM[sample(seq_len(length(dataM)),length(dataM)/4)]<-0

dataM<-data.frame(cbind(ID,dataM))

dataC_add<-data.frame(x4=rnorm(5),x5=rnorm(5),x6=rnorm(5))
dataC<-cbind(dataC,dataC_add)

## The following steps are to create a SummarizedExperiment object.
## If you already have a SummarizedExperiment format data, you can
## ignore the following steps and directly feed it to the IFAA function.

## Merge two dataset by ID variable
data_merged<-merge(dataM,dataC,by="ID",all=FALSE)

## Seperate microbiome data and covariate data, drop ID variable from microbiome data
dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("ID")]]
dataC_sub<-data_merged[,colnames(dataC)]

## Create SummarizedExperiment object
test_dat<-SummarizedExperiment::SummarizedExperiment(
  assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)

## Create a SummarizedExperiment object
test_dat<-SummarizedExperiment::SummarizedExperiment(assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)


test_IFAA_nz <- function() {
  results <- IFAA(experiment_dat = test_dat,
                  testCov = c("x1","x2"),
                  ctrlCov = c("x3","x4","x5","x6"),
                  fdrRate = 0.1,
                  nRef = 4,
                  paraJobs = 2)
  summary_res<-results$full_result
  sig_results<-subset(summary_res,sig_ind==TRUE)
  checkEquals(sig_results$taxon, c("rawCount1","rawCount2"))
  
}