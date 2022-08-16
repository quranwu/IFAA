set.seed(1)
## generate three covariates x1, x2, and x3, with x2 binary
x1<-rnorm(40)
x2<-rbinom(40,1,0.5)
x3<-rnorm(40)
dataC<-cbind(x1,x2,x3)

## Coefficients for x1, x2, and x3 among 60 taxa.
beta_1<-c(0.1,rep(0,59))
beta_2<-c(0,0.2,rep(0,58))
beta_3<-rnorm(60)
beta_mat<-cbind(beta_1,beta_2,beta_3)

## Generate absolute abundance for 60 taxa in ecosystem.
dataM_eco<-floor(exp(10+dataC %*% t(beta_mat) + rnorm(2400,sd=0.05)))

## Generate sequence depth and generate observed abundance
Ci<-runif(40,0.01,0.05)
dataM<-floor(apply(dataM_eco,2,function(x) x*Ci))
colnames(dataM)<-paste0("rawCount",1:60)

## Randomly introduce 0 to make 25% sparsity level.
dataM[sample(seq_len(length(dataM)),length(dataM)/4)]<-0

test_dat<-SummarizedExperiment::SummarizedExperiment(
  assays=list(MicrobData=t(dataM)), colData=dataC)


test_IFAA <- function() {
  results <- IFAA(experiment_dat = test_dat,
                  testCov = c("x1","x2"),
                  ctrlCov = c("x3"),
                  fdrRate = 0.1,
                  nRef = 20,
                  paraJobs = 2)
  summary_res<-results$full_result
  sig_results<-subset(summary_res,sig_ind==TRUE)
  checkEquals(sig_results$taxon, c("rawCount1","rawCount2"))

}