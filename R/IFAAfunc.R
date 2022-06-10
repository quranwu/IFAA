##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' Make inference on the association of microbiome with covariates
##'
##' Most of the time, users just need to feed the first five inputs to the function: `MicrobData`, `CovData`, `linkIDname`, `testCov` and `ctrlCov`. All other inputs can just take their default values.
##' To model the association, the following equation is used:
##'
##' \loadmathjax
##'
#' \mjdeqn{\log(\mathcal{Y}_i^k)|\mathcal{Y}_i^k>0=\beta^{0k}+X_i^T\beta^k+W_i^T\gamma^k+Z_i^Tb_i+\epsilon_i^k,\hspace{0.2cm}k=1,...,K+1}{}
##' where
##' - \mjeqn{\mathcal{Y}_i^k}{} is the AA of taxa \mjeqn{k}{} in subject \mjeqn{i}{} in the entire ecosystem.
##' - \mjeqn{X_i}{} is the covariate matrix.
##' - \mjeqn{W_i}{} is the confounder matrix.
##' - \mjeqn{Z_i}{} is the design matrix for random effects.
##' - \mjeqn{\beta^k}{} is the regression coefficients that will be estimated and tested with the `IFAA()` function.
##'
##' The challenge in microbiome analysis is that \mjeqn{\mathcal{Y}_i^k}{} can not be observed. What is observed is its small proportion:
##' \mjeqn{Y_i^k=C_i\mathcal{Y}^k_i}{}, where \mjeqn{C_i}{} is an unknown number between 0
##' and 1 that denote the observed proportion.
##'
##' The IFAA method can successfully addressed this challenge. The `IFAA()` will estimate the parameter
##' \mjeqn{\beta^k}{} and their 95% confidence intervals. High-dimensional \mjeqn{X_i}{} is handled by
##' regularization.
##'
##' @param experiment_dat A SummarizedExperiment object containing countData and colData. The countData contains microbiome absolute abundance or relative abundance with each column per
##' sample and each row per taxon/OTU/ASV (or any other unit). The colData contains covariates and confounders with each row
##' per sample and each column per variable. Note that the variables in colData has to be numeric or binary.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `FALSE`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nRefMaxForEsti The maximum number of final reference taxa used in phase 2. The default is `2`.
##' @param refTaxa A vector of taxa or OTU or ASV names. These are reference taxa specified by the user to be used in phase 1. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly.
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `testCov`. Default is `0.15`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. The default is `500`.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, the design matrix for X will be standardized in the analyses and the results. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `FALSE`. This argument could be useful for debug.
##' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
##' @param SDThresh The threshold of standard deviations of sequencing reads for been chosen as the reference taxon in phase 2. The default is `0.05` which means the standard deviation of sequencing reads should be at least `0.05` in order to be chosen as reference taxon.
##' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as reference taxon. The default is `0`.
##' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default number is `0.2` which means at least 20% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
##' @param seed Random seed for reproducibility. Default is `1`. It can be set to be NULL to remove seeding.
##' @return An SummarizedExperiment object with assay and metadata.
##'
##' - `assay`: The assay is a dataset for mean results from two reference taxon with each row representing each taxon, columns as "taxon", "cov", "estimate",
##' "SE.est", "CI.low", "CI.up", "adj.p.value", and "sig_ind", describing the taxon name, covariate name, parameter estimates, standard error estimates, lower bound and upper bound of the 95% confidence interval, adjusted p value, and the indicator showing whether the effect of corresponding taxon is significant, respectively.
##'
##' - `metadata`: The metadata is a list. `full_result_sep` contains estimates for each reference taxon separately, `final_ref_taxon` shows the final 2 reference taxon used for analysis, `ref_taxon_count` and `ref_taxon_est` shows how the final reference taxon were selected, and the other elements show the total time used in minutes, random seed used, the p-value cutoff used, and p-value adjustment method used.
##'
##' @examples
##' library(SummarizedExperiment)
##' data(dataM)
##' dim(dataM)
##' dataM[1:5, 1:8]
##' data(dataC)
##' dim(dataC)
##' dataC[1:5, ]
##' \donttest{
##' test_dat<-SummarizedExperiment(assays=list(counts=dataM), colData=dataC)
##' results <- IFAA(experiment_dat = test_dat,
##'                 testCov = c("v1", "v2"),
##'                 ctrlCov = c("v3"),
##'                 fdrRate = 0.15)
##'
##' ## to extract full results:
##' summary_res<-assay(results)
##' ## to extract significant results:
##' summary_res[summary_res$sig_ind==1,,drop=FALSE]
##' ## to extract metadata
##' meta_dat<-metadata(results)
##' ## to extract separate estimates
##' meta_dat$full_result_sep
##'}
##'
##'
##'
##' @references Li et al.(2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Liu et al.(2020) A bootstrap lasso + partial ridge method to construct confidence intervals for parameters in high-dimensional sparse linear models. Statistica Sinica
##' @importFrom foreach `%dopar%` foreach
##' @importFrom methods as
##' @importFrom future availableCores
##' @importFrom Matrix Diagonal Matrix
##' @importFrom HDCI bootLOPR
##' @importFrom qlcMatrix corSparse
##' @import expm
##' @import mathjaxr
##' @import glmnet
##' @import stats
##' @import utils
##' @importFrom SummarizedExperiment assays colData SummarizedExperiment
##' @importFrom gtools mixedorder
##' @export
##' @md

## arg deleted: nperm, x1permt, reguMethod, method, bootLassoAlpha


IFAA=function(
  experiment_dat,
  testCov=NULL,
  ctrlCov=NULL,
  testMany=TRUE,
  ctrlMany=FALSE,
  nRef=40,
  nRefMaxForEsti=2,
  refTaxa=NULL,
  adjust_method="BY",
  fdrRate=0.15,
  paraJobs=NULL,
  bootB=500,
  standardize=FALSE,
  sequentialRun=FALSE,
  refReadsThresh=0.2,
  taxDropThresh=0,
  SDThresh=0.05,
  SDquantilThresh=0,
  balanceCut=0.2,
  seed=1
){
  allFunc=allUserFunc()

  results=list()
  start.time = proc.time()[3]
  MicrobData<-data.frame(t(assays(experiment_dat)$counts))
  MicrobData$id<-rownames(MicrobData)
  CovData<-data.frame(colData(experiment_dat))
  CovData$id<-rownames(CovData)
  linkIDname<-"id"

  runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                   linkIDname=linkIDname,testCov=testCov,
                   ctrlCov=ctrlCov,testMany=testMany,
                   ctrlMany=ctrlMany,taxDropThresh=taxDropThresh,
                   standardize=standardize)
  data=runMeta$data
  results$covariatesData=runMeta$covariatesData
  binaryInd=runMeta$binaryInd
  covsPrefix=runMeta$covsPrefix
  Mprefix=runMeta$Mprefix
  testCovInd=runMeta$testCovInd
  testCovInOrder=runMeta$testCovInOrder
  testCovInNewNam=runMeta$testCovInNewNam
  ctrlCov=runMeta$ctrlCov
  microbName=runMeta$microbName
  newMicrobNames=runMeta$newMicrobNames
  results$covriateNames=runMeta$xNames

  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]


  rm(runMeta)

  if(length(refTaxa)>0){
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa in phase 1 have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels.")
    }
  }


  if (nRefMaxForEsti<2) {
    nRefMaxForEsti<-2
    warning("Warning: Needs at least two final reference taxon for estimation.")
  }


  if(nRef>(length(microbName))){
    stop("Error: number of random reference taxa can not be larger than the total number
           of taxa in the data. Try lower nRef")
  }

  refTaxa_newNam=newMicrobNames[microbName%in%refTaxa]


  results$analysisResults=Regulariz(data=data,testCovInd=testCovInd,
                                    testCovInOrder=testCovInOrder,
                                    testCovInNewNam=testCovInNewNam,
                                    microbName=microbName,nRef=nRef,
                                    nRefMaxForEsti=nRefMaxForEsti,
                                    binaryInd=binaryInd,
                                    binaryInd_test=binaryInd_test,
                                    covsPrefix=covsPrefix,Mprefix=Mprefix,
                                    refTaxa=refTaxa_newNam,
                                    paraJobs=paraJobs,
                                    adjust_method=adjust_method,
                                    fwerRate=fdrRate,
                                    bootB=bootB,
                                    sequentialRun=sequentialRun,
                                    allFunc=allFunc,refReadsThresh=refReadsThresh,
                                    SDThresh=SDThresh,
                                    SDquantilThresh=SDquantilThresh,
                                    balanceCut=balanceCut,seed=seed
  )
  rm(data)





  totalTimeMins = (proc.time()[3] - start.time)/60
  message("The entire analysis took ",round(totalTimeMins,2), " minutes")

  if(length(seed)==1){
    results$seed=seed
  }else{
    results$seed="No seed used."
  }

  output_se_obj<-SummarizedExperiment(assays  = list(results$analysisResults$full_results),
                                      metadata = list(full_result_sep=results$analysisResults$all_cov_list_sep,
                                                      final_ref_taxon=results$analysisResults$finalizedBootRefTaxon,
                                                      ref_taxon_count=results$analysisResults$goodIndpRefTaxWithCount,
                                                      ref_taxon_est=results$analysisResults$goodIndpRefTaxWithEst,
                                                      totalTimeMins=totalTimeMins,
                                                      seed=seed,
                                                      fdrRate=fdrRate,
                                                      adjust_method=adjust_method))



  return(output_se_obj)
}


