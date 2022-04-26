## Test environments
* local  x86_64-apple-darwin17.0 (64-bit), release R 4.0.2 (2020-06-22), 
* Windows Server 2022, R-devel, 64 bit
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

## Rhub check results

  `Minor bug fixes`
  `Package has help file(s) containing install/render-stage \Sexpr{} expressions but no prebuilt PDF manual.`
  
  This is caused by the mathjaxr package which shows Latex equations in IFAA() and MZILN() documentation.
  
  `'qpdf' is needed for checks on size reduction of PDFs`
  
  It seems to be something wrong with Rhub windows 2022 server, since qpdf should be installed on that end. 
  
  `Maintainer: 'Zhigang Li <zhigang.li@ufl.edu>'`
  
## Downstream dependencies
There are currently no downstream dependencies for this package
