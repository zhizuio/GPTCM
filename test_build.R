library(styler)
style_pkg(pkg = "/Users/zhiz/zz/github/GPTCM/R")
library(lintr)
lint_dir(path = "/Users/zhiz/zz/github/GPTCM/R")
style_dir(path="/Users/zhiz/zz/github/GPTCM/vignettes", filetype=c("R", "Rmd"))
lintr::lint("/Users/zhiz/zz/github/GPTCM/vignettes/GPTCM.Rmd")

## Build a new version of the package
remove.packages("GPTCM")
rm(list = ls())
Rcpp::compileAttributes(pkgdir="/Users/zhiz/zz/github/GPTCM")
tools::resaveRdaFiles("/Users/zhiz/zz/github/GPTCM/data/simData.rda", compress="xz")
devtools::document("/Users/zhiz/zz/github/GPTCM")
#devtools::build("/Users/zhiz/zz/github/GPTCM",vignettes=TRUE)
#devtools::build("/Users/zhiz/zz/github/GPTCM")
#tools::package_native_routine_registration_skeleton('/Users/zhiz/zz/github/GPTCM', character_only = FALSE)
devtools::build("/Users/zhiz/zz/github/GPTCM", vignettes=T, args="--compact-vignettes=both")

# 1x-193:GPTCM zhiz$ R CMD build /Users/zhiz/zz/github/GPTCM --compact-vignettes=both
# 1x-193:GPTCM zhiz$ R CMD check /Users/zhiz/zz/github/GPTCM_0.4.tar.gz --as-cran --run-donttest
# devtools::check("/Users/zhiz/zz/github/GPTCM", cran=TRUE, vignettes = TRUE, build_args="--compact-vignettes=both")
# _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_
# library(rhub)
# rhub::check("/Users/zhiz/zz/github/GPTCM_1.1-1.tar.gz")
# pkg_path <- "/Users/zhiz/zz/github/GPTCM"
# chk <- local_check_solaris(pkg_path)

## Install the package
#install.packages("GPTCM_0.7.tar.gz",repos = NULL,type = "source")
install.packages("/Users/zhiz/zz/github/GPTCM_0.0.1.tar.gz",repos = NULL,type = "source")
unloadNamespace("GPTCM")
