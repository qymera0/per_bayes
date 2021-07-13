packageurl <- "http://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.19.3.tar.gz"

install.packages(packageurl, repos = NULL, type = "source")

renv::install(packageurl)