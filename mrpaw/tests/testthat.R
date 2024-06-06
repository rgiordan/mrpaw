#!/usr/bin/env Rscript

library(devtools)
#devtools::load_all()
library(testthat)
library(mrpaw)

test_check("mrpaw")
#test_file("testthat/test_mrpaw.R")
