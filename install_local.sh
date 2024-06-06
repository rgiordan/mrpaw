#!/usr/bin/env bash

PKG_NAME="mrpaw"
R -e 'library(devtools); devtools::document("'$PKG_NAME'"); devtools::install_local("'$PKG_NAME'", force=TRUE, upgrade="never")'
