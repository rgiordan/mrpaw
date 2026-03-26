#!/usr/bin/env bash

PKG_NAME="mrpaw"
R -e '
library(devtools); 
devtools::document("."); 
devtools::install_local(".", force=TRUE, upgrade="never", dependencies=TRUE)
'
