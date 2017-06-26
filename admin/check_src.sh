#!/bin/bash

#-------------------------------#
# SET PATH TO DIRECTORY
WORKDIR=$(git rev-parse --show-toplevel)
#-------------------------------#

## current version
VERSION=$(cat $WORKDIR/Version)

## check
R CMD check $WORKDIR/sources/pCMF_${VERSION}.tar.gz
