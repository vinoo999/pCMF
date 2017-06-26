#!/bin/bash

#-------------------------------#
# SET PATH TO DIRECTORY
WORKDIR=$(git rev-parse --show-toplevel)
#-------------------------------#

## install directory
INSTALLDIR="$WORKDIR/installDir"

if [ ! -d $INSTALLDIR ]; then mkdir -p $INSTALLDIR; fi

## source
VERSION=$(cat $WORKDIR/Version)
SRCDIR="$WORKDIR/sources"
SRC="$SRCDIR/pCMF_${VERSION}.tar.gz"

## remove former version
rm -rf $INSTALLDIR/pCMF

## INSTALLATION
R CMD INSTALL --library=$INSTALLDIR $SRC
