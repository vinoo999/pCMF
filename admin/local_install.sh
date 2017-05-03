#!/bin/bash

#-------------------------------#
# SET PATH TO DIRECTORY
WORKDIR="$HOME/source_code/cmf"
#-------------------------------#

## install directory
INSTALLDIR="$WORKDIR/installDir"

if [ ! -d $INSTALLDIR ]; then mkdir -p $INSTALLDIR; fi

## source
VERSION=`cat $WORKDIR/Version`
SRCDIR="$WORKDIR/sources"
SRC="$SRCDIR/cmf_${VERSION}.tar.gz"

## INSTALLATION
R CMD INSTALL --library=$INSTALLDIR $SRC
