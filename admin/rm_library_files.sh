#!/bin/bash

#-------------------------------#
# SET PATH TO DIRECTORY
WORKDIR=$HOME/source_code/pCMF
#-------------------------------#

## delete *.o files if exists (generated by sourceCpp)
find $WORKDIR -type f -name "*.o" | xargs -I {} rm {}

## delete *.o files if exists (generated by sourceCpp)
find $WORKDIR -type f -name "*.so" | xargs -I {} rm {}
