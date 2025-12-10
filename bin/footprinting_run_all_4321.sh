
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Footrpints Detector v1.0
# Copyright (c) 2008 by Xiaoyu Chen and William Noble
# All rights reserved.
# Redistribution is not permitted without the permission of the authors.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#!/bin/bash

set -e

if [ ! $# = 6 ]
then
    echo "USAGE:"
    echo "footprinting_run_all.sh  yeast.dnasI.tagCounts.bed  yeast.unmappableBase.bed  yeast.intergenic.bed  yeast.footprints.bed tag scripts_dir"
    exit 1
fi


cutCountsFile=$1
mappableFile=$2
DHSBedFile=$3
ftprtsFile=$4


## ****** change parameters as needed ******

# the directory keeping generated files
resDir=./

# the size of flanking region considered for each intergenic region
# this number should >= (localBgrWidth-1) / 2 
# (localBgrWidth=151 see below) 
flankSize=75

# tag for files generated
tag=$5
scripts_dir=$6
## ******


cutCountsMtx=$tag.cutCounts.mtx
mappableMtx=$tag.mappable.mtx
finalDHSBed=$tag.flanking.bed



cd $resDir

if [ ! -f $finalDHSBed ]
then
  awk -v fs=$flankSize '{print $1 "\t" $2-fs "\t" $3+fs}' $DHSBedFile > $finalDHSBed
fi

echo "extracting tag counts ..."

# Compute a matrix. Entries are cutpoint counts 
# each row is a DHS, columns are postions
if [ ! -f $cutCountsMtx ]
then
  footprinting_extract_signal.py \
      -no-region-ids -no-region-scores \
      -varied-width -none-label NaN \
      $finalDHSBed $cutCountsFile \
      > $cutCountsMtx
fi


echo "extracting mappability ..."

# Compute a matrix. Entries are binary.
# each row is a DHS, columns are positions 
# "1" means that the position are mappable.
if [ ! -f $mappableMtx ]
then
  footprinting_extract_signal.py \
      -no-region-ids -no-region-scores \
      -varied-width -none-label 1 \
      $finalDHSBed $mappableFile \
      > $mappableMtx
fi

# Check if either cutCountsMtx or finalDHSBed has zero lines
echo "checking for empty input files..."
cutCountsLines=$(wc -l < "$cutCountsMtx" 2>/dev/null || echo "0")
finalDHSLines=$(wc -l < "$finalDHSBed" 2>/dev/null || echo "0")

if [ "$cutCountsLines" -eq 0 ] || [ "$finalDHSLines" -eq 0 ]; then
    echo "Warning: Either $cutCountsMtx (${cutCountsLines} lines) or $finalDHSBed (${finalDHSLines} lines) has zero lines."
    
    # Create empty final output file
    touch "$resDir/$ftprtsFile"
    
    # Clean up intermediate files that may exist
    [ -f "$resDir/$finalDHSBed" ] && rm "$resDir/$finalDHSBed"
    [ -f "$resDir/$cutCountsMtx" ] && rm "$resDir/$cutCountsMtx"
    [ -f "$resDir/$mappableMtx" ] && rm "$resDir/$mappableMtx"
    [ -f "$resDir/$outPrefix.datatemp" ] && rm "$resDir/$outPrefix.datatemp"
    [ -f "$resDir/$outPrefix.nulltemp" ] && rm "$resDir/$outPrefix.nulltemp"
    
    echo "Empty output created at $resDir/$ftprtsFile"
    exit 0
fi



# max size of an intergenic region
# for PC adjustment
# maxColNum=15000
# for non-PC (qnorm by Li Dubai)
maxColNum=30000
# q-value threshold for output footprints
qValueThresh=1
# min size of footprints 
minK=8
# max size of footprints
maxK=30
# prefix for the file of generated footprints
outPrefix=$tag.footprints
# size of the local background window
localBgrWidth=151

# some other parameters
transformType=1
isDHS=0
KStep=1
isOneSideBkg=0



echo "detecting footprints..."
set -x
octave --eval "pkg load statistics; addpath('$scripts_dir'); footprinting_rank('$cutCountsMtx', '$mappableMtx', $minK, $maxK, $KStep, $isDHS, $localBgrWidth, $flankSize, $maxColNum, 1, $qValueThresh, $transformType, $isOneSideBkg, '$resDir', '$outPrefix');"
set +x

footprinting_generate_bed.py \
    $resDir/$finalDHSBed $resDir/$outPrefix.qv$qValueThresh \
    > $resDir/$outPrefix.temp.bed

grep -v 'chrM' $resDir/$outPrefix.temp.bed \
    | cut -f 1,2,3,4 \
    > $resDir/$ftprtsFile


rm $resDir/$outPrefix.datatemp
rm $resDir/$outPrefix.nulltemp
rm $resDir/$outPrefix.qv$qValueThresh
rm $resDir/$outPrefix.temp.bed
rm $resDir/$finalDHSBed
rm $resDir/$cutCountsMtx
rm $resDir/$mappableMtx
