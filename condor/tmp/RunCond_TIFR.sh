#!/bin/bash
#REFERENCE
#https://github.com/florez/CONDOR

ntupleT2Paths=$1

count=0
cat $ntupleT2Paths | while read ntupleT2Path
do
  ((count++))
  echo -e "\033[01;32m input ntuple=\033[00m" $count": " $ntupleT2Path
  IFS='/ ' read -r -a array <<< "$ntupleT2Path"
  len=${#array[@]}
  fifth_last=`expr $len - 4`
  sec_last=`expr $len - 1`
  ntuple=${array[$sec_last]}
  ##ntuple=${array[$fifth_last]}${array[$sec_last]}
  #echo $ntuple
  iFile=${ntuple/.root/""}
  #----------------------------------------------
  #copy condor scripts to each input ntuple dir
  #replace Cond.sub arguments, as per input
  #submit the condor jobs, for each ntuple
  #----------------------------------------------
  
  mkdir -p $iFile
  cp Cond.sub $iFile
  cp Analyzer_TIFR.sh $iFile
  cd $iFile 
  sed -i "s:FNAME:$ntupleT2Path:g" Cond.sub
  sed -i "s:OUTPUTFILE:$iFile:g" Cond.sub
  sed -i "s:OUTPUTDIR:$iFile:g" Cond.sub
  condor_submit Cond.sub
  cd ../
done
