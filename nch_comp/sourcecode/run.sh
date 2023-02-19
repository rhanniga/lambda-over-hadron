#!/bin/bash

###########################
## Author : Beomkyu Kim  ##
## email  : kimb@cern.ch ##
###########################
gridwork="full"

if [[ "$#" != 3 ]]
then
  echo "Wrong usage"
  echo "usage : $0 <taskname> <Period> <full|terminate|download|merge>"
  ##echo "usage : $0 <full|terminate>"
  exit 0;
fi

##taskname : task name
taskname=$1
##periods : remove some periods if you don't want to run all of them
periods=$2
method=$3

basepath=/alice/cern.ch/user/r/rhanniga/
DOWN_DIR=~/Desktop
currentdir=$PWD
DATA_DIR=$DOWN_DIR/${basepath}/${taskname}${periods};


function download {
alien-token-init 
cd $DOWN_DIR 
perl ${ALICE_PHYSICS}/PWGUD/DIFFRACTIVE/macros/alien_cp.pl ${basepath}/${taskname}${periods}/ root_archive.zip AnalysisResults.root
}

function merge_list {
outname=${1:-temp_out.root}
np=${2:-1} # Number Of Process
tag=tmp-$outname-$(date +%s)-$RANDOM
mkdir -p $tag
xargs -n25 | perl -ne'print "'$tag'/${.}.root $_"' | xargs  -P$np -L1   hadd -f   
#xargs -n25 | perl -ne'print "${.}.root $_"' | xargs -P$np -I% bash -c 'echo hadd -f '$tag/'%'  
hadd -f $outname $tag/*.root
echo $outname
rm $tag/*.root && rmdir $tag
ls
pwd
}
export -f merge_list


######################
#  MERGE Run by Run
#####################
function merge_RunByRun {
cd $DATA_DIR
rm AnalysisResults*.root
ls -d out/* | xargs -P4 -I% bash -c 'hadd -f    AnalysisResults_'${taskname}${periods}'_$(basename %).root $(find % -name AnalysisResults.root)'
hadd -f   AnalysisResults_${taskname}${periods}.root AnalysisResults_${taskname}${periods}_*.root
echo "#------------------------------"
echo "#      Merged Files "
cd - > /dev/null
rm -r  ~/Desktop/${taskname}${periods}
mkdir -p ~/Desktop/${taskname}${periods}
mv $DATA_DIR/AnalysisResults*.root  ~/Desktop/${taskname}${periods}/
echo "#   ./merged/${taskname}${periods}/AnalysisResults_${taskname}${periods}.root"
echo "#------------------------------"
}


if [ $method = "full" ]
then
  		for i in ${periods}
  		do 
    		root -l -b -q runROOT6.C\(\"${taskname}\",\"${i}\",\"${gridwork}\",\"grid\"\)
    		rm ${taskname}*
    		rm *.d *.so
    		rm *.root
    		rm myAnalysis.C
				rm AutoDict*
  		done
		##sleep 300
		##./resubmit_alien.sh

elif [ $method = "download" ]
then
  download
elif [ $method = "merge" ]
then
  merge_RunByRun
elif [ $method = "terminate" ]
then
  download
  merge_RunByRun
fi

##cd $currentdir

