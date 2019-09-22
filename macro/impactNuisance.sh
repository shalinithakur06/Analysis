t2wDataCard=$1
mass=$2

impactDir='outImpact'
mkdir -p $impactDir
cp $t2wDataCard $PWD/$impactDir/
t2wFile=$PWD/$impactDir/$t2wDataCard

cd $impactDir
#combineTool.py -M Impacts -d $t2wFile -m $mass --doInitialFit --robustFit 1 --setParameterRanges r=0.0,0.10 
#combineTool.py -M Impacts -d $t2wFile -m $mass --robustFit 1 --setParameterRanges r=0.0,0.10  --doFit  --task-name nuisImpact --parallel 10 | tee impactLog.txt
combineTool.py -M Impacts -d $t2wFile -m $mass --doInitialFit --robustFit 1 
combineTool.py -M Impacts -d $t2wFile -m $mass --robustFit 1 --doFit  --task-name nuisImpact --parallel 10 | tee impactLog.txt
combineTool.py -M Impacts -d $t2wFile -m $mass -o nuisImpactJSON 
plotImpacts.py --cms-label "   Internal" -i nuisImpactJSON -o nuisImpact
