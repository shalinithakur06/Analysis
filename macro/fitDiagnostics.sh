t2wDataCard=$1
mass=$2

impactDir='outFitDiag'
mkdir -p $impactDir
cp $t2wDataCard $PWD/$impactDir/
t2wFile=$PWD/$impactDir/$t2wDataCard

cd $impactDir
combine $t2wFile -M FitDiagnostics -t -1 --expectSignal 1 --rMin 0.00001
diffNuFilePath=${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py
python $diffNuFilePath -a $PWD/fitDiagnostics.root -g fitDiagnostics
