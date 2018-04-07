echo  "Calculating limit for mass of Excited Lepton $1"
#combineCards.py  datacard_ele_csbar_bbb__mH150.txt datacard_mu_csbar_bbb__mH150.txt  > comb_mH150.txt
text2workspace.py datacard_csbar_mu_mjj_max_13TeV_mH$1.txt -P HiggsAnalysis.CombinedLimit.ChargedHiggs:brChargedHiggs -o comb_mH$1.root
combine --rAbsAcc 0.00001 comb_mH$1.root -M Asymptotic --mass $1 --name ExitedLepton_mu
