#!/bin/bash
if [ -z $1 ]
then
    mainout='AMass_gauss'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi 

# echo "running combine for all mass points on channel "
# combine -M AsymptoticLimits -m 20 -n _mass_a20_mmtt --run blind datacard_full_mass_a20_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 25 -n _mass_a25_mmtt --run blind datacard_full_mass_a25_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 30 -n _mass_a30_mmtt --run blind datacard_full_mass_a30_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 35 -n _mass_a35_mmtt --run blind datacard_full_mass_a35_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 40 -n _mass_a40_mmtt --run blind datacard_full_mass_a40_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 45 -n _mass_a45_mmtt --run blind datacard_full_mass_a45_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 50 -n _mass_a50_mmtt --run blind datacard_full_mass_a50_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 55 -n _mass_a55_mmtt --run blind datacard_full_mass_a55_${mainout}_mmtt.txt
# combine -M AsymptoticLimits -m 60 -n _mass_a60_mmtt --run blind datacard_full_mass_a60_${mainout}_mmtt.txt


# echo "hadding the files"
# rm higgsCombine_aa_mmtt_best.root
# hadd higgsCombine_aa_mmtt_best.root higgsCombine_mass_a*_mmtt.AsymptoticLimits.mH*.root
#
combine -M AsymptoticLimits -m 20 -n _mass_a20_mmmt --run blind datacard_full_mass_a20_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 25 -n _mass_a25_mmmt --run blind datacard_full_mass_a25_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 30 -n _mass_a30_mmmt --run blind datacard_full_mass_a30_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 35 -n _mass_a35_mmmt --run blind datacard_full_mass_a35_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 40 -n _mass_a40_mmmt --run blind datacard_full_mass_a40_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 45 -n _mass_a45_mmmt --run blind datacard_full_mass_a45_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 50 -n _mass_a50_mmmt --run blind datacard_full_mass_a50_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 55 -n _mass_a55_mmmt --run blind datacard_full_mass_a55_${mainout}_mmmt.txt
combine -M AsymptoticLimits -m 60 -n _mass_a60_mmmt --run blind datacard_full_mass_a60_${mainout}_mmmt.txt


echo "hadding the files"
rm higgsCombine_aa_mmmt_best.root
hadd higgsCombine_aa_mmmt_best.root higgsCombine_mass_a*_mmmt.AsymptoticLimits.mH*.root
root -l -b -q 'plotLimit.C+("aa","mmmt",2)'

#
#combine -M AsymptoticLimits -m 20 -n _mass_a20_mmet --run blind datacard_full_mass_a20_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 25 -n _mass_a25_mmet --run blind datacard_full_mass_a25_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 30 -n _mass_a30_mmet --run blind datacard_full_mass_a30_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 35 -n _mass_a35_mmet --run blind datacard_full_mass_a35_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 40 -n _mass_a40_mmet --run blind datacard_full_mass_a40_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 45 -n _mass_a45_mmet --run blind datacard_full_mass_a45_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 50 -n _mass_a50_mmet --run blind datacard_full_mass_a50_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 55 -n _mass_a55_mmet --run blind datacard_full_mass_a55_${mainout}_mmet.txt
#combine -M AsymptoticLimits -m 60 -n _mass_a60_mmet --run blind datacard_full_mass_a60_${mainout}_mmet.txt


#echo "hadding the files"
#rm higgsCombine_aa_mmet_best.root
#hadd higgsCombine_aa_mmet_best.root higgsCombine_mass_a*_mmet.AsymptoticLimits.mH*.root

# combine -M AsymptoticLimits -m 20 -n _mass_a20_mmem --run blind datacard_full_mass_a20_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 25 -n _mass_a25_mmem --run blind datacard_full_mass_a25_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 30 -n _mass_a30_mmem --run blind datacard_full_mass_a30_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 35 -n _mass_a35_mmem --run blind datacard_full_mass_a35_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 40 -n _mass_a40_mmem --run blind datacard_full_mass_a40_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 45 -n _mass_a45_mmem --run blind datacard_full_mass_a45_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 50 -n _mass_a50_mmem --run blind datacard_full_mass_a50_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 55 -n _mass_a55_mmem --run blind datacard_full_mass_a55_${mainout}_mmem.txt
# combine -M AsymptoticLimits -m 60 -n _mass_a60_mmem --run blind datacard_full_mass_a60_${mainout}_mmem.txt


# #echo "hadding the files"
# rm higgsCombine_aa_mmem_best.root
# hadd higgsCombine_aa_mmem_best.root higgsCombine_mass_a*_mmem.AsymptoticLimits.mH*.root

