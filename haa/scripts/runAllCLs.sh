for dir in `ls -d 0 201*_mm*_bern` ; 
do 
echo year and channel $(echo $dir)
yrch=`echo $dir | sed -e 's/_bern//'`
cd $dir 
echo running impacts on $yrch
#cp -r ../runImpacts.sh .
cp -r ../runAsymptoticLimits_extrpolated.sh .
#bash runImpacts.sh $yrch
bash runAsymptoticLimits_extrpolated.sh $yrch
cd ../.;
#; 
done

