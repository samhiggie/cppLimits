#python hToaaFitter_optimize_py2.py -i skimmed_2016_paperex_mmmt_combined.root -id mmmt_inclusive -o 2016_mmmt -ch mmmt -cards -od 2016_mmmt_post_sq  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 2
#python hToaaFitter_optimize_py2.py -i skimmed_2016_paperex_mmtt_combined.root -id mmtt_inclusive -o 2016_mmtt -ch mmtt -cards -od 2016_mmtt_post_sq  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 1
#python hToaaFitter_optimize_py2.py -i skimmed_2016_paperex_mmet_combined.root -id mmet_inclusive -o 2016_mmet -ch mmet -cards -od 2016_mmet_post_sq  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 1
#python hToaaFitter_optimize_py2.py -i skimmed_2016_paperex_mmem_combined.root -id mmem_inclusive -o 2016_mmem -ch mmem -cards -od 2016_mmem_post_sq  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 3
#
#python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmmt_combined.root -id mmmt_inclusive -o 2017_mmmt -ch mmmt -cards -od 2017_mmmt_post_sq  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 2
#python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmtt_combined.root -id mmtt_inclusive -o 2017_mmtt -ch mmtt -cards -od 2017_mmtt_post_sq  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 4
##python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmet_combined.root -id mmet_inclusive -o 2017_mmet -ch mmet -cards -od 2017_mmet_pos_sqt  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 1
#python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmet_combined.root -id mmet_inclusive -o 2017_mmet -ch mmet -cards -od 2017_mmet_post_sq  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 1
##python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmem_combined.root -id mmem_inclusive -o 2017_mmem -ch mmem -cards -od 2017_mmem_pos_sqt  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 2
##python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmem_combined.root -id mmem_inclusive -o 2017_mmem -ch mmem -cards -od 2017_mmem_post_sq  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 1
#python hToaaFitter_optimize_py2.py -i skimmed_2017_paperex_mmem_combined.root -id mmem_inclusive -o 2017_mmem -ch mmem -cards -od 2017_mmem_post_sq  -bkgt bernstein -bkgo 4 -irbkgt bernstein -irbkgo 2
#
#python hToaaFitter_optimize_py2.py -i skimmed_2018_paperex_mmmt_combined.root -id mmmt_inclusive -o 2018_mmmt -ch mmmt -cards -od 2018_mmmt_post_sq  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 2
#python hToaaFitter_optimize_py2.py -i skimmed_2018_paperex_mmtt_combined.root -id mmtt_inclusive -o 2018_mmtt -ch mmtt -cards -od 2018_mmtt_post_sq  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 1
##python hToaaFitter_optimize_py2.py -i skimmed_2018_paperex_mmet_combined.root -id mmet_inclusive -o 2018_mmet -ch mmet -cards -od 2018_mmet_pos_sqt  -bkgt bernstein -bkgo 2 -irbkgt bernstein -irbkgo 1
#python hToaaFitter_optimize_py2.py -i skimmed_2018_paperex_mmet_combined.root -id mmet_inclusive -o 2018_mmet -ch mmet -cards -od 2018_mmet_post_sq  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 3
##python hToaaFitter_optimize_py2.py -i skimmed_2018_paperex_mmem_combined.root -id mmem_inclusive -o 2018_mmem -ch mmem -cards -od 2018_mmem_pos_sqt  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 2
#python hToaaFitter_optimize_py2.py -i skimmed_2018_paperex_mmem_combined.root -id mmem_inclusive -o 2018_mmem -ch mmem -cards -od 2018_mmem_post_sq  -bkgt bernstein -bkgo 1 -irbkgt bernstein -irbkgo 1


## find the sort of principle value then double the coeff range on the command line 
#limit --input-file skimmed_2016_paperex_mmmt_combined.root  -o 2016_mmmt --channel mmmt -d 2016_mmmt_final  --bkgt poly --bkgo 2 --irbkgt poly --irbkgo 2  --runs 0 --coeRangeIr 26500 --coeRange 12400
#limit --input-file skimmed_2016_paperex_mmtt_combined.root  -o 2016_mmtt --channel mmtt -d 2016_mmtt_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0 --sfr 1.25 
#limit --input-file skimmed_2016_paperex_mmet_combined.root  -o 2016_mmet --channel mmet -d 2016_mmet_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0 
##limit --input-file skimmed_2016_paperex_mmem_combined.root  -o 2016_mmem --channel mmem -d 2016_mmem_final  --bkgt poly --bkgo 2 --irbkgt poly --irbkgo 3  --runs 0
#limit --input-file skimmed_2016_paperex_mmem_combined.root  -o 2016_mmem --channel mmem -d 2016_mmem_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0 --coeRangeIr 9500
#
#
##limit --input-file skimmed_2017_paperex_mmmt_combined.root  -o 2017_mmmt --channel mmmt -d 2017_mmmt_final  --bkgt poly --bkgo 2 --irbkgt poly --irbkgo 2  --runs 0 --coeRangeIr 1000000 --coeRange 18300
#limit --input-file skimmed_2017_paperex_mmmt_combined.root  -o 2017_mmmt --channel mmmt -d 2017_mmmt_final  --bkgt poly --bkgo 2 --irbkgt poly --irbkgo 2  --runs 0 --coeRangeIr 200 --coeRange 100
##limit --input-file skimmed_2017_paperex_mmtt_combined.root  -o 2017_mmtt --channel mmtt -d 2017_mmtt_final  --bkgt poly --bkgo 2 --irbkgt poly --irbkgo 4  --runs 0
#limit --input-file skimmed_2017_paperex_mmtt_combined.root  -o 2017_mmtt --channel mmtt -d 2017_mmtt_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0 --sfr 1.25
#limit --input-file skimmed_2017_paperex_mmet_combined.root  -o 2017_mmet --channel mmet -d 2017_mmet_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0
##limit --input-file skimmed_2017_paperex_mmem_combined.root  -o 2017_mmem --channel mmem -d 2017_mmem_final  --bkgt poly --bkgo 4 --irbkgt poly --irbkgo 2  --runs 0
#limit --input-file skimmed_2017_paperex_mmem_combined.root  -o 2017_mmem --channel mmem -d 2017_mmem_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0 --sfr 1.3 --coeRangeIr 11000
#
#
#limit --input-file skimmed_2018_paperex_mmmt_combined.root  -o 2018_mmmt --channel mmmt -d 2018_mmmt_final  --bkgt poly --bkgo 2 --irbkgt poly --irbkgo 2  --runs 0 --coeRangeIr 200 --coeRange 700
#
#limit --input-file skimmed_2018_paperex_mmtt_combined.root  -o 2018_mmtt --channel mmtt -d 2018_mmtt_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 1  --runs 0 --coeRangeIr 25
#limit --input-file skimmed_2018_paperex_mmet_combined.root  -o 2018_mmet --channel mmet -d 2018_mmet_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 2  --runs 0 --coeRange 1000.0
#limit --input-file skimmed_2018_paperex_mmem_combined.root  -o 2018_mmem --channel mmem -d 2018_mmem_final  --bkgt poly --bkgo 1 --irbkgt poly --irbkgo 2  --runs 0 --coeRange 60 --coeRangeIr 95



# find the sort of principle value then double the coeff range on the command line 
limit --input-file skimmed_2016_paperex_mmmt_combined.root  -o 2016_mmmt --channel mmmt -d 2016_mmmt_bern  --bkgt bernstein --bkgo 2 --irbkgt bernstein --irbkgo 2  --runs 0
limit --input-file skimmed_2016_paperex_mmtt_combined.root  -o 2016_mmtt --channel mmtt -d 2016_mmtt_bern  --bkgt bernstein --bkgo 1 --irbkgt bernstein --irbkgo 1  --runs 0 --coeRangeIr 50000
limit --input-file skimmed_2016_paperex_mmet_combined.root  -o 2016_mmet --channel mmet -d 2016_mmet_bern  --bkgt bernstein --bkgo 1 --irbkgt bernstein --irbkgo 1  --runs 0 --coeRange 100000 --coeRangeIr 100000
limit --input-file skimmed_2016_paperex_mmem_combined.root  -o 2016_mmem --channel mmem -d 2016_mmem_bern  --bkgt bernstein --bkgo 2 --irbkgt bernstein --irbkgo 2  --runs 0 


limit --input-file skimmed_2017_paperex_mmmt_combined.root  -o 2017_mmmt --channel mmmt -d 2017_mmmt_bern  --bkgt bernstein --bkgo 2 --irbkgt bernstein --irbkgo 2  --runs 0
limit --input-file skimmed_2017_paperex_mmtt_combined.root  -o 2017_mmtt --channel mmtt -d 2017_mmtt_bern  --bkgt bernstein --bkgo 1 --irbkgt bernstein --irbkgo 1  --runs 0 --coeRangeIr 50000
limit --input-file skimmed_2017_paperex_mmet_combined.root  -o 2017_mmet --channel mmet -d 2017_mmet_bern  --bkgt bernstein --bkgo 1 --irbkgt bernstein --irbkgo 1  --runs 0 --coeRange 50000 --coeRangeIr 100000
limit --input-file skimmed_2017_paperex_mmem_combined.root  -o 2017_mmem --channel mmem -d 2017_mmem_bern  --bkgt bernstein --bkgo 2 --irbkgt bernstein --irbkgo 2  --runs 0 --coeRange 18000


limit --input-file skimmed_2018_paperex_mmmt_combined.root  -o 2018_mmmt --channel mmmt -d 2018_mmmt_bern  --bkgt bernstein --bkgo 2 --irbkgt bernstein --irbkgo 2  --runs 0
limit --input-file skimmed_2018_paperex_mmtt_combined.root  -o 2018_mmtt --channel mmtt -d 2018_mmtt_plots  --bkgt bernstein --bkgo 1 --irbkgt bernstein --irbkgo 1  --runs 0 --coeRange 65000 --coeRangeIr 45000
limit --input-file skimmed_2018_paperex_mmet_combined.root  -o 2018_mmet --channel mmet -d 2018_mmet_plots  --bkgt bernstein --bkgo 1 --irbkgt bernstein --irbkgo 1  --runs 0 --coeRange 13000 --sfr 2 --coeRangeIr 90000 

limit --input-file skimmed_2018_paperex_mmem_combined.root  -o 2018_mmem --channel mmem -d 2018_mmem_bern  --bkgt bernstein --bkgo 2 --irbkgt bernstein --irbkgo 2  --runs 0
