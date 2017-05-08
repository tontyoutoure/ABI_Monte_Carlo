#!/bin/bash

MD5LENGTH=8
NUMBER_OF_THREADS=7
THICKNESS=0.0041

#./light -x 0.002 0.012 -y 0.002 0.012 -n 10000 --mode-pp -o input.pho
#rm output.pho
#./REF_fast_exp -i input.pho -t $TISSUE_PATH/1b.tissue -o output.pho -a 0.01 -b 0.01 -c ${THICKNESS} -r 0.002 -l run.log --number-of-threads $NUMBER_OF_THREADS -e 20 --index-count 1

rm output.pho
for ((i=0; i<500; i=i+1)); do
	echo The ${i}th file is proceeding
	./light -x 0.004 0.014 -y 0.004 0.014 -n 10000000 --mode-pp -o input.pho
	./REF -i input.pho -t ../tissues/1b.tissue -o output.pho -a 0.01 -b 0.01 -c 0.1 -r 0.002 -l run.log --number-of-threads $NUMBER_OF_THREADS -n 0.999999963960949 -L -0.085681550167100 --index-count 1
	MD5=$(md5sum output.pho|head -c $MD5LENGTH)
	rm input.pho
	mv output.pho /media/tontyoutoure/MobileDisk1/pho/1b/output_$MD5.pho
done

#for ((i=0; i<1000; i=i+1)); do
#	echo The ${i}th file is proceeding
#	./light -x 0.0002 0.0102 -y 0.0002 0.0102 -n 10000000 --mode-pp -o input.pho
#	./REF_fast_exp -i input.pho -t tissues/*tissue -o output.pho -a 0.01 -b 0.01 -c ${THICKNESS} -r 0.0001 -l run.log --number-of-threads $NUMBER_OF_THREADS -e 60 --index-count 500
#	MD5=$(md5sum output.pho|head -c $MD5LENGTH)
#	rm input.pho
#	mv output.pho output_$MD5.pho
#	mv output_$MD5.pho /media/tontyoutoure/MobileDisk1/pho/one_layer
#done
