if [ $# -eq 0 ]  
then
echo "$0 file1 block-size  :  Coarse-grain-average of the coarsening block size"
echo "file1 = original file"
exit 1
fi
gawk -v binsize=$2 -v sumx=0 -v sumy=0 -v count=0 \
'{sumx+=$1; sumy+=$2; if(NR%binsize==0) {\
	print sumx/binsize, sumy/binsize; \
	sumx=0; sumy=0; \
	}}' $1 
