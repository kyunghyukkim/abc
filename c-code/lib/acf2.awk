if [ $# -eq 0 ]  
then
echo "$0 file block_size: autocorrelation file ("file") is convered to the autocorrelation of blocked averaged signals with the given block size.  "
exit 1
fi
gawk -v binsize=$2 -v sumx=0 -v sumy=0 -v count=0 \
'{
	x=(NR-1)%binsize;\
	if(x==0) sumy+=$2*binsize;\
	else sumy+=$2*(binsize-x)*2;\
	sumx+=$1;\
	if(NR%binsize==0) {\
		print sumx/binsize, sumy; \
		sumx=0; sumy=0; \
	}}' $1 

