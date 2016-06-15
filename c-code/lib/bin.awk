if [ $# -eq 0 ]  
then
echo "$0 file1 file2 binsize  :  Resize the binsize and generate a new histogram"
echo "file1 = original file"
echo "file2 = output file corresponding to the new given binsize"
exit 1
fi
gawk -v binsize=$3 -v sumx=0 -v sumy=0 -v count=0 \
'{sumx+=$1; sumy+=$2; if(NR%binsize==0) {\
	print sumx/binsize, sumy; \
	sumx=0; sumy=0; \
	}}' $1 > $2
