if [ $# -eq 0 ]  
then
echo "$0 file1 file2 binsize  :  Generate a hsitogram from a raw data (one colume)"
echo "file1 = original file"
echo "file2 = output histogram file"
exit 1
fi
gawk -v binsize=$3 -v max=0 \
'{ x[int($1/binsize)]++; if (max < $1) max=$1; } \
 END{\
	maxbin = int(max/binsize);\
	for(i=0; i<=maxbin; i++){\
		if (x[i]==NULL) x[i]=0;\
		print binsize*i*1.5, x[i];\
	};\
    }' $1 > $2
