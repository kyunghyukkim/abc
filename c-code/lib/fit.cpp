#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include "lsf.cpp"
#include <string.h>


int blank_line(FILE *pFile);
void  move_to_next_line(FILE *pFile);
void num_of_line(char *file, int index, double init, double fin, int *ndata, int *nint, int *nfin);

void help_message(void)
{
	printf("Usage fit.o number_of_column input_file index [i] [j] \n");
	printf("This function returns the least square fit.");
	printf("Fit function is l(x) = a+b*x.\n");
	printf("If the number of column is 2 then we consider ther data does not have errorbars and it 3, there are.\n");
	printf("index and [i:j] are the same meaning used in gnuplot fit function.\n");
	printf("If i and j are not given, the whole data in the indexed block will be fitted.");
	printf("fit.o prints out a dela b delb chi-square ......\n");
	exit(-1);
}

int main(int argc, char *argv[])
{
	int ndata;
	int n;
	int index;
	double init, fin;
	int ninit, nfin;
	double *x,*y, *sig;
	double a, b, siga, sigb;
	double chi2, q;
	int mwt;
	char c;
	int t;

	double ave;
	double error;
	int i,j;
	FILE *ifile;

	if(strcmp(argv[1],"--h")==0) help_message();
	
	if(strcmp(argv[1], "3")==0) {
		mwt=1;
	}
	else 
	{
		mwt=0;
	}


	
	sscanf(argv[1], "%d", &n);
	sscanf(argv[3], "%d", &index);
	if(argc == 6) {
		sscanf(argv[4], "%lf", &init);
		sscanf(argv[5], "%lf", &fin);
	}
	else {
		init=-1.0e308;
		fin=1.0e308;
	}
		
	num_of_line(argv[2], index, init, fin, &ndata, &ninit, &nfin);
	
	if((ifile=fopen( argv[2], "r"))==NULL) {
		printf( "File open error : %s\n\a", argv[2]);
		exit(-1);
	}

	x=(double *) calloc(ndata+1, sizeof(double));
	y=(double *) calloc(ndata+1, sizeof(double));
	sig = (double *) calloc(ndata+1, sizeof(double));
	
	t=1;
	
	while(t!=index){
		if(blank_line(ifile)) {t++;}
		else move_to_next_line(ifile);
	}
	for(i=0;i<ninit-1;i++) move_to_next_line(ifile);
	for (i=1;i<=ndata;i++){
		if(fscanf(ifile, "%lf", x+i)==EOF) {
			ndata=i-1;
			break;
		}
		if(fscanf(ifile, "%lf", y+i)==EOF) {
			ndata=i-1;
			break;
		}
		if(mwt) 	if(fscanf(ifile, "%lf", sig+i)==EOF) {
			ndata=i-1;
			break;
		}
	}
	
	fit(x, y, ndata, sig, mwt, &a, &b, &siga, &sigb, &chi2, &q); 

	printf("%.15e\t %.15e\t %.15e\t %.15e\t %lf\t %lf\n", a,siga, b,sigb,chi2,q);

	free(x);
	free(y);
	free(sig);
}

//Check the current line is blank or not and then move the file pointer to the initial starting point.
int blank_line(FILE *pFile)
{
	fpos_t position;
	char c;
	
	fgetpos(pFile, &position);
	do{
		c=fgetc(pFile);
	} while (c==' ' || c=='\t');
	if(c=='\n' || c==EOF ) return 1;
	else {
		fsetpos(pFile, &position);
		return 0;
	}
}

//move the next line.
void move_to_next_line(FILE *pFile)
{
	int c;
	do{
		c=getc(pFile);
	}while (c!='\n' && c!=EOF);
}
	
//index is the same for gnuplot
//using [init:fin] in fit function of gnuplot
//*ndata is # of data for fitting
//*ninit is the row # of the first data for fitting in the indexed subset data.
//*nfin is the row # of the last data.
void num_of_line(char *file, int index, double init, double fin, int *ndata, int *ninit, int *nfin)
{
	FILE *pFile;
	char c;
	int m = 0;
	int t=1;
	int i;
	double read;

	pFile=fopen (file,"r");
	while(t!=index){
		if(blank_line(pFile)) {t++;}
		else move_to_next_line(pFile);
	}
	do{
		if(blank_line(pFile)==0) {
			m++;
			fscanf(pFile, "%lf", &read);
			if(m==1){
				if(read>=init) {
					*ninit=1;
				}
				else *ninit=2;
			}
			
			if(m!=1) 
				if(read<init) *ninit=m+1;
			
			if(read<=fin) *nfin=m;
			move_to_next_line(pFile);
			
		}
		else {
			*ndata=*nfin-*ninit+1;
			return;
		}
	}while (1);
	fclose(pFile);
}
