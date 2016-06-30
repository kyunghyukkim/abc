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

void initialize(FILE *ifile, double *x, double *y, int ndata);
void autocorr(double *S, double *acf, int, int);
void mean(double *S, double *, int);
void help_message(void)
{
	printf("Usage acf.o input_file number_of_autocorrelation_data_points\n");
	printf("This function returns the autocorrelatin of the second column of input_file, where the first column is time.\n");
	exit(-1);
}



int main(int argc, char *argv[])
{
	int ndata;
	int n;
	int index;
	double init, fin;
	int ninit, nfin;
	double *x, *acf;
	double *y;
	char c;
	int t;
	double meanvalue;

	int i,j;
	FILE *ifile;
	int ncor;

	if(strcmp(argv[1],"--h")==0) help_message();


	init=-1.0e308;
	fin=1.0e308;

	index =1;

	num_of_line(argv[1], index , init, fin, &ndata, &ninit, &nfin);
//	printf("%d", ndata);
	
	ifile=fopen(argv[1], "r");
	sscanf(argv[2], "%d", &ncor);
	
	x=(double *) calloc(ndata+1, sizeof(double));
	y=(double*) calloc(ndata+1, sizeof(double));
	acf = (double *) calloc(ndata+1, sizeof(double));
	
	for (i=1;i<=ndata; i++){
		fscanf(ifile, "%lf %lf", x+i, y+i); 
	}

//	for(int i=1;i<=ndata;i++){
//		printf("%f\t%d\n", x[i], y[i]);}
	
	mean(y, &meanvalue, ndata);
	//printf ("\n%f\n\n\n\n", meanvalue);
	autocorr(y, acf, ncor, ndata);
	
	double  delt = x[2]-x[1];
	for(i=0;i<ncor;i++){
		printf("%f\t%f\n", delt*i, acf[i]-meanvalue*meanvalue);
	}

	free(x);
	free(y);
	free(acf);
}

void mean(double *S, double *m, int ndata)
{
	*m=0;
	int counter=0;

	for(int i=0;i<=ndata;i++){
		*m+=S[i];
		counter++;
	}
	*m=*m/(double)counter;
}

void autocorr(double *S, double *acf, int ncor, int ndata)
{
	int t;
	double a0;
	int Qhead;
	double *x_c;
	int *counter;
	int i,j;

	x_c = (double *) calloc(ndata+1, sizeof(double));
	counter = (int *) calloc(ndata+1, sizeof(int));

	for (i=0;i<=ndata;i++) counter[i] = 0;

	for(int i=1;i<=ndata; i++){	
		a0  = S[i];
		for(j=0;j<=ncor && i+j<=ndata;j++){
		       	acf[j]+=a0*S[i+j];
			counter[j]++;
		}

	}
	
	for (i=0;i<=ncor;i++) {
		acf[i] = acf[i]/(double)counter[i];
	}
	

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
