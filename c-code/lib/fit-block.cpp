#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include "lsf.cpp"
#include <string.h>
#include "average.cpp"
int blank_line(FILE *pFile);
void  move_to_next_line(FILE *pFile);
void num_of_line(char *file, int width, double init, double fin, int *ndata, int *nint, int *nfin);

void help_message(void)
{
	printf("Usage fit-block.o number_of_column input_file init_value_of_1st_column block_width \n");
	printf("This function returns the least square fit of blocks with their width given by 'block-width'.\n");
	printf("The first block starts from init_value_of_1st_column.\n");
	printf("Fit function is l(x) = a+b*x.\n");
	printf("If the number of column is 2 then we consider ther data do not have errorbars and if 3, there are.\n");
	printf("fit-block.o prints out mean_of_1st_column_block a dela b delb chi-square ......\n");
	exit(-1);
}

int main(int argc, char *argv[])
{
	int ndata;
	int n;
	double width;
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
	double minx;
	double errx;

	//There is no blank line separating data in the data file.
	int index=1;
	
	if(strcmp(argv[1],"--h")==0) help_message();
	
	if(strcmp(argv[1], "3")==0) {
		mwt=1;
	}
	else 
	{
		mwt=0;
	}


	
	sscanf(argv[1], "%d", &n);
	sscanf(argv[3], "%lf", &width);
	sscanf(argv[4], "%lf", &init);
	fin=init+width;
	
	printf("#minx\t y_int\t y_int_err\t slope\t slope_err\t chi2\t q\n");
	while(1){
		num_of_line(argv[2], index, init, fin, &ndata, &ninit, &nfin);
		//Terminate when there is no data to fit, because the access point has already reached EOF or the block size width is too small.
		if(ndata==0) exit(-1);


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
		ave_error1(ndata, x, &minx, &errx);

		//printf("%f\t %.15e\t %.15e\t %.15e\t %.15e\t %lf\t %lf\n", minx, a,siga, b,sigb,chi2,q);
		printf("%f\t %.15e\t %.15e\t %.15e\t %.15e\t %lf\t %lf\n", minx, a,siga, b,sigb,chi2,q);

		free(x);
		free(y);
		free(sig);
		init+=width;
		fin+=width;
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
