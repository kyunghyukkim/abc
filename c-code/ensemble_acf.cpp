//Calculate covariance of J.
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdlib.h>
//#include <malloc.h>
#include <stdio.h>
#include <time.h>
#include <cstring>
#include "./lib/mt19937ar-cok.cpp"
#include "./lib/corfun.cpp"
#include "./lib/average.cpp"

//using namespace std;
#define MAX(c1,c2) ((c1>c2)? c1:c2)

long  *S, *INIT_S;
double  *C;


double TIME;
double DELT; //time increment in each MC step.
double  *aveS, *stdS;
double  *sumS, *sumS2;
double  *aveJ, **covJ;
double *sumJ, **sumJ2;
long counter;
long **counterS;

double *J;
double *v;
double *counterv;
int **stoi;
int N;

int NUM_SPEC;
int NUM_SPEC2;
int NUM_PARAM;
int NUM_REACT;

double GRID_CONST;
long NGRID_BF_STEADY;
long NGRID_AT_STEADY;
int NENSEMBLE;
int TCORGRID;  // = TCOR/GRID_CONST
double TCOR;
int **x_c;  //counter for autocorrelation with a given time distance
double ***time_acf; //non-normalized acf
double ***time_acf_ave;
int Qhead;

long INDEX_GRID;
long *SMAX;

void model_def(void);
void init_state(void);
bool MC_not_finished(int);
void time_ave(void);
void time_evol(void);
void MC_results(void);
void resize_s_counter(int index_species, int new_copy_numb_of_moleculesi);
void propensity(void);
int determine_reaction(double *sv, double rand);
void create_output_files(void);
void cross_pair_sum(int head, int  *a, int *b, double *acf);
unsigned long  TCORGRID_mult_of_2(void);
void fourier(double *acf, double *onepsd, unsigned long tcorgrid);

std::ofstream FILE_TEVOL;
std::ofstream FILE_FLUX;








int main(void)
{
	long nextseed=0;
	model_def();
	create_output_files();
	
	time_evol();
	time_ave();
	MC_results();
	
}

void create_output_files(void)
{
	FILE *temp;
	std::ofstream myfile;


	char filename[20];
	std::string mystr;

	for(int i=1;i<=NUM_SPEC;i++){
		std::ostringstream oss;
		oss<<"prob-s"<<i;
		mystr= oss.str();
		strcpy(filename, mystr.c_str());
		myfile.open(filename);
		myfile<<"#S\t Probability of S"<<i<<"\n";
		myfile.close();
	}

	myfile.open("ave-std");
	myfile<<"#S[i]\t ave \t std\n";
	myfile.close();

	myfile.open("tevol");
	myfile<<"#Time\t";
	for(int i=1;i<=NUM_SPEC;i++)
		myfile<<"S["<<i<<"]\t";
	myfile<<std::endl;
	myfile.close();
}
	
void model_def(void)
{
	int i,j ;

	NUM_SPEC=3;
	NUM_SPEC2=NUM_SPEC*NUM_SPEC;
	NUM_PARAM=4;
	NUM_REACT=4;

	FILE *temp;
	S=(long *) calloc(NUM_SPEC+1, sizeof(long));
	INIT_S=(long *)calloc(NUM_SPEC+1, sizeof(long));
	aveS=(double *)calloc(NUM_SPEC+1, sizeof(double));
	stdS=(double *)calloc(NUM_SPEC+1, sizeof(double));
	sumS=(double *)calloc(NUM_SPEC+1, sizeof(double));
	sumS2=(double *)calloc(NUM_SPEC+1, sizeof(double));
	J=(double *)calloc(NUM_REACT+1, sizeof(double));
	aveJ=(double *)calloc(NUM_REACT+1, sizeof(double));
	sumJ=(double *)calloc(NUM_REACT+1, sizeof(double));
	//	flux covariance matrix
	covJ=(double **)calloc(NUM_REACT+1, sizeof(double *));
	sumJ2=(double **)calloc(NUM_REACT+1, sizeof(double *));
	for(i=1;i<=NUM_REACT;i++) {
		covJ[i]=(double *)calloc(NUM_REACT+1, sizeof(double));
		sumJ2[i]=(double *)calloc(NUM_REACT+1, sizeof(double));
	}
	C=(double *)calloc(NUM_PARAM+1, sizeof(double));
	counterS=(long **)calloc(NUM_SPEC+1, sizeof(long *));
	counterv=(double *)calloc(NUM_REACT+1, sizeof(double));
	SMAX=(long *)calloc(NUM_SPEC+1, sizeof(long));
	FILE *fp;
	int d1[10];
	fp=fopen("info-acf","r");
	for(i=1;i<=NUM_SPEC;i++) {
		fscanf(fp, "%s %d", d1, &INIT_S[i]);
	}

	for(i=1;i<=NUM_PARAM;i++) fscanf(fp, "%s %lf", d1, &C[i]);

	fscanf(fp, "%s %lf", d1, &GRID_CONST);
	fscanf(fp, "%s %ld", d1, &NGRID_BF_STEADY);
	fscanf(fp, "%s %ld", d1, &NGRID_AT_STEADY);
	fscanf(fp, "%s %lf", d1, &TCOR);

	fclose(fp);
	TCORGRID=(int)(TCOR/GRID_CONST)+1;
	x_c=(int **) calloc(NUM_SPEC+1, sizeof(int *));
	for(i=1;i<=NUM_SPEC;i++)
		x_c[i]=(int *) calloc(TCORGRID, sizeof(int));
	time_acf=(double ***) calloc(NUM_SPEC+1, sizeof(double **));//Not normalized
	time_acf_ave=(double ***) calloc(NUM_SPEC+1, sizeof(double **));
	for(i=1;i<=NUM_SPEC;i++){
		SMAX[i]=INIT_S[i];
		counterS[i]=(long *) calloc( SMAX[i]+1, sizeof(long));
		if(SMAX[i]<S[i]) resize_s_counter(i, S[i]);
	}
	for(i=1;i<=NUM_SPEC;i++){
		time_acf[i]=(double **) calloc( NUM_SPEC+1, sizeof(double *));
		time_acf_ave[i]=(double **)calloc (NUM_SPEC+1, sizeof(double *));
	}
	for(i=1;i<=NUM_SPEC;i++){
		for(j=1;j<=NUM_SPEC;j++){
			time_acf[i][j]=(double *) calloc(TCORGRID, sizeof(double));
			time_acf_ave[i][j]=(double *)calloc (TCORGRID, sizeof(double));
		}
	}
	//	propensity function
	v=(double *)calloc(NUM_REACT+1, sizeof(double));

	//	Stoichiometry matrix
	stoi=(int **)calloc(NUM_SPEC+1, sizeof(int *));
	for(i=1;i<=NUM_SPEC;i++) 
		stoi[i]=(int *)calloc(NUM_REACT+1, sizeof(int));

	for(i=1;i<=NUM_SPEC;i++)for(j=1;j<=NUM_REACT;j++) stoi[i][j]=0;


	propensity();

	stoi[1][1]=1;
	stoi[1][2]=-1;
	stoi[1][3]=-1;
	stoi[1][4]=1;

	stoi[2][3]=-1;
	stoi[2][4]=1;
	stoi[3][3]=1;
	stoi[3][4]=-1;

}

void propensity(void)
{
	v[1]=C[1];
	v[2]=C[2]*S[1];
	v[3]=C[3]*S[1]*S[2];
	v[4]=C[4]*S[3];
}




void time_evol(void)
{
	int i;
	
	init_genrand(time(NULL)+genrand_int32());
	init_state();
		

	FILE_TEVOL.open("tevol", std::ios::out);
	FILE_FLUX.open("v_t", std::ios::out);
	
	
	for(i=1;i<=NUM_REACT;i++) counterv[i]=0;
	while(INDEX_GRID < NGRID_BF_STEADY )  MC_not_finished(1);
	FILE_TEVOL.close();
	FILE_FLUX.close();
}





void time_ave(void)
{
	long l;
	int i,j,k;
	double weight;
	long ngrid_total=NGRID_BF_STEADY+NGRID_AT_STEADY;


	init_genrand(time(NULL)+genrand_int32());
	init_state();




	for(i=1;i<=NUM_SPEC;i++) {
		sumS[i]=0;
		sumS2[i]=0;
	}
	for(i=1;i<=NUM_REACT;i++) {
		sumJ[i]=0;
		for(j=1;j<=NUM_REACT;j++){
			sumJ2[i][j]=0;
		}
	}
	for(i=1;i<=NUM_SPEC;i++)
		for(k=0;k<TCORGRID;k++)
			x_c[i][k]=0;
	for(i=1;i<=NUM_SPEC;i++){
		for(j=1;j<=NUM_SPEC;j++){
			for(k=0;k<TCORGRID;k++){
				time_acf[i][j][k]=0;
				time_acf_ave[i][j][k]=0;
			}
		}
	}
	counter=0;
	for(i=1;i<=NUM_SPEC;i++) 
		for(j=0;j<=SMAX[i];j++)	counterS[i][j]=0;

	for(i=1;i<=NUM_REACT;i++) counterv[i]=0;

	for(i=1;i<=NUM_SPEC;i++) 
		for(j=1;j<=NUM_SPEC;j++) 
			Qhead=0;
	while(INDEX_GRID < ngrid_total )  {
		MC_not_finished(0);

	}

}






void MC_results(void)
{
	int i, j, k;

	for(i=1;i<=NUM_SPEC;i++){
		aveS[i]=(double)sumS[i]/counter;
		stdS[i]=sqrt((double)sumS2[i]/counter-aveS[i]*aveS[i]);
	}

	for(i=1;i<=NUM_REACT;i++){
		aveJ[i]=(double)sumJ[i]/counter;
	}
	for(i=1;i<=NUM_REACT;i++)
		for(j=1;j<=NUM_REACT;j++){
			covJ[i][j]=(double)sumJ2[i][j]/counter-aveJ[i]*aveJ[j];
		}

	char filename[20];	
	std::string mystr;
	std::ofstream myfile;	
	for(int i=1;i<=NUM_SPEC;i++){
		std::ostringstream oss;
		oss<<"prob-s"<<i;
		mystr= oss.str();		
		strcpy(filename, mystr.c_str());
		myfile.open(filename, std::ios::out);
		for(j=0;j<=SMAX[i];j++) 
			myfile<<j<<"\t"<<(double)counterS[i][j]/counter<<std::endl;
		myfile.close();
	}


	for(int i=1;i<=NUM_SPEC;i++){
		std::ostringstream oss;
		oss<<"s"<<i;
		mystr= oss.str();		
		strcpy(filename, mystr.c_str());
		myfile.open(filename, std::ios::out);
		for(j=0;j<=SMAX[i];j++) 
			myfile<<j<<"\t"<<(double)counterS[i][j]/counter<<std::endl;
		myfile.close();
	}

	myfile.open("s-ave-std", std::ios::out);
	for(i=1;i<=NUM_SPEC;i++)
		myfile<<i<<"\t"<<aveS[i]<<"\t"<<stdS[i]<<std::endl;
	myfile.close();

	myfile.open("J-ave-std", std::ios::out);
	for(i=1;i<=NUM_REACT;i++)
		myfile<<i<<"\t"<<aveJ[i]<<"\t"<<sqrt(covJ[i][i])<<std::endl;
	myfile.close();

	myfile.open("J-cov", std::ios::out);
	myfile<<"#row_index column_index J-cov"<<std::endl;
	for(i=1;i<=NUM_REACT;i++)
		for(j=1;j<=NUM_REACT;j++)
			myfile<<i<<"\t"<<j<<"\t"<<covJ[i][j]<<std::endl;
	myfile.close();

	for(i=1;i<=NUM_SPEC;i++){
		for(j=1;j<=NUM_SPEC;j++){
			std::ostringstream oss;
			oss<<"acf"<<i<<"_"<<j;
			mystr= oss.str();		
			strcpy(filename, mystr.c_str());
			myfile.open(filename, std::ios::out);
			myfile<<"#time"<<"\t"<<"acf"<<std::endl;
			for(k=0;k<TCORGRID;k++){
				time_acf_ave[i][j][k]=(double)time_acf[i][j][k]/(counter-k)-aveS[i]*aveS[j];
				myfile<<k*GRID_CONST<<"\t"<<time_acf_ave[i][j][k]<<std::endl;
			}
			myfile.close();
		}
	}

	system("paste acf1_1 acf1_3 acf3_1 acf3_3 | awk '{print $1, $2+$4+$6+$8}' > acfy_y");
/*
	unsigned long n_point_fourier;
	n_point_fourier	= TCORGRID_mult_of_2();
	double ***onepsd;
	onepsd=(double ***) calloc(NUM_SPEC+1, sizeof(double **));
	for(i=1;i<=NUM_SPEC;i++) 
		onepsd[i] = (double **) calloc(NUM_SPEC+1, sizeof(double *));
	for(i=1;i<=NUM_SPEC;i++) 
		for(j=i;j<=NUM_SPEC;j++) 
			onepsd[i][j]=(double *) calloc(n_point_fourier/2+1, sizeof(double));

	for(i=1;i<=NUM_SPEC;i++){
		for(j=i;j<=NUM_SPEC;j++){
			std::ostringstream oss;
			oss<<"psd"<<i<<"_"<<j;
			mystr= oss.str();		
			strcpy(filename, mystr.c_str());
			myfile.open(filename, std::ios::out);

			fourier(time_acf_ave[i][j], onepsd[i][j], n_point_fourier);
			myfile<<"#frequency"<<"\t"<<"psd"<<std::endl;
			for(k=0;k<=n_point_fourier/2;k++)
				myfile<<1./n_point_fourier/GRID_CONST*k<<"\t"<<onepsd[i][j][k]<<std::endl;
			myfile.close();
		}

		}
*/
}

unsigned long  TCORGRID_mult_of_2(void)
{
	unsigned long tcorgrid;

	tcorgrid=(unsigned long) TCORGRID;
	if(tcorgrid<2) tcorgrid=0;
	else if (tcorgrid<4) tcorgrid=2;
	else if (tcorgrid<8) tcorgrid=4;
	else if (tcorgrid<16) tcorgrid=8;
	else if (tcorgrid<32) tcorgrid=16;
	else if (tcorgrid<64) tcorgrid=32;
	else if (tcorgrid<128) tcorgrid=64;
	else if (tcorgrid<256) tcorgrid=128;
	else if (tcorgrid<512) tcorgrid=256;
	else if (tcorgrid<1024) tcorgrid=512;
	else if (tcorgrid<2048) tcorgrid=1024;
	else if (tcorgrid<4096) tcorgrid=2048;
	else if (tcorgrid<8192) tcorgrid=4096;
	else if (tcorgrid<8192*2) tcorgrid=4096*2;
	else if (tcorgrid<8192*2*2) tcorgrid=4096*2*2;
	else {
		std::cout<< "TCORGRID is too large so increase GRID_CONT"<<std::endl;
		exit(-1);	
	}
	return tcorgrid;
}



	



bool  MC_not_finished(int time_evolution)
{
	double sv[NUM_REACT];
	double sum_rate=0;
	double rand;
	int index;
	int i,j;

	propensity();
	for(i=1;i<=NUM_REACT;i++) sum_rate+=v[i];

	if(sum_rate!=0) {
		DELT=-log(genrand_real3())/sum_rate;
		for(i=1;i<=NUM_REACT;i++) sv[i]=v[i]/sum_rate;
	}
	//if reaction finishes, then new DELT is going to be the last DELT just before reaction finishes.
	TIME+=DELT;
	
	if(time_evolution==1){
		while(TIME>INDEX_GRID*GRID_CONST) {
			FILE_TEVOL<<INDEX_GRID*GRID_CONST<<"\t";
			for(i=1;i<=NUM_SPEC;i++)
				FILE_TEVOL<<S[i]<<"\t";
			FILE_TEVOL<<std::endl;

			FILE_FLUX << INDEX_GRID* GRID_CONST <<"\t";
			for(i=1;i<=NUM_REACT;i++) 
				FILE_FLUX << counterv[i]/GRID_CONST << "\t";
			FILE_FLUX<<std::endl;
			for(i=1;i<=NUM_REACT;i++) counterv[i]=0;


			if(++INDEX_GRID>NGRID_BF_STEADY) return 0;
		}

		//If reaction is not finished, do update. otherwise leave it unchanged.
		rand=genrand_real3();
		if(sum_rate!=0){
			index=determine_reaction(sv, rand);
			counterv[index]++;
			for(i=1;i<=NUM_SPEC;i++) {
				S[i]+=stoi[i][index];
			}
		}
		return 1;  //Enzymatic reactions are not finished.


	}
	else if(time_evolution==0){
		while(TIME>INDEX_GRID*GRID_CONST) {
			for(i=1;i<=NUM_SPEC;i++) {
				sumS[i]+=S[i];
				sumS2[i]+=S[i]*S[i];
				if(S[i]<=SMAX[i])	counterS[i][S[i]]++;
				else resize_s_counter(i, S[i]);
			}

			for(i=1;i<=NUM_REACT;i++)
				J[i]= counterv[i]/GRID_CONST;
			for(i=1;i<=NUM_REACT;i++){
				sumJ[i]+=J[i];
				for(j=1;j<=NUM_REACT;j++)
					sumJ2[i][j]+= J[i]*J[j];
				
			}	       
			counter++;
			for(i=1;i<=NUM_REACT;i++) counterv[i]=0;

			Qhead=Qhead%TCORGRID;
			for(i=1;i<=NUM_SPEC;i++)  x_c[i][Qhead] = S[i];
			for(i=1;i<=NUM_SPEC;i++)
				for(j=1;j<=NUM_SPEC;j++)
					cross_pair_sum(Qhead, x_c[i], x_c[j], time_acf[i][j]);
			Qhead++;
//			std::cout<<INDEX_GRID<<"\t"<<counter<<"\t"<<S[1]<<"\t"<<S[2]<<std::endl;

			if(++INDEX_GRID>NGRID_AT_STEADY) return 0;
		}


		//If reaction is not finished, do update. otherwise leave it unchanged.
		rand=genrand_real3();
		if(sum_rate!=0){
			index=determine_reaction(sv, rand);
			counterv[index]++;
			for(i=1;i<=NUM_SPEC;i++) {
				S[i]+=stoi[i][index];
			}
		}
		return 1;  //Enzymatic reactions are not finished.
	}
}


void cross_pair_sum(int Qhead, int  *a, int *b, double *acf)
{
	int t;
	int a0=a[Qhead];
	int i;

	for(t=0;t<=Qhead;t++)           acf[t]+=a0*b[Qhead-t];
	for(i=TCORGRID-1;t<TCORGRID;t++,i--)    acf[t]+=a0*b[i];
}

			


int determine_reaction(double *sv, double rand)
{
	int i;
	double  sum_sv=0;
//	std::cout<<sv[1]<<"\t"<<sv[2]<<"\t"<<sv[3]<<"\t"<<sv[4]<<"\t"<<rand<<"\n"; 
	for(i=1;i<=NUM_REACT;i++){
		sum_sv+=sv[i];
		if(rand<sum_sv) return i;
	}
	printf("Something wrong with SSA code");
	exit(-1);
}



void resize_s_counter(int index, int new_SMAX)
{
	long *m;
	int i,j;
	int old_SMAX;

	old_SMAX=SMAX[index];
	SMAX[index]=new_SMAX;

	m=(long *) calloc( SMAX[index]+10, sizeof(long));
	for(j=0;j<=old_SMAX;j++) m[j]=counterS[index][j];
	//The first counter occurs and thus s_counter is set to 1.
	for(;j<SMAX[index];j++) m[j]=0;
	m[SMAX[index]]=1;
	free(counterS[index]);
	counterS[index]=m;
}


	
void init_state(void)
{
	int i;
	for(i=1;i<=NUM_SPEC;i++) S[i]=INIT_S[i];
	TIME=0;
	DELT=0;
}
	
                                                                                                                                                                                                         








/*
void print_array(int *data, int number)
{
int i;
for (i=0;i<number;i++){
	printf("%d ", data[i]);
	}
printf("\n");
}
*/
