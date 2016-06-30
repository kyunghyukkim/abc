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


extern "C" {
	const char* gillespie(char*);
}

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

long INDEX_GRID;
long *SMAX;
std::stringstream STR_OUTPUT;

void model_def(char *);
void init_state(void);
void time_evol(void);
void resize_s_counter(int index_species, int new_copy_numb_of_moleculesi);
void propensity(void);
int determine_reaction(double *sv, double rand);
void create_output_files(void);
bool  MC_not_finished(int time_evolution);

std::ofstream FILE_TEVOL;
std::ofstream FILE_FLUX;






const char* gillespie(char *str)
{
	long nextseed=0;
		
	model_def(str);
//	create_output_files();
	
	time_evol();
	
	const std::string& tmp = STR_OUTPUT.str();
	const char* cstr = tmp.c_str();
	return cstr;
//	time_ave();
}

	
void model_def(char *str_c)
{
	int i,j ;
	
	std::stringstream str(str_c);

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


	int d1[10];
	for(i=1;i<=NUM_SPEC;i++) {
		str >> INIT_S[i];
		std::cout << INIT_S[i] << std::endl;
	}

	for(i=1;i<=NUM_PARAM;i++) {
		str >> C[i];
		std::cout << C[i] << std::endl;
	}

	str >> GRID_CONST >> NGRID_BF_STEADY >> NGRID_AT_STEADY;
	std::cout << GRID_CONST << " " << NGRID_BF_STEADY<< " " << NGRID_AT_STEADY << std::endl;
	std::cout <<  "NGRID_AT_STEADY = " <<NGRID_AT_STEADY << std::endl;
	// output string memory allocation
	// time and S1 with NGRID_BF_STEADY items. 
	
	for(i=1;i<=NUM_SPEC;i++){
		SMAX[i]=INIT_S[i];
		counterS[i]=(long *) calloc( SMAX[i]+1, sizeof(long));
		if(SMAX[i]<S[i]) resize_s_counter(i, S[i]);
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
		

	for(i=1;i<=NUM_REACT;i++) counterv[i]=0;
	while(INDEX_GRID < NGRID_BF_STEADY )  MC_not_finished(1);
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
			STR_OUTPUT <<  INDEX_GRID*GRID_CONST << " ";
			for(i=1;i<=NUM_SPEC;i++){
				STR_OUTPUT << S[i] << " ";
			}
			STR_OUTPUT << std::endl;
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
	
                                                                                                                                                                                                         






