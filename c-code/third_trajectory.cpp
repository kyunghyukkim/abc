#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cstring>
#include "./lib/mt19937ar-cok.cpp"
#include "./lib/corfun.cpp"
#include "./lib/average.cpp"

#include<iostream>
using namespace std;

//using namespace std;
#define MAX(c1,c2) ((c1>c2)? c1:c2)


extern "C" {
	char * gillespie(char*);
}

long  *S, *INIT_S;
double  *C;


double TIME;
double DELT; //time increment in each MC step.
double  *aveS, *stdS;
double  *sumS, *sumS2;
long **counterS;
double *v;

int **stoi;
int N;

int NUM_SPEC;
int NUM_SPEC2;
int NUM_PARAM;
int NUM_REACT;

double GRID_CONST;
long NGRID_BF_STEADY;
int NENSEMBLE;

long INDEX_GRID;
long *SMAX;

void model_def(char *, std::stringstream& str);
void init_state(void);
void time_evol(std::stringstream& str);
void resize_s_counter(int index_species, int new_copy_numb_of_moleculesi);
void propensity(void);
int determine_reaction(double *sv, double rand);
// void create_output_files(void);
bool  MC_not_finished(std::stringstream& str);

// std::ofstream FILE_TEVOL;
// std::ofstream FILE_FLUX;






char * gillespie(char *str)
{
	cout << "str\n";
   cout << str;
   cout << "\n";
   cout << "*str\n";
   cout << *str;
   cout << "\n";
   
   std::stringstream STR_OUTPUT;
	
   cout << "STR_OUTPUT\n";
   cout << STR_OUTPUT;
   cout << "\n";
   cout << "str before\n";
   cout << str;
   cout << "\n";
   
	model_def(str, STR_OUTPUT);
	time_evol(STR_OUTPUT);
	
	const std::string& tmp = STR_OUTPUT.str();
   cout << "tmp\n";
   cout << tmp;
   cout << "\n";
   
	STR_OUTPUT.str("");
	STR_OUTPUT.clear();
	
	const char* cstr = tmp.c_str();
   
   cout << "cstr\n";
   cout << cstr;
   cout << "\n";

	int i=0;
	char * cstr2 = static_cast<char *>(malloc(strlen(cstr) * sizeof(char)));
	while (cstr[i]){
		cstr2[i] = cstr[i];
		i++;
		}
   cout << "Hello You!\n";
   cout << cstr2;
   cout << "\n";
	return cstr2;
	//return cstr;
}

	
void model_def(char *str_c, std::stringstream& STR_OUTPUT)
{
	int i,j ;
	
   cout << "str_c\n";
   cout << str_c;
   cout << "\n";
   
	std::stringstream str(str_c);

   cout << "str STREAMED\n";
   cout << str.str();
   cout << "\n";
   
   cout << "DOUBLE DOWN\n";
   
   //double Test = 2.2;
   
   cout << "DOUBLE UP\n";
   
   //cout << Test;
   //cout << "\n";

	NUM_SPEC=3;
	NUM_SPEC2=NUM_SPEC*NUM_SPEC;
	NUM_PARAM=4;
	NUM_REACT=4;

	S=(long *) calloc(NUM_SPEC+1, sizeof(long));
	INIT_S=(long *)calloc(NUM_SPEC+1, sizeof(long));
	aveS=(double *)calloc(NUM_SPEC+1, sizeof(double));
	stdS=(double *)calloc(NUM_SPEC+1, sizeof(double));
	sumS=(double *)calloc(NUM_SPEC+1, sizeof(double));
	sumS2=(double *)calloc(NUM_SPEC+1, sizeof(double));
   C=(double *)calloc(NUM_PARAM+1, sizeof(double));
	counterS=(long **)calloc(NUM_SPEC+1, sizeof(long *));
	SMAX=(long *)calloc(NUM_SPEC+1, sizeof(long));
   
   //cout << "Expect Failure\n";
   //cout << aveS[1];
   //cout << "\n";
   
   cout << "INIT_S\n";
   cout << INIT_S;
   cout << "\n";
   
   cout << "C\n";
   cout << C;
   cout << "\n";

   cout << "STREAMING IN 1\n";
	for(i=1;i<=NUM_SPEC;i++) {
		str >> INIT_S[i];
      cout << "INIT_S burh\n";
      cout << INIT_S[i];
      cout << "\n";
	}

   cout << "STREAMING IN 2\n";
	for(i=1;i<=NUM_PARAM;i++) {
		str >> C[i];
      //cout << "C burh\n";
      //cout << C[i];
      //cout << "\n";
	}

	 str >> GRID_CONST >> NGRID_BF_STEADY ;
    
    //cout << "\nGRID_CONST\n";
    //cout << GRID_CONST;
    cout << "\nNGRID_BF_STEADY\n";
    cout << NGRID_BF_STEADY;
    cout << "\n";
    
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
   
   cout << "MODEL\n";



}

void propensity(void)
{
	v[1]=C[1];
	v[2]=C[2]*S[1];
	v[3]=C[3]*S[1]*S[2];
	v[4]=C[4]*S[3];
}



void time_evol(std::stringstream& STR_OUTPUT)
{
	
   cout << "TIME START\n";
   
	init_genrand(time(NULL)+genrand_int32());
   
   cout << "Geriatric\n";
   
	init_state();
   
   cout << "State\n";		

	//STR_OUTPUT << std::endl;	
	//STR_OUTPUT << "Time" << " ";
	//for (i=1;i<=NUM_SPEC;i++){
	//	STR_OUTPUT << "S" <<i<< " ";
	//}
	//STR_OUTPUT << std::endl;	
	while(INDEX_GRID < NGRID_BF_STEADY )  MC_not_finished(STR_OUTPUT);
   
   cout << "TIME\n";
}

bool  MC_not_finished(std::stringstream& STR_OUTPUT)
{
	double sv[NUM_REACT];
	double sum_rate=0;
	double rand;
	int index;
	int i;

   cout << "Going to Propensity\n";

	propensity();
	for(i=1;i<=NUM_REACT;i++) sum_rate+=v[i];

   cout << "Propensity For-d\n";

	if(sum_rate!=0) {
		DELT=-log(genrand_real3())/sum_rate;
		for(i=1;i<=NUM_REACT;i++) sv[i]=v[i]/sum_rate;
	}
   
   cout << "IF?\n";
   
	//if reaction finishes, then new DELT is going to be the last DELT just before reaction finishes.
	TIME+=DELT;
   cout << "YOU\n";
	while(TIME>INDEX_GRID*GRID_CONST) {
      cout << "SAY\n";
      cout << "INDEX GRID\n";
      cout << INDEX_GRID;
      cout << "\n";
      //cout << "GRID_CONST\n";
      //cout << GRID_CONST;
      //cout << "\n";
		
      cout << "Can I do numbers?\n";
      long INBETWEEN = INDEX_GRID*GRID_CONST;
      double DOUBLEINBETWEEN = INDEX_GRID*GRID_CONST;
      cout << "NOW THE SIN\n";
      cout << INBETWEEN;
      cout << "\n";
      
      cout << "Can I do Strings\n";
      cout << "INDEX_GRID\n";
      cout << INDEX_GRID << " ";
      cout << "\n";
      //cout << "GRID_CONST\n";
      //cout << GRID_CONST << " ";
      //cout << "\n";
      
      cout << "Can I do String Stream?\n";
      cout << "INDEX_GRID\n";
      STR_OUTPUT << INDEX_GRID << " ";
      cout << "\n";
      cout << "GRID_CONST\n";
      STR_OUTPUT << GRID_CONST;
      cout << "\n";

      
      //STR_OUTPUT <<  INDEX_GRID*GRID_CONST << " ";
      cout << "SO\n";
		for(i=1;i<=NUM_SPEC;i++){
			STR_OUTPUT << S[i] << " ";
		}
		cout << "FUCK\n";
		STR_OUTPUT << std::endl;
		cout << S[1];
		if(++INDEX_GRID>NGRID_BF_STEADY) return 0;
	}

	//If reaction is not finished, do update. otherwise leave it unchanged.
	rand=genrand_real3();
	if(sum_rate!=0){
		index=determine_reaction(sv, rand);
		for(i=1;i<=NUM_SPEC;i++) {
			S[i]+=stoi[i][index];
		}
	}
	return 1;  //Enzymatic reactions are not finished.
}



int determine_reaction(double *sv, double rand)
{
	int i;
	double  sum_sv=0;
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
	int j;
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
	INDEX_GRID=0;
}

