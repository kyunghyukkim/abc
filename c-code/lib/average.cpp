#include <math.h>

//block average: number of block is N, the size of each block is M
//array a[i] with i=0,1,2,....
void average(int N, int M, double *a, double *a_sum)
{
	int i,j;
	
	for(i=0;i<N;i++) *(a_sum+i) = 0.;

	for(i=0;i<N;i++){
	       for(j=0;j<M;j++)	{
		       *(a_sum+i)+=*(a+j*N + i);
	       }
	       *(a_sum+i) = (double) *(a_sum+i)/M;
	}
}
 
//standard deviation and average of a data set is calculated. 
//array a[i] with i=0,1,2,....
void ave_error(int M, double *a, double *a_ave, double *a_error)
{
	int i;
	double ai;
	
	*a_ave=0.;
	*a_error=0.;
	for(i=0;i<M;i++)	{
		ai=*(a+i);
	       (*a_ave)+=ai;
	       (*a_error)+=ai*ai;
	}
	*a_ave = (double) *a_ave/ M;
	
	double temp; 
	temp = (double)(*a_error/M)-*a_ave* *a_ave; 
	if (temp >=0) *a_error= sqrt(temp);
	else *a_error = 0;

}


//standard deviation and average of a data set is calculated. 
//array a[i] with i=1,2,....,M
void ave_error1(int M, double *a, double *a_ave, double *a_error)
{
	int i;
	double ai;
	
	*a_ave=0.;
	*a_error=0.;
	for(i=1;i<=M;i++)	{
		ai=*(a+i);
	       (*a_ave)+=ai;
	       (*a_error)+=ai*ai;
	}
	*a_ave = (double) *a_ave/M;
	
	double temp; 
	temp = (double)(*a_error/M)-*a_ave* *a_ave; 
	if (temp >=0) *a_error= sqrt(temp);
	else *a_error = 0;

}
//standard deviation and average of a data set is calculated. 
//array a[i] with i=0,1,2,....
//array a[i] is integer one.
void ave_error(int M, int *a, double *a_ave, double *a_error)
{
	int i;
	double ai;
	
	*a_ave=0.;
	*a_error=0.;
	for(i=0;i<M;i++)	{
		ai=(double)*(a+i);
	       (*a_ave)+=ai;
	       (*a_error)+=ai*ai;
	}
	*a_ave = (double) *a_ave/M;
	double temp; 
	temp = (double)(*a_error/M)-*a_ave* *a_ave; 
	if (temp >=0) *a_error= sqrt(temp);
	else *a_error = 0;

}
//standard deviation and average of a data set with weight is calculated. 
//array a[i] with i=0,1,2,....
//array a[i] is double-type.
void ave_error(int M, double *a, double *w, double *a_ave, double *a_error)
{
	int i;
	double ai;
	double wi; 
	double sumwi;

	*a_ave=0.;
	*a_error=0.;
	sumwi=0;

	for(i=0;i<M;i++)	{
		ai=*(a+i);
		wi=*(w+i);
		sumwi+=wi;
	       (*a_ave)+=ai*wi;
	       (*a_error)+=ai*ai*wi;
	}
	*a_ave = (double) *a_ave/sumwi;

	double temp; 
	temp = (double)(*a_error/sumwi)-*a_ave* *a_ave; 
	if (temp >=0) *a_error= sqrt(temp);
	else *a_error = 0;
}

