#include<iostream>
#include<stdio.h>
#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}



double gammln(double xx);
void nrerror(char error_text[]);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);

void nrerror(char text[])
{
	printf("%s", text);
	exit(-1);
}

void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if(x<=0.0){
		if(x<0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for(n=1;n<=ITMAX;n++){
			++ap;
			del *= x/ap;
			sum += del;
			if(fabs(del) <fabs(sum)*EPS){
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}

void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for(i=1;i<=ITMAX;i++){
		an=-i*(i-a);
		b+=2.0;
		d=an*d+b;
		if(fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if(fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if(fabs(del-1.0) <EPS) break;
	}
	*gammcf=exp(-x+a*log(x) - (*gln))*h;
}
			

double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp-=(x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double square(double a)
{
	return a*a;
}

double gammq(double a, double x)
{
	double gamser, gammcf, gln;

	if(x<0.0 || a<=0.0) nrerror("Invalid arguments in routine gamma");
	if(x< (a+1.0)){
		gser(&gamser, a,x, &gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf, a,x,&gln);
		return gammcf;
	}
}




	
void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a, double *b, double *siga, double *sigb, double *chi2, double *q)
{
	float gammaq(float a, float x);
	int i;
	double wt, t, sxoss,sx=0.0, sy=0.0, st2=0.0,ss,sigdat;
	
	*b=0.0;
	if(mwt) {
		ss=0.0;
		for(i=1;i<=ndata;i++){
			wt=1.0/square(sig[i]);
			ss+=wt;
			sx+=x[i]*wt;
			sy+=y[i]*wt;
		}
	} else {
		for(i=1;i<=ndata;i++){
			sx+=x[i];
			sy+=y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if(mwt){
		for(i=1;i<=ndata;i++){
			t=(x[i]-sxoss)/sig[i];
			st2+=t*t;
			*b+=t*y[i]/sig[i];
		}
	} else {
		for(i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2+=t*t;
			*b+=t*y[i];
		}
	}
	*b/=st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	*q=1.0;
	
	
	if(mwt==0) {
		for(i=1;i<=ndata;i++)	
			*chi2+=square(y[i]-(*a)-(*b)*x[i]);
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *=sigdat;
		*sigb *=sigdat;
	} else{
		for(i=1;i<=ndata;i++)	
			*chi2+=square((y[i]-(*a)-(*b)*x[i])/sig[i]);
		if(ndata>2) *q=gammq(0.5*(ndata-2), 0.5*(*chi2));
	}
	
}

/*
void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[], int ma, double **covar, double *chisq, double (*funcs)(double, double [], int))
{
	void covsrt(double **covar, int ma, int ia[], int mfit);
	void gaussj(double **a, int n, double **b, int m);
	int i,j,k,l,m,mfit=0;
	double ym,wt,sum,sig2i,**beta,*afunc;

	beta=matrix(1,ma,1,1);
	afunc=vector(1,ma);
	for(j=1;j<=ma;j++)
		if(ia[j]) mfit++;
	if(mfit==0) nrerror("lfit: no parameters to be fitted");
	for(j=1;j<=mfit;j++){
		for(k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for(i=1;i<=ndat;i++){
		(*funcs)(x[i], afunc,ma);
		ym=y[i];
		if(mfit<ma){
			for(j=1;j<=ma;j++)
				if(!ia[j]) ym-=a[j]*afunc[j];
		}
		sig2i=1.0/SQR(sig[i]);
		for(j=0,l=1;l<=ma;l++){
			if(ia[l]){
				wt=afunc[l]*sig2i;
				for(j++,k=0,m=1;m<=l;m++)
					if(ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1]+=ym*wt;
			}
		}
	}
	for(j=2;j<=mfit;j++)
		for(k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	gaussj(covar,mfit,beta,1);
	for(j=0,l=1;l<=ma;l++)
		if(ia[l]) a[l]=beta[++j][1];
	*chisq=0.0;
	for(i=1;i<=dat;i++){
		(*funcs)(x[i],afunc,ma);
		for(sum=0.0,j=1;j<=ma;j++) sum+= a[j]*afunc[j];
		*chisq+= SQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,ia,mfit);
	free_vector(afunc,1,ma);
	free_matrix(beta,1,ma,1,1);
}

void covsrt(double **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	double swap;

	for(i=mfit+1;i<=ma;i++)
		for(j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for(j=ma;j>=1;j--){
		if(ia[j]){
			for(i=1;i<=ma;i++) SWAP(covar[i][k], covar[i][j]);
			for(i=1;i<=ma;i++) SWAP(covar[k][i], covar[j][i]);
			k--;
		}
	}
}
			
			
*/
