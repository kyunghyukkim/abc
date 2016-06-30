#define MIN(c1,c2) ((c1<=c2)? c1:c2)

void corfun(int trun, int tcor, int *a, double *acf, int *norm);
void corfun_norm(int tcor, double *acf, int *norm);
void corfun_init(int tcor, double *acf, int *norm);

void corfun_init(int tcor, double *acf, int *norm)
{
	int t;
	for(t=0;t<tcor;t++){
		acf[t] = 0.;
		norm[t]=0;
	}
}
	
void corfun(int trun, int tcor, int *a, double *acf, int *norm)
{
	int t, t0, tt0, tt0max;
	int a0;

	for(t0=0;t0<trun;t0++){
		a0=a[t0];
		tt0max = MIN(trun, t0+tcor);

		for(tt0=t0;tt0<tt0max;tt0++){
			t=tt0-t0;
			acf[t]+=(double) a0*a[tt0];
			norm[t]++;
		}

	}
}

void corfun_norm(int tcor, double *acf, int *norm)
{
	int t;
	for(t=0;t<tcor;t++)
		acf[t]=acf[t]/norm[t];
}

	
