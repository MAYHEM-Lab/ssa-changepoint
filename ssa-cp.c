#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "redblack.h"

#include "mioarray.h"
#include "mioregress.h"

#define BAD_VAL (-1000000.0)

Array2D *LagVarArray2D(Array2D *x)
{
	/*
	 * assumes that x is a trajectory matrix
	 *
	 * each column is a set of rows containing ydim "lags" which correspond to variables
	 * in a multi-variate setting
	 *
	 * the number of columns defines xdim observations of the ydim
	 * dimenstional variable
	 */
	Array2D *lcv;
	Array2D *xt;

	xt = TransposeArray2D(x);
	if(xt == NULL) {
		return(NULL);
	}

	lcv = MultiplyArray2D(x,xt);
	if(lcv == NULL) {
		FreeArray2D(xt);
		return(NULL);
	}

	FreeArray2D(xt);

	return(lcv);
}

Array2D *TrajectoryMatrix(Array1D *x, int start, int m, int k)
{
	Array2D *tr;
	int i;
	int j;
	int t_i;
	int t_j;

	tr = MakeArray2D(m,k);
	if(tr == NULL) {
		return(NULL);
	}

	
	i = start; // start at start
	for(t_i=0; t_i < m; t_i++) {
		j = i; // j walks down the row from current i
		for(t_j=0; t_j < k; t_j++) {
			tr->data[(t_i*tr->xdim)+t_j] = x->data[j];
			j++;
		}
		i++; // there are m different starting points, one for each row
	}

	return(tr);
}

int SortEigenVectors(Array1D *ev, Array2D *ea)
{
	RB *map;
	RB *rb;
	int i;
	int j;
	Array2D *sea;
	Array1D *sev;
	int s_j;

	map = RBInitD();
	if(map == NULL) {
		return(-1);
	}

	sea = MakeArray2D(ea->ydim,ea->xdim);
	if(sea == NULL) {
		RBDestroyD(map);
		return(-1);
	}

	sev = MakeArray1D(ev->ydim);
	if(sev == NULL) {
		RBDestroyD(map);
		FreeArray2D(sea);
		return(-1);
	}

	for(i=0; i < ev->ydim; i++) {
		RBInsertD(map,ev->data[i],(Hval)i);
	}

	s_j = 0;
	/*
	 * move the columns and create sorted ev
	 */
	RB_BACKWARD(map,rb) {
		j = rb->value.i;
		for(i=0; i < ea->ydim; i++) {
			sea->data[i * sea->xdim + s_j] =
				ea->data[i * ea->xdim + j];
		}
		sev->data[s_j] = ev->data[j];
		s_j += 1;
	}

	/*
	 * now overwrite ev with sorted data
	 */
	for(j=0; j < ev->ydim; j++) {
		ev->data[j] = sev->data[j];
	}

	/*
	 * and overwrite ea with sorted ea
	 */
	for(i=0; i < ea->ydim; i++) {
		for(j=0; j < ea->xdim; j++) {
			ea->data[i*ea->xdim+j] = 
				sea->data[i*sea->xdim+j];
		}
	}

	RBDestroyD(map);
	FreeArray2D(sea);
	FreeArray1D(sev);

	return(1);
}
	

double DStat(Array2D *trajectory, Array2D *eigenvectors)
{
	Array2D *ut;
	Array2D *u_ut;
	double d;
	double d_i;
	Array1D *x;
	Array1D *xt;
	Array2D *xt_x;
	Array2D *xt_u_ut;
	Array2D *xt_u_ut_x;
	int i;
	int j;

	ut = TransposeArray2D(eigenvectors);
	if(ut == NULL) {
		return(-1);
	}

	u_ut = MultiplyArray2D(eigenvectors,ut);
	if(u_ut == NULL) {
		FreeArray2D(ut);
		return(-1);
	}

	x = MakeArray1D(trajectory->ydim);
	if(x == NULL) {
		FreeArray2D(ut);
		FreeArray2D(u_ut);
		return(-1);
	}

	d = 0;
	for(j=0; j < trajectory->xdim; j++) {
		for(i=0; i < trajectory->ydim; i++) {
			x->data[i] = trajectory->data[i*trajectory->xdim+j];
		}
		xt = TransposeArray1D(x);
		if(xt == NULL) {
			FreeArray2D(ut);
			FreeArray2D(u_ut);
			FreeArray1D(x);
			return(-1);
		}
		xt_x = MultiplyArray2D(xt,x);
		if(xt_x == NULL) {
			FreeArray2D(ut);
			FreeArray2D(u_ut);
			FreeArray1D(x);
			FreeArray1D(xt);
			return(-1);
		}
		d_i = xt_x->data[0];
		FreeArray2D(xt_x);
		xt_u_ut = MultiplyArray2D(xt,u_ut);
		if(xt_u_ut == NULL) {
			FreeArray2D(ut);
			FreeArray2D(u_ut);
			FreeArray1D(x);
			FreeArray1D(xt);
			return(-1);
		}
		xt_u_ut_x = MultiplyArray2D(xt_u_ut,x);
		if(xt_u_ut_x == NULL) {
			FreeArray2D(ut);
			FreeArray2D(u_ut);
			FreeArray1D(x);
			FreeArray1D(xt);
			FreeArray2D(xt_u_ut);
			return(-1);
		}
		d_i = d_i - xt_u_ut_x->data[0];
		d += d_i;
		FreeArray2D(xt_u_ut_x);
		FreeArray2D(xt_u_ut);
		FreeArray1D(xt);
	}

	return(d);
}


/*
 * x is trajectory matrix for test_matrix
 */
double ComputeMu(Array2D *x, 
		 int lags, 
		 int q,
		 Array1D *ea, 
		 double cv) 
{
	Array2D *tr;
	Array2D *sub;
	int i;
	int j;
	int s;
	int m;
	int j_sub;
	int max_m;
	int max_s;
	double mu;
	double max_mu = BAD_VAL;
	double d;
	double d_tilde;
	double max_d_tilde;

	/*
	 * walk down from m as the number of cols
	 * in test matrix
	 */
	for(m=x->xdim; m > 0; m--) {
		for(s=0; s < x->xdim; s++) {
			if((s+m) <= x->xdim) {
				sub = MakeArray2D(lags,m);
				if(sub == NULL) {
					return(max_mu);
				}
				for(i=0; i < lags; i++) {
					j_sub = 0;
					for(j=s; j < (s+m); j++) {
						sub->data[i*sub->xdim+j_sub] =
							x->data[i*x->xdim+j];
						j_sub++;
					}
				}
				d = DStat(sub,ea);
				if(d == -1.0) {
					fprintf(stderr,
		"ComputeMu failed dstat for m: %d, s: %d\n",
					m,s);
					fflush(stderr);
					FreeArray2D(sub);
					return(max_mu);
				}
				d_tilde = (1.0/(lags*q)) * d;
				if((d_tilde < cv) && 
				   (d_tilde > max_mu)) {
					max_mu = d_tilde;
					max_m = m;
					max_s = s;
					FreeArray2D(sub);
					return(max_mu);
				}

				FreeArray2D(sub);
			}
		}
	}

	return(max_mu);
		
}

int ChangePointIndex(Array2D *x, int lags, int p, int q, Array1D *ea, double mu,
		double cv)
{
	int n;
	int i;
	double h;
	double kappa;
	double W_n;
	double S_n;
	double S_np1;
	Array2D *test_tr;
	Array1D *test_x;

	/*
	 * threshold value
	 *
	 * cv is 1-a upper tail quantile of standard normal
	 */
	h = ((2*cv) / (lags*q)) * sqrt((1.0/3.0) * q * ((3*lags*q) - (q*q) +1));
	kappa = 1.0 / (3 * sqrt(lags * q));
	if(mu == BAD_VAL) {
		return(-1);
	}

	/*
	 * sweep across trajectory matrix 
	 */
	for(n=p+1; n < x->ydim - (lags*q-1); n++) {
		test_x = MakeArray1D(lags*q);
		if(test_x == NULL) {
			exit(1);
		}
		for(i=0; i < lags*q; i++) {
			test_x->data[i] = x->data[i+n];
		}
		test_tr = TrajectoryMatrix(test_x,0,lags,q);
		if(test_tr == NULL) {
			fprintf(stderr,"ChangePointIndex failed to get tr\n");
			fflush(stderr);
			return(-1);
		}
		if(n == p) {
			S_n = (DStat(test_tr,ea) / (lags * q)) / mu;
			W_n = S_n;
		} else {
			S_np1 = (DStat(test_tr,ea) / (lags * q))/ mu;
			W_n = W_n + 
			  S_np1 -
			  S_n -
			  (kappa/sqrt(lags*q));
			if(W_n < 0) {
				W_n = 0;
			}
			S_n = S_np1;
		}
printf("n: %d, h: %f, W_n: %f\n",n,h,W_n);
		FreeArray2D(test_tr);
		FreeArray1D(test_x);
		if(W_n > h) {
			return(n);
		}
	}

	return(-1);
}

int OldChangePointIndex(Array2D *test_x, int lags, int q, Array1D *ea, double mu,
		double cv)
{
	int n;
	double h;
	double kappa;
	double W_n;
	double S_n;
	double S_np1;
	Array2D *test_tr;

	/*
	 * threshold value
	 *
	 * cv is 1-a upper tail quantile of standard normal
	 */
	h = ((2*cv) / (lags*q)) * sqrt((1.0/3.0) * q * ((3*lags*q) - (q*q) +1));
	kappa = 1.0 / (3 * sqrt(lags * q));
	if(mu == BAD_VAL) {
		return(-1);
	}

	/*
	 * sweep across trajectory matrix 
	 */
	for(n=1; n < (lags*q-1); n++) {
		test_tr = TrajectoryMatrix(test_x,n,lags,q);
		if(test_tr == NULL) {
			fprintf(stderr,"ChangePointIndex failed to get tr\n");
			fflush(stderr);
			return(-1);
		}
		if(n == 0) {
			S_n = (DStat(test_tr,ea) / (lags * q)) / mu;
			W_n = S_n;
		} else {
			S_np1 = (DStat(test_tr,ea) / (lags * q))/ mu;
			W_n = W_n + 
			  S_np1 -
			  S_n -
			  (kappa/sqrt(lags*q));
			if(W_n < 0) {
				W_n = 0;
			}
			S_n = S_np1;
		}
printf("n: %d, h: %f, W_n: %f\n",n,h,W_n);
		FreeArray2D(test_tr);
		if(W_n > h) {
			return(n);
		}
	}

	return(-1);
}
#ifdef STANDALONE

#define ARGS "x:l:p:q:K:C:"
char *Usage = "usage: ssa-cp -x xfile\n\
\t-C critical value from N(0,1) for hyp. test\n\
\t-l lags\n\
\t-K number of samples in base matrix\n\
\t-p starting index of test matrix\n\
\t-q number of samples (columns) of of test matrix\n";

char Xfile[4096];

int main(int argc, char *argv[])
{
	int c;
	int size;
	MIO *d_mio;
	MIO *xmio;
	Array1D *ev;
	Array2D *ea;
	Array2D *l_ea;
	Array2D *lcv;
	Array2D *tr_base;
	Array2D *tr_test;
	Array1D *x;
	Array1D *base_x;
	Array1D *test_x;
	int lags;
	double d_base;
	double d_test;
	int p;
	int q;
	int i;
	int j;
	int K;
	int N;
	double cv;
	int start;
	int end;
	double mu;

	p = 0;
	q = 0;
	K = 0;
	cv = 1.96;
	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'l':
				lags = atoi(optarg);
				break;
			case 'p':
				p = atoi(optarg);
				break;
			case 'q':
				q = atoi(optarg);
				break;
			case 'K':
				K = atoi(optarg);
				break;
			case 'C':
				cv = atof(optarg);
				break;
			default:
				fprintf(stderr,
			"unrecognized command: %c\n",(char)c);
				fprintf(stderr,"%s",Usage);
				exit(1);
		}
	}

	if(Xfile[0] == 0) {
		fprintf(stderr,"must specify xfile\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(K == 0) {
		fprintf(stderr,"must enter window size\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(p == 0) {
		fprintf(stderr,"must specify starting point for test matrix\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(q == 0) {
		fprintf(stderr,"must specify sample size for test matrix\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	/*
	 * N is actual window size
	 */
	N = lags * K;

	if(p < N) {
		fprintf(stderr,"window size must have at least p elements ");
		fprintf(stderr,"or base and test matrix overlap\n");
		exit(1);
	}

	if((p - N) <= 0) {
		fprintf(stderr,"p too small for base matrix size N\n");
		exit(1);
	}

	if(lags > N/2) {
		fprintf(stderr,"lags should be less than N / 2\n");
		exit(1);
	}

	size = MIOSize(Xfile);
	d_mio = MIOOpenText(Xfile,"r",size);
	if(d_mio == NULL) {
		fprintf(stderr,"couldn't open %s\n",Xfile);
		exit(1);
	}
	xmio = MIODoubleFromText(d_mio,NULL);
	if(xmio == NULL) {
		fprintf(stderr,"no valid data in %s\n",Xfile);
		exit(1);
	}

	x = MakeArray2DFromMIO(xmio);

	if(x == NULL) {
		exit(1);
	}

	if((p+(lags*q)) > x->ydim) {
		fprintf(stderr,"data has %d elements, but p+(lags*q) is %d\n",
			x->ydim,p+(lags*q));
		exit(1);
	}

	base_x = MakeArray1D(N);	// base is everything up to p-1
	if(base_x == NULL) {
		exit(1);
	}

	test_x = MakeArray1D(lags*q);		// test array is p through q
	if(test_x == NULL) {
		exit(1);
	}


	start = p - N;
	for(i=0; i < N; i++) {
		base_x->data[i] = x->data[start+i];
	}

	end = p + (lags * q);
	for(i=p; i < end; i++) {
		test_x->data[i-p] = x->data[i];
	}


	tr_base = TrajectoryMatrix(base_x,0,lags,K);
	if(tr_base == NULL) {
		fprintf(stderr, 
			"couldn't create trajctory matrix for %d lags\n",
				lags);
		exit(1);
	}
	lcv = LagVarArray2D(tr_base);
	if(lcv == NULL) {
		fprintf(stderr,
			"couldn't compute lag-co-var matrix\n");
		exit(1);
	}
	ev = EigenValueArray2D(lcv);
	if(ev == NULL) {
		fprintf(stderr,"couldn't get eigen values\n");
		exit(1);
	}

	ea = EigenVectorArray2D(lcv);
	if(ea == NULL) {
		fprintf(stderr,"couldn't get eigen vectors\n");
		exit(1);
	}
#if 0
	printf("before sort\n");
	PrintArray1D(ev);
	printf("\n");
	printf("eigen vectors of lag-co-var-matrix\n");
	PrintArray1D(ea);
	printf("\n");
#endif

	SortEigenVectors(ev,ea);

	printf("eigen values of lag-co-var-matrix\n");
	PrintArray1D(ev);
	printf("\n");

	printf("eigen vectors of lag-co-var-matrix\n");
	PrintArray1D(ea);

	/*
	 * for now, drop the last term as being the "noise" series
	 *
	 * ev and ea define the l = M-1 dimensional subspace
	 *
	 * need a way ro choose how many to exclude
	 */

	/*
	 * make trajectory matrix for test
	 */
	tr_test = TrajectoryMatrix(test_x,0,lags,q);
	if(tr_test == NULL) {
		fprintf(stderr,"couldn't mae test trajectory matrix\n");
		exit(1);
	}

	/*
	 * trim rightmost column of eigenvectors
	 */
	l_ea = MakeArray2D(ea->ydim,ea->xdim-1);
	if(l_ea == NULL) {
		exit(1);
	}

	for(i=0; i < l_ea->xdim; i++) {
		for(j=0; j < l_ea->ydim; j++) {
			l_ea->data[i*l_ea->xdim+j] =
				ea->data[i*ea->xdim+j];
		}
	}

	/*
	 * these are d_tilde values due to division
	 */
	d_base = DStat(tr_base,l_ea) / (lags * q);
	d_test = DStat(tr_test,l_ea) / (lags * q);

	mu = ComputeMu(tr_base,lags,K,l_ea,cv);


	printf("D stats: base: %f test: %f\n",d_base,d_test);
	printf("mu: %f\n",mu);

	i = ChangePointIndex(x,lags,p,q,l_ea,mu,cv);
	if(i != -1) {
		printf("change point found at index %d\n",i);
	} else {
		printf("no change point found\n");
	}

	FreeArray2D(lcv);
	FreeArray2D(tr_base);
	FreeArray2D(tr_test);
	FreeArray2D(base_x);
	FreeArray2D(test_x);
	FreeArray2D(x);
	FreeArray1D(ev);
	FreeArray2D(ea);
	FreeArray2D(l_ea);

	return(0);
}

#endif
