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
//	sea = NormalizeColsArray2D(ea);

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

	FreeArray2D(ut);
	FreeArray2D(u_ut);
	FreeArray2D(x);
	return(d);
}


/*
 * x is data matrix for test_matrix
 */
double NewComputeMu(Array2D *x, 
		 int start,
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
	int sub_i;

	sub = MakeArray2D(lags,q);
	if(sub == NULL) {
		exit(1);
	}

	/*
	 * walk down from m as the number of cols
	 * in test matrix
	 */
	for(m=start; m > lags+q; m--) {
		sub_i=0;
		for(i=m-(q+lags); i < start; i++) {
			sub->data[sub_i] = x->data[i];
			sub_i++;
		}
		tr = TrajectoryMatrix(sub,0,lags,q);
		if(tr == NULL) {
			exit(0);
		}
		d = DStat(tr,ea);
		if(d == -1.0) {
			fprintf(stderr,
	"ComputeMu failed dstat for m: %d, s: %d\n",
			m,s);
			fflush(stderr);
			FreeArray2D(sub);
			return(max_mu);
		}
		d_tilde = (1.0/(lags*q)) * d;
#if 0
		if((d_tilde < cv) && 
		   (d_tilde > (-1.0*cv)) &&
		   (d_tilde > max_mu)) {
#endif
		if((d_tilde < cv) &&
		   (d_tilde > max_mu)) {
			max_mu = d_tilde;
			max_m = m;
			max_s = s;
			FreeArray2D(sub);
			FreeArray2D(tr);
			return(max_mu);
		}

		FreeArray2D(tr);
	}

	FreeArray2D(sub);
	fprintf(stderr,"no mu\n");
	fflush(stderr);
	return(max_mu);
		
}

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
				   (d_tilde > (-1.0*cv)) &&
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

	fprintf(stderr,"ComputeMu: no mu\n");
	return(max_mu);
		
}

int ChangePointSweep(Array2D *x, int lags, int K, 
		int q, double cv, int signal)
{
	Array1D *ev;
	Array2D *ea;
	Array2D *l_ea;
	Array2D *tea;
	Array2D *lcv;
	Array2D *tr_base;
	Array2D *tr_test;
	Array2D *tr_before;
	Array1D *base_x;
	Array1D *test_x;
	Array1D *before_x;
	int start;
	int end;
	int p;
	int i;
	int j;
	int x_n;
	int N;
	double mu;
	double mu_np1;
	double h;
	double kappa;
	double W_n;
	double W_np1;
	double S_n;
	double S_np1;
	double d;

	N = lags + K;
	if(x->ydim - (N+(lags+q)) <= 0) {
		fprintf(stderr,
			"data has %d values but requires %d for sweep\n",
			x->ydim,N+(lags*q));
		fflush(stderr);
		return(-1);
	}

	base_x = MakeArray1D(N);	// base is everything up to p-1
	if(base_x == NULL) {
		exit(1);
	}

	test_x = MakeArray1D(lags+q);		// test array is p through q
	if(test_x == NULL) {
		exit(1);
	}

	h = ((2*cv) / (lags*q)) * sqrt((1.0/3.0) * q * ((3*lags*q) - (q*q) +1));
	kappa = 1.0 / (3 * sqrt(lags * q));


	for(start = 0; start < (x->ydim - (N+(lags+q))); start++) {

		before_x = MakeArray1D(start+N+lags+q);
		if(before_x == NULL) {
			exit(1);
		}

		for(i=0; i < start+N+lags+q; i++) {
			before_x->data[i] = x->data[i];
		}
		for(i=0; i < N; i++) {
			base_x->data[i] = x->data[start + i];
		}

		p = start + N;	// p starts immediately after the base array
		for(i=0; i < lags+q; i++) {
			test_x->data[i] = x->data[p+i];
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
		SortEigenVectors(ev,ea);

		printf("eigen values of lag-co-var-matrix\n");
		PrintArray1D(ev);
		printf("\n");

		printf("eigen vectors of lag-co-var-matrix\n");
		PrintArray1D(ea);

		/*
		 * for now, drop the last term as being the "signal" series
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
		l_ea = MakeArray2D(ea->ydim,signal);
		if(l_ea == NULL) {
			exit(1);
		}
		for(i=0; i < l_ea->xdim; i++) {
			for(j=0; j < l_ea->ydim; j++) {
				l_ea->data[i*l_ea->xdim+j] =
					ea->data[i*ea->xdim+j];
			}
		}

		tr_before = TrajectoryMatrix(before_x,0,lags,K);
		if(tr_before == NULL) {
			exit(1);
		}


		if(start == N) { // we are computing S_1
//			mu = ComputeMu(tr_test,lags,K,l_ea,cv);
			mu = ComputeMu(tr_base,lags,K,l_ea,cv);
//			mu = ComputeMu(tr_before,lags,K,l_ea,cv);
//			mu = ComputeMu(x,start+N,lags,K,l_ea,cv);
			S_n = (DStat(tr_test,l_ea)/(lags*q))/mu;
			W_n = S_n;
			d = DStat(tr_test,l_ea) / (lags*q);
		} else { // we have S_N from previous iteration
//			mu_np1 = ComputeMu(tr_test,lags,K,l_ea,cv);
			mu_np1 = ComputeMu(tr_base,lags,K,l_ea,cv);
//			mu_np1 = ComputeMu(tr_before,lags,K,l_ea,cv);
//			mu_np1 = ComputeMu(x,start+N,lags,K,l_ea,cv);
			S_np1 = (DStat(tr_test,l_ea)/(lags*q))/mu_np1;
			W_np1 = W_n +
				S_np1 -
				S_n -
				(kappa/sqrt(lags*q));
			if(W_np1 < 0) {
				W_np1 = 0;
			}
printf("S_n: %f S_np1: %f mu: %f m_np1: %f kappa: %f\n",S_n,S_np1,mu,mu_np1,kappa);
fflush(stdout);

			S_n = S_np1; // for next iteration
			W_n = W_np1;
			mu = mu_np1;
			d = DStat(tr_test,l_ea) / (lags*q);
		}


		FreeArray2D(lcv);
		FreeArray2D(tr_base);
		FreeArray2D(tr_test);
		FreeArray1D(ev);
		FreeArray2D(ea);
		FreeArray2D(l_ea);
		FreeArray1D(before_x);
		FreeArray2D(tr_before);
printf("start: %d, target: %d, h: %f, W_n: %f, d: %f mu: %f\n",start,start+N,h,W_n,d,mu);
fflush(stdout);

	}

	FreeArray2D(base_x);
	FreeArray2D(test_x);

	return(-1);
}
#ifdef STANDALONE

#define ARGS "x:l:p:q:K:C:e:"
char *Usage = "usage: ssa-cp -x xfile\n\
\t-e number of signal series\n\
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
	Array1D *x;
	Array1D *x_sub;
	int lags;
	int p;
	int q;
	int i;
	int K;
	int N;
	double cv;
	int signal;
	

	p = 0;
	q = 0;
	K = 0;
	cv = 1.96;
	signal = 1;
	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'e':
				signal = atoi(optarg);
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

	if(q == 0) {
		fprintf(stderr,"must specify sample size for test matrix\n");
		fprintf(stderr,"%s",Usage);
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

	if(p == 0) { // if p == 0, start from beginning of x
		i = ChangePointSweep(x,lags,K,q,cv,signal);
	} else {
		x_sub = MakeArray1D(x->ydim - p);
		if(x_sub == NULL) {
			exit(1);
		}
		for(i=p; i < x->ydim; i++) {
			x_sub->data[i-p] = x->data[i];
		}
		i = ChangePointSweep(x_sub,lags,K,q,cv,signal);
		FreeArray1D(x_sub);
		if(i != -1) {
			i = i + p; // quote i from original series
		}
	}

	if(i == -1) {
		printf("no change point found\n");
	} else {
		printf("change point found at %d\n",i);
	}

	return(0);
}

#endif
