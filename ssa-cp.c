#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "redblack.h"

#include "mioarray.h"
#include "mioregress.h"

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

Array2D *TrajectoryMatrix(Array1D *x, int start, int end, int m, int k)
{
	Array2D *tr;
	int i;
	int j;
	int t_i;
	int t_j;

	int N;

	N = (end - start);
	if(m > N / 2) {
		fprintf(stderr, "m must be less than difference between start and end\n");
		fflush(stderr);
		return(NULL);
	}

	if(m > k) {
		fprintf(stderr,"m must be less than k\n");
		fflush(stderr);
	}

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


double ComputeMu(Array2D *x, 
		 int n, 
		 int lags, 
		 int K, 
		 int Q, 
		 Array1D *ea, 
		 double cv) 
{
	Array2D *tr;
	int i;
	double mu;
	double max_mu = -1000000.0;
	double d;
	double d_tilde;
	double max_d_tilde;

	tr = TrajectoryMatrix(x,0,n,lags,K);
	if(tr == NULL) {
		return(-1);
	}
	d = DStat(tr,ea);
	if(d == -1.0) {
		FreeArray2D(tr);
		return(-1);
	}
	d_tilde = (1.0/(lags*Q)) * d;
	if(d_tilde > max_mu) {
		max_mu = mu;
		max_d_tilde = d_tilde;
	}

	FreeArray2D(tr);

	for(i=1; i < n; i++) {
		if(((n - i)/2) < lags) {
			break;
		}
		tr = TrajectoryMatrix(x,i,n,lags,K);
		if(tr == NULL) {
			return(-1);
		}
		d_tilde = (1.0/(lags*Q)) * d;
		if((d_tilde > -1.0*cv) && (d_tilde < cv)) {
			if(d_tilde > max_mu) {
				max_mu = d_tilde;
			}
		}
		FreeArray2D(tr);
	}

	return(max_mu);
		
}
#ifdef STANDALONE

#define ARGS "x:l:p:q:N:C:"
char *Usage = "usage: ssa-cp -x xfile\n\
\t-C critical value from N(0,1) for hyp. test\n\
\t-l lags\n\
\t-N window size\n\
\t-p starting index of test matrix\n\
\t-q ending index of test matrix\n";

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
	int Q;
	int i;
	int j;
	int N;
	double cv;
	int start;
	double mu;

	lags = 1;
	Q = 0;
	N = 0;
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
				Q = atoi(optarg);
				break;
			case 'N':
				N = atoi(optarg);
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

	if(p == 0) {
		fprintf(stderr,"must specify starting point for test matrix\n");
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


	start = 0;
	if(N > 0) {
		start = p - N;
	} else if(N > p) {
		fprintf(stderr,"p must be >= N\n");
		exit(1);
	} else {
		N = p / 2; // no N => use halfway to starting point
		start = p - N;
		if(start < 0) {
			fprintf(stderr,"start: %d, N: %d, p:%d\n",
					start,N,p);
			exit(1);
		}
	}

	if(Q == 0) {
		Q = p+N;
	}

	if(x->ydim < Q) {
		fprintf(stderr,"data has %d elements, but p+q is %d\n",
			x->ydim,p+Q);
		exit(1);
	}

	base_x = MakeArray1D(N);	// base is everything up to p-1
	if(base_x == NULL) {
		exit(1);
	}

	test_x = MakeArray1D(Q);		// test array is p through q
	if(test_x == NULL) {
		exit(1);
	}


	for(i=0; i < N; i++) {
		base_x->data[i] = x->data[start+i];
	}
	printf("base range: %d %d\n",start,start+N);

	for(i=start+N; i < Q; i++) {
		test_x->data[i-(start+N)] = x->data[i];
	}
	printf("test range: %d %d\n",start+N,Q);
	fflush(stdout);


	/*
	 * for now, do whole matrix
	 */
	tr_base = TrajectoryMatrix(base_x,0,base_x->ydim,lags,base_x->ydim - lags);
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
	tr_test = TrajectoryMatrix(test_x,0,test_x->ydim,lags,Q);
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
	d_base = DStat(tr_base,l_ea) / (lags * Q);
	d_test = DStat(tr_test,l_ea) / (lags * Q);

	mu = ComputeMu(x,p-N,lags,N-lags,Q,l_ea,cv);


	printf("D stats: base: %f test: %f\n",d_base,d_test);
	printf("mu: %f\n",mu);

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
