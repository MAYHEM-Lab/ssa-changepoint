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

	sea = MakeArray2D(ea->xdim,ea->ydim);
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
	for(i=0; i < trajectory->ydim; i++) {
		for(j=0; j < trajectory->xdim; j++) {
			x->data[j] = trajectory->data[i*trajectory->xdim+j];
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
		

#ifdef STANDALONE

#define ARGS "x:l:"
char *Usage = "usage: ssa-cp -x xfile\n\
\t-l lags\n";

char Xfile[4096];

int main(int argc, char *argv[])
{
	int c;
	int size;
	MIO *d_mio;
	MIO *xmio;
	Array1D *ev;
	Array2D *ea;
	Array2D *lcv;
	Array2D *tr;
	Array1D *x;
	int lags;
	double d;

	lags = 1;
	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'l':
				lags = atoi(optarg);
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

	/*
	 * for now, do whole matrix
	 */
	tr = TrajectoryMatrix(x,0,x->ydim,lags,x->ydim - lags);
	if(tr == NULL) {
		fprintf(stderr, 
			"couldn't create trajctory matrix for %d lags\n",
				lags);
		exit(1);
	}
	lcv = LagVarArray2D(tr);
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
	printf("before sort\n");
	PrintArray1D(ev);
	printf("\n");
	printf("eigen vectors of lag-co-var-matrix\n");
	PrintArray1D(ea);
	printf("\n");



	SortEigenVectors(ev,ea);

	printf("eigen values of lag-co-var-matrix\n");
	PrintArray1D(ev);
	printf("\n");

	printf("eigen vectors of lag-co-var-matrix\n");
	PrintArray1D(ea);

	d = DStat(tr,ea);

	printf("D stat: %f\n",d);

	FreeArray2D(lcv);
	FreeArray2D(tr);
	FreeArray1D(ev);
	FreeArray2D(ea);

	return(0);
}

#endif
