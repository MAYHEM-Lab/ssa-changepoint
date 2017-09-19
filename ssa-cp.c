#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

#ifdef STANDALONE

#define ARGS "x"
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
	tr = TrajectoryMatrix(x,0,x->xdim,lags,x->xdim - lags);
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
	printf("eigen values of lag-co-var-matrix\n");
	PrintArray1D(ev);
	printf("\n");

	ea = EigenVectorArray2D(lcv);
	if(ea == NULL) {
		fprintf(stderr,"couldn't get eigen vectors\n");
		exit(1);
	}
	printf("eigen vectors of lag-co-var-matrix\n");
	PrintArray1D(ea);

	FreeArray2D(lcv);
	FreeArray2D(tr);
	FreeArray1D(ev);
	FreeArray2D(ea);

	return(0);
}

#endif
