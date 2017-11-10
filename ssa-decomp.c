#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>

#include "redblack.h"
#include "dlist.h"

#include "mioarray.h"
#include "mioregress.h"

#define RAND() (drand48())

#define BAD_VAL (-1000000.0)

struct object_stc
{
	int id;
	Array1D *series;
	int L;
	int K;
	char *str;
};

typedef struct object_stc Object;

struct bin_stc
{
	Array1D *centroid;
	Dlist *list;
	int count;
};

typedef struct bin_stc Bin;


Array2D *UnitizeArray2D(Array2D *u)
{
	Array2D *on_u;

	on_u = NormalizeColsArray2D(u);
	return(on_u);
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



Array2D *DiagonalAverage(Array2D *X)
{
	Array2D *Y;
	int i;
	int j;
	int k;
	int N;
	double acc;
	double tot;

	/*
	 * first row and last column
	 */
	N = X->xdim + (X->ydim - 1);
	Y = MakeArray1D(N);
	if(Y == NULL) {
		exit(1);
	}

	for(k = 0; k < N ; k++) {
		acc = 0;
		tot = 0;
		for(i=0; i < X->ydim; i++) {
			for(j=0; j < X->xdim; j++) {
				if((i + j) == k) {
					acc += X->data[i*X->xdim+j];
					tot++;
				}
			}
		}
		Y->data[k] = acc / tot;
	}

	return(Y);
}


Array1D *SSADecomposition(Array1D *series, int L, int start, int end)
{
	Array1D *ev;
	Array2D *ea;
	Array2D *U;
	Array2D *V;
	Array2D *S;
	Array2D *tr;
	Array2D *D;
	Array2D *D2;
	Array2D *Dm2;
	Array2D *X;
	Array2D *X_t;
	Array2D *V_t;
	Array2D *rank_1;
	Array2D *t_u;
	Array2D *t_x;
	Array2D *t_y;
	Array2D *diff;
	Array2D *sum;
	int K;
	int N;
	Array1D *U_col;
	Array1D *V_row;
	int d;
	int i;
	int j;
	int k;
	Array1D **dcomp;
	Array1D *Y;
	int y_i;
	Array2D *I;
	double acc;


	N = series->ydim;
	K = N - L + 1;


	X = TrajectoryMatrix(series,0,L,K);
	if(X == NULL) {
		fprintf(stderr, 
			"couldn't create trajctory matrix for %d lags\n",
				L);
		exit(1);
	}
	X_t = TransposeArray2D(X);
	if(X_t == NULL) {
		fprintf(stderr,"couldn't transpose array X\n");
		exit(1);
	}

	S = MultiplyArray2D(X_t,X);
	if(S == NULL) {
		fprintf(stderr,
			"couldn't compute S\n");
		exit(1);
	}
	FreeArray2D(X_t);

	ev = EigenValueArray2D(S);
	if(ev == NULL) {
		fprintf(stderr,"couldn't get eigen values\n");
		exit(1);
	}

	ea = EigenVectorArray2D(S);
	if(ea == NULL) {
		fprintf(stderr,"couldn't get eigen vectors\n");
		exit(1);
	}
	FreeArray2D(S);

	V = UnitizeArray2D(ea);
	if(V == NULL) {
		fprintf(stderr,"couldn't unitize eigen vectors\n");
		exit(1);
	}
	FreeArray2D(ea);
	SortEigenVectors(ev,V);


	t_u = MultiplyArray2D(X,V);
	if(t_u == NULL) {
		exit(1);
	}

	/*
	 * need number of cols from V_t
	 */
	Dm2 = MakeArray2D(t_u->xdim,V->xdim);
	if(Dm2 == NULL) {
		exit(1);
	}

	D2 = MakeArray2D(t_u->xdim,V->xdim);
	if(D2 == NULL) {
		exit(1);
	}

	D = MakeArray2D(t_u->xdim,V->xdim);
	if(D == NULL) {
		exit(1);
	}

	for(i=0; i < Dm2->ydim; i++) {
		for(j=0; j < Dm2->xdim; j++) {
			Dm2->data[i*D->xdim+j] = 0;
			D2->data[i*D->xdim+j] = 0;
			D->data[i*D->xdim+j] = 0;
		}
	}

	/*
	 * diagonal is 1/sqrt(eigenvalue) for non-zero entries
	 *
	 * if L < K, leave others zero
	 */
	for(i=0; i < L; i++) {
		if(i < Dm2->xdim) {
			Dm2->data[i*Dm2->xdim+i] = 1.0/sqrt(ev->data[i]);
			D2->data[i*D2->xdim+i] = sqrt(ev->data[i]);
			D->data[i*D->xdim+i] = ev->data[i];
		}
	}

	U = MultiplyArray2D(t_u,Dm2);
	if(U == NULL) {
		fprintf(stderr,"can't get U\n");
		exit(1);
	}
	FreeArray2D(t_u);

#if 0
printf("L: %d, K: %d\n",L,K);
t_u = TransposeArray2D(U);
I = MultiplyArray2D(t_u,U);
PrintArray2D(I);
exit(1);


for(i=0; i < I->ydim; i++) {
	for(j=0; j < I->xdim; j++) {
		if(fabs(I->data[i*I->xdim+j]) < (1.0/1000000.0)) {
			I->data[i*I->xdim+j] = 0.0;
		}
		printf("%1.0f ",I->data[i*I->xdim+j]);
	}
	printf("\n");
}
t_u = TransposeArray2D(U);
I = MultiplyArray2D(t_u,U);
printf("U is %d x %d\n",U->ydim,U->xdim);
printf("U_t * U\n");
for(i=0; i < I->ydim; i++) {
	for(j=0; j < I->xdim; j++) {
		if(fabs(I->data[i*I->xdim+j]) < (1.0/1000000.0)) {
			I->data[i*I->xdim+j] = 0.0;
		}
		printf("%1.1f ",I->data[i*I->xdim+j]);
	}
	printf("\n");
}

printf("D\n");
PrintArray2D(D);
t_u = TransposeArray2D(D);
printf("D^t\n");
PrintArray2D(t_u);
exit(1);
for(i=0; i < K; i++) {
	printf("eigenvalue_%d: %f\n",i+1,ev->data[i]);
}
printf("VT\n");
PrintArray2D(V_t);
V_t = InvertArray2D(V);
printf("V-1\n");
PrintArray2D(V_t);

exit(1);

I = MultiplyArray2D(V,V_t);
PrintArray2D(I);
exit(1);
t_u = MultiplyArray2D(V,D);
t_x = MultiplyArray2D(t_u,V_t);
printf("S\n");
PrintArray2D(S);
printf("XT * X\n");
PrintArray2D(t_x);
exit(1);
for(i=0; i < U->ydim; i++) {
	printf("col0: %f\n",U->data[i*U->xdim+0]);
}
for(i=0; i < U->ydim; i++) {
	printf("col1: %f\n",U->data[i*U->xdim+1]);
}
for(i=0; i < U->ydim; i++) {
	printf("col2: %f\n",U->data[i*U->xdim+2]);
}
exit(1);
#endif
	V_t = TransposeArray2D(V);
//	V_t = InvertArray2D(V);
	if(V_t == NULL) {
		fprintf(stderr,"couldn't get V_t\n");
		exit(1);
	}
	FreeArray2D(V);
	FreeArray2D(X);


	U_col = MakeArray1D(U->ydim);
	if(U_col == NULL) {
		exit(1);
	}

	V_row = MakeArray2D(1,V_t->xdim);
	if(V_row == NULL) {
		exit(1);
	
	}

	sum = MakeArray2D(U_col->ydim,V_row->xdim);
	if(sum == NULL) {
		exit(1);
	}

	for(i=0; i < sum->ydim; i++) {
		for(j=0; j < sum->xdim; j++) {
			sum->data[i*sum->xdim+j] = 0;
		}
	}

	for(j = start; j < end; j++) {
		for(i=0; i < U_col->ydim; i++) {
			U_col->data[i] = U->data[i*U->xdim+j] * sqrt(ev->data[j]);
		}
		for(i=0; i < V_row->xdim; i++) {
			V_row->data[0*V_row->xdim + i] =
				V_t->data[j*V_t->xdim+i];
		}

		rank_1 = MultiplyArray2D(U_col,V_row);
		if(rank_1 == NULL) {
			exit(1);
		}

		for(i=0; i < sum->ydim; i++) {
			for(k=0; k < sum->xdim; k++) {
				sum->data[i*sum->xdim+k] =
					sum->data[i*sum->xdim+k] +
				        rank_1->data[i*rank_1->xdim+k];
			}
		}
#if 0
acc = 0;
for(i=0; i < rank_1->ydim; i++) {
	for(k=0; k < rank_1->xdim; k++) {
		acc += (rank_1->data[i*rank_1->xdim+k] *
		       rank_1->data[i*rank_1->xdim+k]); 
	}
}
printf("rank_1-%d: %f\n",j+1,acc);
#endif


		FreeArray2D(rank_1);

	}

	Y = DiagonalAverage(sum);
	FreeArray2D(sum);

	FreeArray1D(U_col);
	FreeArray2D(V_row);
	FreeArray2D(U);
	FreeArray2D(V_t);
	FreeArray1D(ev);

	return(Y);
		
}

Array1D **SSARank1(Array1D *series, int L, int start, int end)
{
	Array1D *ev;
	Array2D *ea;
	Array2D *U;
	Array2D *V;
	Array2D *S;
	Array2D *tr;
	Array2D *D;
	Array2D *D2;
	Array2D *Dm2;
	Array2D *X;
	Array2D *X_t;
	Array2D *V_t;
	Array2D *rank_1;
	Array2D *t_u;
	Array2D *t_x;
	Array2D *t_y;
	Array2D *diff;
	Array2D *sum;
	int K;
	int N;
	Array1D *U_col;
	Array1D *V_row;
	int d;
	int i;
	int j;
	int k;
	Array1D **dcomp;
	Array1D *Y;
	int y_i;
	Array2D *I;
	double acc;
	Array1D **rank_arrays;


	N = series->ydim;
	K = N - L + 1;


	X = TrajectoryMatrix(series,0,L,K);
	if(X == NULL) {
		fprintf(stderr, 
			"couldn't create trajctory matrix for %d lags\n",
				L);
		exit(1);
	}
	X_t = TransposeArray2D(X);
	if(X_t == NULL) {
		fprintf(stderr,"couldn't transpose array X\n");
		exit(1);
	}

	S = MultiplyArray2D(X_t,X);
	if(S == NULL) {
		fprintf(stderr,
			"couldn't compute S\n");
		exit(1);
	}
	FreeArray2D(X_t);

	ev = EigenValueArray2D(S);
	if(ev == NULL) {
		fprintf(stderr,"couldn't get eigen values\n");
		exit(1);
	}

	ea = EigenVectorArray2D(S);
	if(ea == NULL) {
		fprintf(stderr,"couldn't get eigen vectors\n");
		exit(1);
	}
	FreeArray2D(S);

	V = UnitizeArray2D(ea);
	if(V == NULL) {
		fprintf(stderr,"couldn't unitize eigen vectors\n");
		exit(1);
	}
	FreeArray2D(ea);
	SortEigenVectors(ev,V);



	t_u = MultiplyArray2D(X,V);
	if(t_u == NULL) {
		exit(1);
	}

	/*
	 * need number of cols from V_t
	 */
	Dm2 = MakeArray2D(t_u->xdim,V->xdim);
	if(Dm2 == NULL) {
		exit(1);
	}

	D2 = MakeArray2D(t_u->xdim,V->xdim);
	if(D2 == NULL) {
		exit(1);
	}

	D = MakeArray2D(t_u->xdim,V->xdim);
	if(D == NULL) {
		exit(1);
	}

	for(i=0; i < Dm2->ydim; i++) {
		for(j=0; j < Dm2->xdim; j++) {
			Dm2->data[i*D->xdim+j] = 0;
			D2->data[i*D->xdim+j] = 0;
			D->data[i*D->xdim+j] = 0;
		}
	}

	/*
	 * diagonal is 1/sqrt(eigenvalue) for non-zero entries
	 *
	 * if L < K, leave others zero
	 */
	for(i=0; i < L; i++) {
		if(i < Dm2->xdim) {
			Dm2->data[i*Dm2->xdim+i] = 1.0/sqrt(ev->data[i]);
			D2->data[i*D2->xdim+i] = sqrt(ev->data[i]);
			D->data[i*D->xdim+i] = ev->data[i];
		}
	}

	U = MultiplyArray2D(t_u,Dm2);
	if(U == NULL) {
		fprintf(stderr,"can't get U\n");
		exit(1);
	}
	FreeArray2D(t_u);

	V_t = TransposeArray2D(V);
	if(V_t == NULL) {
		fprintf(stderr,"couldn't get V_t\n");
		exit(1);
	}
	FreeArray2D(V);
	FreeArray2D(X);


	U_col = MakeArray1D(U->ydim);
	if(U_col == NULL) {
		exit(1);
	}

	V_row = MakeArray2D(1,V_t->xdim);
	if(V_row == NULL) {
		exit(1);
	
	}

	rank_arrays = (Array1D **)malloc((end - start)*sizeof(Array1D *));
	if(rank_arrays == NULL) {
		exit(1);
	}

	for(j = start; j < end; j++) {
		for(i=0; i < U_col->ydim; i++) {
			U_col->data[i] = U->data[i*U->xdim+j] * sqrt(ev->data[j]);
		}
		for(i=0; i < V_row->xdim; i++) {
			V_row->data[0*V_row->xdim + i] =
				V_t->data[j*V_t->xdim+i];
		}

		rank_1 = MultiplyArray2D(U_col,V_row);
		if(rank_1 == NULL) {
			exit(1);
		}

		/*
		Y = DiagonalAverage(rank_1);
		if(Y == NULL) {
			exit(1);
		}
		FreeArray2D(rank_1);
		*/

		rank_arrays[j] = rank_1;

	}


	FreeArray1D(U_col);
	FreeArray2D(V_row);
	FreeArray2D(U);
	FreeArray2D(V_t);
	FreeArray1D(ev);

	return(rank_arrays);
		
}
#ifdef STANDALONE

int SignalRange(char *range, int *start, int *end)
{
	char buf[4096];
	int s;
	int e;
	char *cursor;
	int i;

	
	memset(buf,0,sizeof(buf));
	cursor = range;
	i = 0;
	while(*cursor != 0) {
		if(isspace(*cursor)) {
			cursor++;
			continue;
		}
		if(*cursor == '-') {
			break;
		}
		if(isdigit(*cursor)) {
			buf[i] = *cursor;
			i++;
			if(i == 4096) {
				break;
			}
		}
		cursor++;
	}

	if(i == 4096) {
		return(-1);
	}

	if(i == 0) {
		s = -1;
	} else {
		s = atoi(buf);
	}

	cursor++;
	memset(buf,0,sizeof(buf));
	i = 0;
	while(*cursor != 0) {
		if(isspace(*cursor)) {
			cursor++;
			continue;
		}
		if(isdigit(*cursor)) {
			buf[i] = *cursor;
			i++;
			if(i == 4096) {
				break;
			}
		}
		cursor++;
	}

	if(i == 4096) {
		return(-1);
	}

	/*
	 * only one value signifies end of range -- not start
	 */
	if(i == 0) {
		e = s;
		s = -1;
	} else {
		e = atoi(buf);
	}

	*start = s;
	*end = e;

	return(1);
}
		

Object *InitObject(int id, Array1D *series, int L, int K)
{
	Object *lob;

	lob = (Object *)Malloc(sizeof(Object));
	if(lob == NULL) {
		exit(1);
	}

	lob->series = DiagonalAverage(series);
	lob->id = id;
	lob->L = L;
	lob->K = K;

	return(lob);
}

void FreeObject(Object *ob)
{
	FreeArray1D(ob->series);
	Free(ob);
	return;
}

Bin *InitBin(int len)
{
	Bin *b;
	int i;

	b = (Bin *)malloc(sizeof(Bin));
	if(b == NULL)
	{
		exit(1);
	}
	memset(b,0,sizeof(Bin));

	b->count = 0;
	b->list = DlistInit();
	if(b->list == NULL)
	{
		exit(1);
	}

	b->centroid = NULL;

	return(b);
}

void FreeBin(Bin *b)
{
	Object *ob;

	DlistRemove(b->list);

	FreeArray1D(b->centroid);
	Free(b);

	return;
}

void PrintBin(FILE *fd, Bin *b)
{
        DlistNode *d;
        Object *ob;
        double min;
        double max;

	DLIST_FORWARD(b->list,d)
	{
		ob = (Object *)d->value.v;
		fprintf(fd,"\t%d\n",
				ob->id);
	}
	fprintf(fd,"\n");

        return;
}


void ComputeCentroid(Bin *b, int L, int K)
{
	double dist;
	Array2D *tr_cent;
	Array2D *tr_dst;
	Array2D *r_cent;
	int i;
	int j;
	DlistNode *d;
	Object *ob;

	tr_cent = NULL;
	DLIST_FORWARD(b->list,d)
        {
                ob = (Object *)d->value.v;
		tr_dst = TrajectoryMatrix(ob->series,0,L,K);
		if(tr_dst == NULL) {
			exit(1);
		}
		if(tr_cent == NULL) {
			tr_cent = MakeArray2D(tr_dst->ydim,tr_dst->xdim);
			if(tr_cent == NULL) {
				exit(1);
			}
			for(i=0; i < tr_cent->ydim; i++) {
				for(j=0; j < tr_cent->xdim; j++) {
					tr_cent->data[i*tr_cent->xdim+j] = 0; 
				}
			}
		}
		for(i=0; i < tr_dst->ydim; i++) {
			for(j=0; j < tr_dst->xdim; j++) {
				tr_cent->data[i*tr_cent->xdim+j] +=
					tr_dst->data[i*tr_dst->xdim+j];
			}
		}
		FreeArray2D(tr_dst);
	}


	r_cent = DiagonalAverage(tr_cent);
	if(r_cent == NULL) {
		exit(1);
	}
	FreeArray2D(tr_cent);
	for(i=0; i < r_cent->ydim; i++) {
		for(j=0; j < r_cent->xdim; j++) {
			r_cent->data[i*r_cent->xdim+j] /=
				(double)(b->count);
		}
	}


	if(b->centroid != NULL) {
		FreeArray1D(b->centroid);
	}
	b->centroid = r_cent;

	return;
}

void AddObjectToBin(Bin *b, Object *ob)
{
	int count;
	DlistNode *d;
	int i;
	double old_val;
	double new_val;
	Array2D *r_series;

	count = b->count;
	b->count++;
	DlistAppend(b->list,(Hval)(void *)ob);

	return;
}

double Omega(int i, int L, int N)
{
        int K;
        int L_star;
        int K_star;
        int omega;

        K = N - L + 1;

        if(L < K) {
                L_star = L;
                K_star = K;
        } else {
                L_star = K;
                K_star = L;
        }

        if((i >=0) && (i <= (L_star - 1))) {
                omega = i+1;
        } else if((i <= L_star) && (i <= K_star)) {
                i = L_star;
        } else {
                i = N - 1;
        }

        return(omega);
}

double InnerProduct(Array1D *f_1, Array1D *f_2, int L, int K)
{
        double sum;
        int i;
        int N;

        N = L + K - 1;

        sum = 0;
        for(i=0; i < N; i++) {
                sum += (Omega(i,L,K) * f_1->data[i] * f_2->data[i]);
        }

        return(sum);

}

double Corr(Array1D *x, Array1D *y)
{
	int i;
	double x_bar;
	double y_bar;
	double x2;
	double y2;
	double count;
	double acc;
	double num;
	double r;

	acc = 0;
	count = 0;

	for(i=0; i < x->ydim; i++) {
		acc += x->data[i];
		count++;
	}
	x_bar = acc / count;

	acc = 0;
	count = 0;
	for(i=0; i < y->ydim; i++) {
		acc += y->data[i];
		count++;
	}
	y_bar = acc / count;

	acc = 0;
	count = 0;
	for(i=0; i < x->ydim; i++) {
		acc += ((x->data[i] - x_bar) * (y->data[i] - y_bar));
		count++;
	}
	num = acc / count;

	acc = 0;
	count = 0;
	for(i=0; i < x->ydim; i++) {
		acc += ((x->data[i] - x_bar) * (x->data[i] - x_bar));
		count++;
	}
	x2 = acc / count;

	acc = 0;
	count = 0;
	for(i=0; i < y->ydim; i++) {
		acc += ((y->data[i] - y_bar) * (y->data[i] - y_bar));
		count++;
	}
	y2 = acc / count;

	r = num / sqrt(x2 * y2);

	return(fabs(r));
} 
		


double WCorr(Array2D *tr_x, Array2D *tr_y, int L, int K) {
        double rho;
        double ip1;
        double ip2;
        double ip3;

        ip1 = InnerProduct(tr_x, tr_y, L, K);
        ip2 = InnerProduct(tr_x, tr_x, L, K);
        ip3 = InnerProduct(tr_y, tr_y, L, K);

        rho = ip1 / ((sqrt(ip2) * sqrt(ip3)));

        return(rho);
}

/*
 * src is centrod, dst is tr matrix
 */
double DistanceW(Bin *b, Array1D *dst, int L, int K)
{
	double dist;
	Array2D *tr_cent;
	Array2D *tr_dst;
	int i;
	int j;
	DlistNode *d;
	Object *ob;
	int adj = 0;

	tr_cent = NULL;
	DLIST_FORWARD(b->list,d)
        {
                ob = (Object *)d->value.v;
		if(ob->series == dst) {
			adj = 1;
			continue;
		}
		/*
		tr_dst = TrajectoryMatrix(ob->series,0,L,K);
		if(tr_dst == NULL) {
			exit(1);
		}
		*/
		tr_dst = ob->series;
		if(tr_cent == NULL) {
			tr_cent = MakeArray2D(tr_dst->ydim,tr_dst->xdim);
			if(tr_cent == NULL) {
				exit(1);
			}
			for(i=0; i < tr_cent->ydim;i++) {
				for(j=0; j < tr_cent->xdim; j++) {
					tr_cent->data[i*tr_cent->xdim+j] = 
						tr_dst->data[i*tr_dst->xdim+j];
				}
			}
			
		} else {
			for(i=0; i < tr_dst->ydim; i++) {
				for(j=0; j < tr_dst->xdim; j++) {
					tr_cent->data[i*tr_cent->xdim+j] +=
						tr_dst->data[i*tr_dst->xdim+j];
				}
			}
			/*
			for(i=0; i < tr_dst->ydim; i++) {
				for(j=0; j < tr_dst->xdim; j++) {
					tr_cent->data[i*tr_cent->xdim+j] /=
						(double)(b->count-adj);
				}
			}
			*/
		}
	}

	if(tr_cent == NULL) {
		return(1.0);
	}


//	dist = 1.0 - WCorr(tr_cent,dst,L,K);
	dist = WCorr(tr_cent,dst,L,K);
	FreeArray2D(tr_cent);
	return(fabs(dist));
}

double DistanceC(Bin *b, Array1D *dst, int L, int K)
{
	double dist;
	Array2D *tr_cent;
	Array2D *tr_dst;
	Array2D *r_series;
	Array2D *r_cent;
	int i;
	int j;
	DlistNode *d;
	Object *ob;


	dist = 1.0 - Corr(b->centroid,dst);
//	dist = Corr(r_cent,r_series);
	return(fabs(dist));
}


int IsEqCentroid(Bin *b1, Bin *b2)
{
	int i;

	for(i=0; i < b1->centroid->ydim; i++) {
		if(b1->centroid->data[i] != b2->centroid->data[i]) {
printf("%d %e %e\n",i,b1->centroid->data[i],b2->centroid->data[i]);
			return(0);
		}
	}

	return(1);

}

int IsDone(Bin **bins, Bin **new_bins, int means)
{
	int i;

	for(i=0; i < means; i++)
	{
		if((bins[i]->count == 0) && (new_bins[i]->count == 0))
		{
			continue;
		}
		if(!IsEqCentroid(bins[i],new_bins[i])) {
			return(0);
		}
	}

	return(1);
}

		

Bin **KMeans(Array1D **decomp, int L, int K, int means)
{
	int i;
	int j;
	Bin **bins;
	Bin **new_bins;
	Object *ob;
	DlistNode *d;
	int done;
	double min_dist;
	int min_j;
	double dist;
	struct timeval tm;


	/*
	 * make a set of k bins
	 */
	bins = (Bin **)Malloc(means*sizeof(Bin *));
	if(bins == NULL)
	{
		exit(1);
	}

	for(i=0; i < means; i++)
	{
		bins[i] = InitBin(decomp[0]->ydim);
		if(bins[i] == NULL)
		{
			exit(1);
		}
	}

	/*
	 * initially, scatter randomly
	 */
	gettimeofday(&tm,NULL);
	srand48(tm.tv_sec+tm.tv_usec);
	for(i=0; i < L; i++) {
		ob = InitObject(i,decomp[i],L,K);
		if(ob == NULL) {
			exit(1);
		}
		j = RAND() * means;
		AddObjectToBin(bins[j],ob);
		ComputeCentroid(bins[j],L,K);
	}

	done = 0;
	while(!done)
	{
		new_bins = (Bin **)Malloc(means*sizeof(Bin *));
		if(new_bins == NULL)
		{
			exit(1);
		}

		for(i=0; i < means; i++)
		{
			new_bins[i] = InitBin(decomp[0]->ydim);
			if(new_bins[i] == NULL)
			{
				exit(1);
			}
		}

		/*
		 * traverse the existing data
		 */
		for(i=0; i < means; i++)
		{
			DLIST_FORWARD(bins[i]->list,d)
			{
				ob = (Object *)d->value.v;
				/*
				 * find the destination bin that is
				 * closest
				 */
				min_dist = 999999999999999999999.99;
				min_j = 0;
				for(j=0; j < means; j++)
				{
					dist = DistanceC(bins[j],
						        ob->series,L,K);
					if(dist < min_dist)
					{
						min_dist = dist;
						min_j = j;
					}
				}
printf("moving %d from %d to %d, dist: %f\n",ob->id,i,min_j,dist);
				AddObjectToBin(new_bins[min_j],ob);
				ComputeCentroid(new_bins[min_j],L,K);
			}
		}

		/*
		 * if we have converged, bail out
		 */
		if(IsDone(bins,new_bins,means))
		{
			for(i=0; i < means; i++)
			{
				FreeBin(bins[i]);
			}
			free(bins);
			return(new_bins);
		}
		/*
		 * otherwise, make new_bin into bins and repeat
		 */
		for(i=0; i < means; i++)
		{
			FreeBin(bins[i]);
		}
		Free(bins);
		bins = new_bins;
				
	}

	return(new_bins);
}

	

#define ARGS "x:l:N:e:p:"
char *Usage = "usage: ssa-decomp -x xfile\n\
\t-e range of signal series (start - end || end)\n\
\t-l lags\n\
\t-N number of time series elements in base matrix\n\
\t-p starting index\n";

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
	int i;
	int j;
	double rho;
	int K;
	int N;
	int signal;
	Array1D *series_signal;
	double x_n;
	int err;
	int start;
	int end;
	Bin **bins;
	Array1D **rank_arrays;
	

	N = 0;
	p = 0;
	start = -1;
	end = -1;
	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'e':
				err = SignalRange(optarg,&start,&end);
				break;
			case 'l':
				lags = atoi(optarg);
				break;
			case 'p':
				p = atoi(optarg);
				break;
			case 'N':
				N = atoi(optarg);
				break;
			default:
				fprintf(stderr,
			"unrecognized command: %c\n",(char)c);
				fprintf(stderr,"%s",Usage);
				exit(1);
		}
	}

	if(err < 0) {
		fprintf(stderr,"error in -e argument\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(Xfile[0] == 0) {
		fprintf(stderr,"must specify xfile\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(N == 0) {
		fprintf(stderr,"must enter window size\n");
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

	if(start == -1) {
		start = 0;
		rank_arrays = SSARank1(x,lags,0,lags);
		if(rank_arrays == NULL) {
			exit(1);
		}
		K = x->ydim - lags + 1;
		bins = KMeans(rank_arrays,lags,K,3);
		for(i=0; i < 3; i++) {
			printf("Bin: %d\n",i);
			PrintBin(stdout,bins[i]);
		}
		free(rank_arrays);
		free(bins);
	} else {
		if(start > 0) {
			start = start - 1;
		}
	}

	if(end == -1) {
		end = lags - 1;
	}


	if(p == 0) { // if p == 0, start from beginning of x
		series_signal = SSADecomposition(x,lags,start,end);
	} else {
		x_sub = MakeArray1D(x->ydim - p);
		if(x_sub == NULL) {
			exit(1);
		}
		for(i=p; i < x->ydim; i++) {
			x_sub->data[i-p] = x->data[i];
		}
		series_signal = SSADecomposition(x_sub,lags,start,end);
		FreeArray1D(x_sub);
	}

	for(i=0; i < series_signal->ydim; i++) {
		printf("%f\n",series_signal->data[i]);
	}


	return(0);
}

#endif
