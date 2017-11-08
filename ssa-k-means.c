#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "dlist.h"
#include "mio.h"
#include "mioarray.h"

#define RAND() (drand48())

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

char Line_buf[1024*1024];
int Object_id;

Object *InitObject(int id, Array1D *series, int L, int K)
{
	Object *lob;

	lob = (Object *)Malloc(sizeof(Object));
	if(lob == NULL) {
		exit(1);
	}

	lob->series = series;
	lob->id = id;
	lob->L = L;
	lob->K = K;

	return(lob);
}

void FreeObject(Object *ob)
{
	Free(ob);
}


void FreeObject(Object *ob)
{
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

	b->centroid = MakeArray1D(len);
	if(b->centroid == NULL) {
		exit(1);
	}

	for(i-0; i < b->centroid->ydim; i++) {
		b->centroid->data[i] = 0;
	}

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

        if(Verbose == 1)
        {
                DLIST_FORWARD(b->list,d)
                {
                        ob = (Object *)d->value.v;
                        fprintf(fd,"\t%d\n",
                                        ob->id);
                }
                fprintf(fd,"\n");
        }

        return;
}


AddObjectToBin(Bin *b, Object *ob)
{
	int count;
	Object *ob;
	DlistNode *d;
	int i;

	count = b->count;
	b->count++;
	DlistAppend(b->list,(Hval)(void *)ob);

	for(i=0; i < b->centroid->ydim; i++) {
		old_val = b->centroid->data[i] * (double)count;
		new_val = old_val + ob->series->data[i];
		b->centroid->data[i] = new_val / (double)(b->count);
	}

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

double WCorr(Array1D x, Array1D *y, int L, int K) {
        double rho;
        double ip1;
        double ip2;
        double ip3;
	Array2D *tr_x;
	Array2D *tr_y;

	tr_x = TrajectoryMatrix(x,0,L,K);
	if(tr_x == NULL) {
		exit(1);
	}

	tr_y = TrajectoryMatrix(y,0,L,K);
	if(tr_y == NULL) {
		exit(1);
	}

        ip1 = InnerProduct(tr_x, tr_y, L, K);
        ip2 = InnerProduct(tr_x, tr_x, L, K);
        ip3 = InnerProduct(tr_y, tr_y, L, K);

        rho = ip1 / ((sqrt(ip2) * sqrt(ip3)));

	FreeArray2D(tr_x);
	FreeArray2D(tr_y);

        return(rho);
}

double Distance(Array1D *src, Array1D *dst, int L, int K)
{
	double dist;

	dist = 1.0 - WCorr(src,dst,L,K);
	return(dist);
}

int IsEqCentroid(Bin *b1, Bin *b2)
{
	int i;

	for(i=0; i < b1->centroid->ydim; i++) {
		if(b1->centroiddata[i] != b2->centroid->data[i]) {
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
	for(i=0; i < decomp->ydim; i++) {
	{
		ob = InitObject(i,decomp[i],L,K);
		if(ob == NULL) {
			exit(1);
		}
		j = RAND() * means;
		AddObjectToBin(bins[i],ob);
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
					dist = Distance(bins[j]->centroid,
						        ob->series,L,K);
					if(dist < min_dist)
					{
						min_dist = dist;
						min_j = j;
					}
				}
				AddObjectToBin(new_bins[min_j],ob);
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
			done = 1;
			break;
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

