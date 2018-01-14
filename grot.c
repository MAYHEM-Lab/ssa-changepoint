
int Eigenize(Array1D *series, int L, int K, Array1D **o_ev, Array2D *o_V)
{
	Array1D *ev;
	Array2D *ea;
	Array2D *V;
	Array2D *S;
	Array2D *X;
	Array2D *X_t;
	int K;
	int N;


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
	FreeArray2D(X);

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

	*o_ev = ev;
	*o_V = V;

	return(1);
}

