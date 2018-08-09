CC=gcc
MPATH=../mio
APATH=../matrix
DPATH=../distributions
EPATH=../euca-cutils
LPATH=${APATH}/lapack-3.8.0/LAPACKE/include/
CFLAGS=-g -I${MPATH} -I${APATH} -I${LPATH} -I${EPATH} -I${DPATH} -I/usr/local/include

#LLIB=-L${APATH}/lapack-3.5.0 -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5/ -llapacke -llapack -lrefblas -lblas -ltmglib -lgfortran
#LLIB=-L./lapack-3.5.0 -L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/ -llapacke -llapack -lrefblas -lblas -ltmglib -lgfortran
LLIB=-L${APATH}/lapack-3.8.0 -L/usr/local/Cellar/gcc/8.1.0/lib/gcc/8/ -llapacke -llapack -lblas -lgfortran -lpthread
ALIB=${APATH}/mioarray.o
LIBS=${MPATH}/mymalloc.o ${MPATH}/mio.o ${EPATH}/libutils.a -lm ${DPATH}/normal.o ${ALIB} ${LLIB}

#for OSX and homebrew
#LLIB=-L/usr/local/opt/lapack/lib -L/usr/local/opt/openblas/lib -llapacke -llapack -lopenblas
#LPATH=/usr/local/opt/lapack/include
#BPATH=/usr/local/opt/openblas/include
#CFLAGS=-g -I${MPATH} -I${APATH} -I${BPATH} -I${LPATH} -I${EPATH} -I${DPATH} -I/usr/local/include


all: ssa-cp ssa-decomp ssa-basis

ssa-cp: ssa-cp.c
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o ssa-cp ssa-cp.c ${LIBS}

ssa-decomp: ssa-decomp.c
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o ssa-decomp ssa-decomp.c ${LIBS}

ssa-basis: ssa-basis.c
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o ssa-basis ssa-basis.c ${LIBS}

clean:
	rm *.o ssa-cp ssa-decomp ssa-basis
