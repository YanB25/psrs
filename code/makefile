all: psrs.out
CC=mpicc
# CCFLAG="-std=C99"
RUN=mpirun
psrs.out: psrs.o
	${CC} -o $@ $^
psrs.o: psrs.c
	${CC} ${CCFLAG} -c -o $@ $< 

run: psrs.out
	${RUN} -n 3 ./$< "../data/input" 26

