all: psrs.out
CC=mpicc
# CCFLAG="-std=C99"
RUN=mpirun
psrs.out: psrs.o
	${CC} -o $@ $^
psrs.o: psrs.c
	${CC} ${CCFLAG} -c -o $@ $< 

run: psrs.out
	${RUN} -n 4 ./$< "../data/input" 26

clean:
	find . -type f -name "*.o" -delete
	find . -type f -name "*.out" -delete
