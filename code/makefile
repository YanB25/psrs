all: psrs.out
CC=mpic++
CCFLAG=-Ofast -std=c++11
RUN=mpirun
psrs.out: psrs.o 
	${CC} -o $@ $^
psrs.o: psrs.cpp
	${CC} ${CCFLAG} -c -o $@ $< 

run: psrs.out
	${RUN} -n 5 ./$< "binary_input"

gendata.out: gendata.c
	gcc -o $@ $<
gen:
	./gendata.out

testread.out: testread.c
	gcc -o $@ $<
test:
	./testread.out

clean:
	find . -type f -name "*.o" -delete
	find . -type f -name "*.out" -delete
