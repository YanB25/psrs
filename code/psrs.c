#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
int* scatterv_size(int size, int world_size) {
    /* 
    * return the scatter size used in mpi_scatterv
    * @param size, length of the array
    * @param world_size, mpi communicator size
    * @return ret, int* of length world_size. remember to free it.
    */
    int* ret = (int*) malloc(sizeof(int) * world_size);
    int base_size = size / world_size;
    int more_size = size % world_size;
    for (int i = 0; i < world_size - more_size; ++i) {
        ret[i] = base_size;
    }
    for (int i = world_size - more_size; i < world_size; ++i) {
        ret[i] = base_size + 1;
    }
    return ret;
}

int* scatterv_dipl(int* scatterv_array, int world_size) {
    int* ret = (int*) malloc(sizeof(int) * world_size);
    ret[0] = 0;
    for (int i = 1; i < world_size; ++i) {
        ret[i] = ret[i-1] + scatterv_array[i-1];
    }
    // for (int i = 0; i < world_size; ++i) {
    //     printf("%d ", ret[i]);
    // }
    // printf("\n");
    return ret;
}

int cmp(const void* lhs, const void* rhs) {
    return *((int*)lhs) - *((int*) rhs);
}

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int* array = NULL;
    int* my_array = NULL;
    int SIZE;
    sscanf(argv[2], "%d", &SIZE);

    int* scatterv_size_array = scatterv_size(SIZE, world_size);
    int* scatterv_dipl_array = NULL;

    if (my_rank == 0) {
        scatterv_dipl_array = scatterv_dipl(scatterv_size_array, world_size);

        FILE* inputfile;
        inputfile = fopen(argv[1], "r");
        if (!inputfile) {
            printf("file not found %s, abort.\n", argv[1]);
            return 1;
        }

        array = (int*) malloc(sizeof(int) * SIZE);

        for (int i = 0; i < SIZE; ++i) {
            fscanf(inputfile, "%d", &array[i]);
        }
    }
    my_array = (int*) malloc(sizeof(int) * scatterv_size_array[my_rank]);

    // step1: divide averagely
    MPI_Scatterv(array, scatterv_size_array,
        scatterv_dipl_array, MPI_INT, my_array, 
        scatterv_size_array[my_rank], MPI_INT, 0, MPI_COMM_WORLD);

    // step2: local sort
    qsort((void*)my_array, scatterv_size_array[my_rank], 
        sizeof(int), cmp);

    for (int i = 0; i < scatterv_size_array[my_rank]; ++i) {
        printf("[%d]: %d\n", my_rank, my_array[i]);
    }

    // step3: choose pivot
    int* pivot_array = (int*) malloc(sizeof(int) * world_size);
    int step_size = scatterv_size_array[my_rank] / world_size;
    for (int i = 0; i < world_size; ++i) {
        pivot_array[i] = my_array[i * step_size];
        printf("[%d] pivot %d\n", my_rank, pivot_array[i]);
    }

    // step4: gather all and sort all the samples from p
    int* all_pivot = (int*) malloc(sizeof(int) * world_size * world_size);
    MPI_Allgather(pivot_array, world_size, MPI_INT, all_pivot, world_size, MPI_INT, MPI_COMM_WORLD);
    qsort((void*)all_pivot, world_size*world_size, sizeof(int), cmp);
    if (my_rank == 0) {
        for (int i = 0; i < world_size * world_size; ++i) {
            printf("all pivot %d\n", all_pivot[i]);
        }
    }
    free(pivot_array);
    free(scatterv_dipl_array);
    free(scatterv_size_array);
    free(array);
    free(my_array);
    MPI_Finalize();
}