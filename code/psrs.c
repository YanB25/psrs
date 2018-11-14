#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define ul_t unsigned long
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
    assert(world_size > more_size);
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
    //     printf("%ud ", ret[i]);
    // }
    // printf("\n");
    return ret;
}

int cmp(const void* lhs, const void* rhs) {
    return *((ul_t*)lhs) - *((ul_t*) rhs);
}

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    assert(world_size > 0);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    ul_t* array = NULL;
    ul_t* my_array = NULL;
    ul_t SIZE;
    sscanf(argv[2], "%lu", &SIZE);

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

        array = (ul_t*) malloc(sizeof(ul_t) * SIZE);

        for (ul_t i = 0; i < SIZE; ++i) {
            fscanf(inputfile, "%lu", &array[i]);
        }
    }
    my_array = (ul_t*) malloc(sizeof(ul_t) * scatterv_size_array[my_rank]);

    // step1: divide averagely
    MPI_Scatterv(array, scatterv_size_array,
        scatterv_dipl_array, MPI_UNSIGNED_LONG, my_array, 
        scatterv_size_array[my_rank], MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    free(array);
    free(scatterv_dipl_array);

    // step2: local sort
    qsort((void*)my_array, scatterv_size_array[my_rank], 
        sizeof(ul_t), cmp);

    for (int i = 0; i < scatterv_size_array[my_rank]; ++i) {
        printf("[%d]: %lu\n", my_rank, my_array[i]);
    }

    // step3: choose pivot
    ul_t* pivot_array = (ul_t*) malloc(sizeof(ul_t) * world_size);
    ul_t step_size = scatterv_size_array[my_rank] / world_size;
    for (ul_t i = 0; i < world_size; ++i) {
        pivot_array[i] = my_array[i * step_size];
        printf("[%d] pivot %lu\n", my_rank, pivot_array[i]);
    }

    // step4: gather all and sort all the samples from p
    ul_t* all_pivot = (ul_t*) malloc(sizeof(ul_t) * world_size * world_size);
        // main privot is for (step 5).
    ul_t* main_pivot = (ul_t*) malloc(sizeof(ul_t) * (world_size-1));
    MPI_Gather(pivot_array, world_size, MPI_UNSIGNED_LONG, all_pivot, world_size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
        qsort((void*)all_pivot, world_size*world_size, sizeof(ul_t), cmp);
        for (ul_t i = 0; i < world_size * world_size; ++i) {
            printf("all pivot %lu\n", all_pivot[i]);
        }

        //step5: send main pivot to all process
        for (int i = 0; i < world_size - 1; ++i) {
            main_pivot[i] = all_pivot[(i + 1) * world_size];
        }
    }
    MPI_Bcast(main_pivot, world_size - 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    for (int i = 0; i < world_size-1; ++i) {
        printf("I am [%d] and get main pivot %lu\n", my_rank, main_pivot[i]);
    }

    // step6: devide by main pivot.
    int* send_count_array = (int*) malloc(sizeof(int) * world_size);
    memset(send_count_array, 0, sizeof(int) * world_size);
    int* sdispls_array = (int*) malloc(sizeof(int) * world_size);
    memset(sdispls_array, 0, sizeof(int) * world_size);
    int mySize = scatterv_size_array[my_rank];
    printf("[%d] my array size is %d\n", my_rank, mySize);
    int myidx = 0;
    sdispls_array[0] = 0;
    int before_count_sum = 0;
    for (int i = 0; i < mySize; ++i) {
        if (myidx < world_size - 1 && my_array[i] >= main_pivot[myidx]) {
            // printf("debug [%d], myarray[%d] = %lu, main_pivot[%d] = %lu\n", my_rank, i, my_array[i], myidx, main_pivot[myidx]);
            send_count_array[myidx] = i - before_count_sum;
            before_count_sum += send_count_array[myidx];

            // sdispls_array[myidx + 1] = i; //not used anymore
            while (my_array[i] >= main_pivot[myidx] && myidx < world_size - 1) {
                myidx++;
            }
        }
    }
    int before_disp_sum = 0;
    for (int i = 0; i < world_size; ++i) {
        sdispls_array[i] = before_disp_sum;
        before_disp_sum += send_count_array[i];
    }
    send_count_array[world_size - 1] = mySize - before_count_sum;
    for (int i = 0; i < world_size; ++i) {
        // printf("[%d] send_count_array[%d] %d\n", my_rank, i, send_count_array[i]);
        // printf("[%d] sdispls array[%d] %d\n", my_rank, i, sdispls_array[i]);
    }

    // step7: global swap.
    // firstly all to all the counts in all processes
    int* all_count_array = (int*) malloc(sizeof(int) * world_size);
    MPI_Alltoall(send_count_array, 1, MPI_INT, all_count_array,
        1, MPI_INT, MPI_COMM_WORLD);
    // get how many numbers self process should get.
    int my_recv_count = 0;
    for (int i = 0; i < world_size; ++i) {
        my_recv_count += all_count_array[i];
    }
    // used in alltoallv function's rdispls param. so need to accumulate.
    int* rdispls = (int*) malloc(sizeof(int) * world_size);
    rdispls[0] = 0;
    for (int i = 1; i < world_size; ++i) {
        rdispls[i] = rdispls[i - 1] + all_count_array[i - 1];
    }
    for (int i = 0; i < world_size; ++i) {
        printf("from [%d]: all_count_array[%d] = %d\n", my_rank, i, all_count_array[i]);
        printf("from [%d]: rdispls[%d] = %d\n", my_rank, i, rdispls[i]);
    }
    printf("[%d] I should receiv size %d\n", my_rank, my_recv_count);
    ul_t* my_swap_array = (ul_t*) malloc(sizeof(ul_t) * my_recv_count);

    MPI_Alltoallv(my_array, send_count_array, sdispls_array, MPI_UNSIGNED_LONG,
        my_swap_array, all_count_array, rdispls, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    free(all_count_array);

    for (int i = 0; i < my_recv_count; ++i) {
        printf("[%d] recv %lu\n", my_rank, my_swap_array[i]);
    }


    free(rdispls);

    free(sdispls_array);
    free(send_count_array);
    
    free(main_pivot);
    free(all_pivot);
    free(pivot_array);
    free(scatterv_size_array);
    free(my_array);
    MPI_Finalize();
}