#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <queue>
using std::priority_queue;
#define ul_t unsigned long
struct State {
    ul_t data;
    int queue_id;
    State (ul_t d, int id): data(d), queue_id(id) {}    
    bool operator<(const State& rhs) const {
        return data > rhs.data;
    }
};

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

    double begin_time = MPI_Wtime();
    if (my_rank == 0) {
        printf("begin.\n");
    }
    printf("[%d] is working\n", my_rank);

    ul_t* array = NULL;
    ul_t* my_array = NULL;
    ul_t SIZE;

    FILE* inputfile;
    inputfile = fopen(argv[1], "rb");
    if (!inputfile) {
        printf("ERROR ! no such file %s\n", argv[1]);
    }
    fread(&SIZE, sizeof(ul_t), 1, inputfile);

    int* scatterv_size_array = scatterv_size(SIZE, world_size);
    int* scatterv_dipl_array = NULL;

    if (my_rank == 0) {
        scatterv_dipl_array = scatterv_dipl(scatterv_size_array, world_size);

        fseek(inputfile, sizeof(ul_t), SEEK_SET);

        array = (ul_t*) malloc(sizeof(ul_t) * SIZE);

        fread(array, sizeof(ul_t), SIZE, inputfile);
    }
    my_array = (ul_t*) malloc(sizeof(ul_t) * scatterv_size_array[my_rank]);

    double step0 = MPI_Wtime();
    if (my_rank == 0) {
        printf("finished step 0 data prepared. elapse %lf\n", step0-begin_time);
    }

    // step1: divide averagely
    MPI_Scatterv(array, scatterv_size_array,
        scatterv_dipl_array, MPI_UNSIGNED_LONG, my_array, 
        scatterv_size_array[my_rank], MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    free(array);
    free(scatterv_dipl_array);

    double step1 = MPI_Wtime();
    if (my_rank == 0) {
        printf("fnished step 1 divide averagely. elapse %lf\n", step1-begin_time);
    }

    // step2: local sort
    qsort((void*)my_array, scatterv_size_array[my_rank], 
        sizeof(ul_t), cmp);

    // for (int i = 0; i < scatterv_size_array[my_rank]; ++i) {
        // printf("[%d]: %lu\n", my_rank, my_array[i]);
    // }

    double step2 = MPI_Wtime();
    if (my_rank == 0) {
        printf("finished step 2 local sort. elapse %lf\n", step2-begin_time);
    }

    // step3: choose pivot
    ul_t* pivot_array = (ul_t*) malloc(sizeof(ul_t) * world_size);
    ul_t step_size = scatterv_size_array[my_rank] / world_size;
    for (ul_t i = 0; i < world_size; ++i) {
        pivot_array[i] = my_array[i * step_size];
        // printf("[%d] pivot %lu\n", my_rank, pivot_array[i]);
    }

    double step3 = MPI_Wtime();
    if (my_rank == 0) {
        printf("finished step 3 choose pivot. elapse %lf\n", step3-begin_time);
    }

    // step4: gather all and sort all the samples from p
    ul_t* all_pivot = (ul_t*) malloc(sizeof(ul_t) * world_size * world_size);
        // main privot is for (step 5).
    ul_t* main_pivot = (ul_t*) malloc(sizeof(ul_t) * (world_size-1));
    MPI_Gather(pivot_array, world_size, MPI_UNSIGNED_LONG, all_pivot, world_size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
        qsort((void*)all_pivot, world_size*world_size, sizeof(ul_t), cmp);
        for (ul_t i = 0; i < world_size * world_size; ++i) {
            // printf("all pivot %lu\n", all_pivot[i]);
        }

        double step4 = MPI_Wtime();
        printf("finished step 4 gather all and sort. elapse %lf\n", step4-begin_time);

        //step5: send main pivot to all process
        for (int i = 0; i < world_size - 1; ++i) {
            main_pivot[i] = all_pivot[(i + 1) * world_size];
        }
    }
    MPI_Bcast(main_pivot, world_size - 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    for (int i = 0; i < world_size-1; ++i) {
        // printf("I am [%d] and get main pivot %lu\n", my_rank, main_pivot[i]);
    }

    double step5 = MPI_Wtime();
    if (my_rank == 0) {
        printf("finished step 5 send main privot. elapse %lf\n", step5 - begin_time);
    }

    // step6: devide by main pivot.
    int* send_count_array = (int*) malloc(sizeof(int) * world_size);
    memset(send_count_array, 0, sizeof(int) * world_size);
    int* sdispls_array = (int*) malloc(sizeof(int) * world_size);
    memset(sdispls_array, 0, sizeof(int) * world_size);
    int mySize = scatterv_size_array[my_rank];
    // printf("[%d] my array size is %d\n", my_rank, mySize);
    int myidx = 0;
    sdispls_array[0] = 0;
    // int before_count_sum = 0;
    for (int i = 0; i < mySize; ++i) {
        if (my_array[i] < main_pivot[myidx]) {
            send_count_array[myidx]++;
            // before_count_sum++;
        } else {
            // array >= index
            while (myidx < world_size && my_array[i] >= main_pivot[myidx]) {
                myidx++;
            }
            send_count_array[myidx]++;
            // before_count_sum++;
        }
    }
    int before_disp_sum = 0;
    for (int i = 0; i < world_size; ++i) {
        sdispls_array[i] = before_disp_sum;
        before_disp_sum += send_count_array[i];
    }
    // send_count_array[world_size - 1] = mySize - before_count_sum;
    // for (int i = 0; i < world_size; ++i) {
    //     // printf("[%d] send_count_array[%d] %d\n", my_rank, i, send_count_array[i]);
    //     // printf("[%d] sdispls array[%d] %d\n", my_rank, i, sdispls_array[i]);
    // }

    double step6 = MPI_Wtime();
    if (my_rank == 0) {
        printf("finihsed step 6 devide by main privot. elapse %lf\n", step6-begin_time);
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
    // for (int i = 0; i < world_size; ++i) {
    //     // printf("from [%d]: all_count_array[%d] = %d\n", my_rank, i, all_count_array[i]);
    //     // printf("from [%d]: rdispls[%d] = %d\n", my_rank, i, rdispls[i]);
    // }
    // printf("[%d] I should receiv size %d\n", my_rank, my_recv_count);
    ul_t* my_swap_array = (ul_t*) malloc(sizeof(ul_t) * my_recv_count);

    MPI_Alltoallv(my_array, send_count_array, sdispls_array, MPI_UNSIGNED_LONG,
        my_swap_array, all_count_array, rdispls, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    printf("[%d] recv count %d\n", my_rank, my_recv_count);
    // for (int i = 0; i < my_recv_count; ++i) {
    //     printf("[%d] recv %lu\n", my_rank, my_swap_array[i]);
    // }

    double step7 = MPI_Wtime();
    if (my_rank == 0) {
        printf("finished step 7 global swap. elapse %lf\n", step7-begin_time);
    }

    // step8: n-way merge sort.
    int* queue_beg = (int*) malloc(sizeof(int) * world_size);
    memcpy(queue_beg, rdispls, sizeof(int) * world_size);

    priority_queue<State> pq;
    int idx = 0;
    for (int i = 0; i < world_size; ++i) {
        if (all_count_array[i] != 0) {
            pq.push(State(my_swap_array[rdispls[i]], i));
        }
    }


    ul_t* result = (ul_t*) malloc(sizeof(ul_t) * my_recv_count);
    idx = 0;
    while (!pq.empty()) {
        State top = pq.top(); pq.pop();
        // printf("-- [%d] %lu : %d\n", my_rank, top.data, top.queue_id);
        ul_t data = top.data;
        int queue_id = top.queue_id;
        result[idx++] = data;
        queue_beg[queue_id]++;
        bool notFull = (queue_id != world_size - 1 && queue_beg[queue_id] < rdispls[queue_id + 1])
            || (queue_id == world_size - 1 && queue_beg[queue_id] < my_recv_count);
        if (notFull) {
            pq.push(State(my_swap_array[queue_beg[queue_id]], queue_id));
        }
    }

    double step8 = MPI_Wtime();
    if (my_rank == 0) {
        printf("fnished step 8 local sort. elapse %lf\n", step8-begin_time);
    }

    for (int i = 0; i < my_recv_count; ++i) {
        printf("[%d] sorted %lu\n", my_rank, result[i]);
    }

    free(all_count_array);

    free(queue_beg);

    free(my_swap_array);

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