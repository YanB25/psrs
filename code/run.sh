#!/bin/sh
path="/GPUFS/nsccgz_yfdu_16/ouyry/shared_dir/psrs_data"
p_array=(1 2 4 8 16 32 64 112)
n_array=(12 14 18 22 26 30 31)
for i in {0..7}
do
    for j in {0..6}
    do
        echo ${p_array[i]} ${n_array[j]}
        I_MPI_FABRICS=shm:ofa mpirun -n ${p_array[i]} -machine machinefile ./psrs.out $path ${n_array[j]} > out/p${p_array[i]}_n${n_array[j]}
    done
done

# 如果你发现代码跑不了，帮我顺手fix一下
# mpirun -n ${线程数} ./psrs.out "path/to/data" ${数据量的n的值}
# 即./psrs.out的第一个参数是输入文件名，第二个参数是2^n的那个n