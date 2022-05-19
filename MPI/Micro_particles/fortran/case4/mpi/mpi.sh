#!/bin/bash
# test mpi by n cores
# run this by typing: ./mpi.sh ./test_mpi 6 12 24
# ./test_mpi is the name of mpi file, 6, 12, 24 are numbr of cores

echo "Start mpi testing, filename: $1"
ag=1
for i in "$@"
do
# escape first argument
    if [[ $ag -gt 1 ]]
    then
        echo "test $i cores and output to log_$i"
        mpirun -np $i $1 -> "log_$i"
        echo "sleep 10 seconds"
        sleep 10
    fi
    let "ag++"
done