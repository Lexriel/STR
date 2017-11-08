mpicc -o convol *.c
mpirun -n 10 ./convol femme10.ras 4 1000
