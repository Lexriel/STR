EXECS=mpi_hello_world
MPICC?=mpicc

all: ${EXECS}

mpi_hello_world: mpi_hello_world.c
	${MPICC} -o ${EXECS} mpi_hello_world.c

clean:
	rm ${EXECS}
