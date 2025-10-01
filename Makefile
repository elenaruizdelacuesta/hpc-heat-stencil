# Compilers
CC = gcc
MPICC = mpicc
ARCH ?= native

CFLAGS = -Iinclude -Wall -march=$(ARCH)
MPI_FLAGS = -Wall -Iinclude -fopenmp -O3 -march=$(ARCH)

SRC_DIR = src
SERIAL_SRC = $(SRC_DIR)/stencil_serial.c
PARALLEL_SRC = $(SRC_DIR)/stencil_parallel.c

SERIAL_OUT = stencil_serial
PARALLEL_OUT = stencil_parallel

# Default parameters for local runs
OMP_THREADS ?= 1
X ?= 100
Y ?= 100
NITER ?= 5

# Targets
all: $(PARALLEL_OUT)
serial: $(SERIAL_OUT)
parallel: $(PARALLEL_OUT)

$(SERIAL_OUT): $(SERIAL_SRC)
	$(CC) $(CFLAGS) $(SERIAL_SRC) -o $(SERIAL_OUT)

$(PARALLEL_OUT): $(PARALLEL_SRC)
	$(MPICC) $(MPI_FLAGS) $(PARALLEL_SRC) -o $(PARALLEL_OUT)

# Run commands for local testing
run-serial: $(SERIAL_OUT)
	./$(SERIAL_OUT) -f 0 -e 4 -o 0 -p 1 -x $(X) -y $(Y) -n $(NITER)

run-parallel: $(PARALLEL_OUT)
	export OMP_NUM_THREADS=$(OMP_THREADS)
	export OMP_PROC_BIND=close
	export OMP_PLACES=cores
	./$(PARALLEL_OUT) -v 0 -o 0 -p 1 -x $(X) -y $(Y) -n $(NITER)
	
# Clean
clean:
	rm -f $(SERIAL_OUT) $(PARALLEL_OUT)
	rm -f output/*.bin
	rm -f output/*.csv
