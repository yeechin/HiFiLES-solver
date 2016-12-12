# \file Makefile
# \brief _____________________________
# \author - Original code: SD++ developed by Patrice Castonguay, Antony Jameson,
#                          Peter Vincent, David Williams (alphabetical by surname).
#         - Current development: Aerospace Computing Laboratory (ACL) directed
#                                by Prof. Jameson. (Aero/Astro Dept. Stanford University).
# \version 1.0.0
#
# HiFiLES (High Fidelity Large Eddy Simulation).
# Copyright (C) 2013 Aerospace Computing Laboratory.

# Settings

# TODO: add config script to build parmetis and look for BLAS, parmetis, CUDA, MPI and TECIO libraries

# Copy the appropriate file from makefiles/makefile.machine.in to HIFILES_HOME as makefile.in
#
include makefile.OSX.in

# Compiler

ifeq ($(COMP),GCC)
	CC      = g++
endif

ifeq ($(COMP),INTEL)
	CC      = icpc
endif

ifeq ($(PARALLEL),MPI)
	CC      = mpicxx
endif

ifeq ($(NODE),GPU)
	NVCC	= nvcc
endif

# Pre-processing macros

OPTS    = -D_$(NODE)

ifeq ($(PARALLEL),MPI)
	OPTS    += -D_$(PARALLEL) 
endif

ifeq ($(MACHINE),ENRICO)
	OPTS    += -D_$(MACHINE)
endif

ifeq ($(MACHINE),YOSEMITESAM)
	OPTS    += -D_$(MACHINE)
endif

# Includes

OPTS    += -I include 

ifeq ($(NODE),GPU)
	OPTS	+= -I $(CUDA_DIR)/include 
endif

ifeq ($(TECIO),YES)
	OPTS += -I $(TECIO_DIR)/tecsrc
endif

ifeq ($(PARALLEL),MPI)
	OPTS += -I $(MPI_DIR)/include
	OPTS += -I $(PARMETIS_DIR)/include
	OPTS += -I $(PARMETIS_DIR)/metis/include
endif

ifeq ($(CODE),DEBUG)
	OPTS	+= -g 
endif

ifeq ($(CODE),RELEASE)
	ifeq ($(COMP),GCC)
		OPTS	+= -O3 
	endif
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif

ifeq ($(CODE),RELEASE)
	ifeq ($(NODE),GPU)
		OPTS	+= -w
	endif
endif

ifeq ($(COMP),INTEL)
	ifeq ($(BLAS),MKL_BLAS)
		OPTS	+= -mkl=parallel
	endif
endif

# Libraries

ifeq ($(BLAS),ACCELERATE_BLAS)
	LIBS	+= -framework Accelerate
	OPTS	+= -flax-vector-conversions -D_$(BLAS)
endif

ifeq ($(BLAS),STANDARD_BLAS)
	LIBS    += -L $(BLAS_DIR)/lib -lcblas
	OPTS    += -I $(BLAS_DIR)/include -D_$(BLAS)
endif

ifeq ($(BLAS),ATLAS_BLAS)
	LIBS    += -L $(BLAS_DIR)/lib -lcblas -latlas
	OPTS    += -I $(BLAS_DIR)/include -D_STANDARD_BLAS
endif

ifeq ($(BLAS),NO_BLAS)
	OPTS	+= -D_$(BLAS)
endif

ifeq ($(NODE),GPU)
	LIBS	+= -L $(CUDA_DIR)/lib64 -lcudart -lcublas -lcusparse -lm
endif

# Source

SRC	= src/
OBJ	= obj/
BIN	= bin/
vpath %.cpp src
vpath %.cu src
vpath %.h include
vpath %.hpp include

# Objects

OBJS    = $(OBJ)HiFiLES.o $(OBJ)geometry.o $(OBJ)mesh.o $(OBJ)matrix_structure.o $(OBJ)vector_structure.o $(OBJ)linear_solvers_structure.o $(OBJ)solver.o $(OBJ)output.o $(OBJ)eles.o $(OBJ)eles_tris.o $(OBJ)eles_quads.o $(OBJ)eles_hexas.o $(OBJ)eles_tets.o $(OBJ)eles_pris.o $(OBJ)inters.o $(OBJ)int_inters.o $(OBJ)bdy_inters.o $(OBJ)funcs.o $(OBJ)flux.o $(OBJ)source.o $(OBJ)global.o $(OBJ)input.o $(OBJ)cubature_1d.o $(OBJ)cubature_tri.o $(OBJ)cubature_quad.o $(OBJ)cubature_hexa.o $(OBJ)cubature_tet.o

ifeq ($(NODE),GPU)
	OBJS	+=  $(OBJ)cuda_kernels.o
endif

ifeq ($(PARALLEL),MPI)
	OBJS += $(OBJ)mpi_inters.o
	OBJS += $(PARMETIS_BUILD_DIR)/libparmetis.a $(PARMETIS_BUILD_DIR)/libmetis.a
endif
	
ifeq ($(TECIO),YES)
	OBJS += $(TECIO_DIR)/tecio.a
endif

# Compile

.PHONY: default help clean

default: HiFiLES

help:
	@echo ' '
	@echo 'You should specify a target to make: make <arg> '
	@echo 'where <arg> is one of the following options: ' 
	@echo '	- HiFiLES :	compiles HiFiLES solver'
	@echo '	- clean :	clean HiFiLES'
	@echo ' '

HiFiLES: $(OBJS)
	$(CC) $(OPTS) -o $(BIN)HiFiLES $(OBJS) ${LIBS}

$(OBJ)HiFiLES.o: HiFiLES.cpp geometry.h input.h flux.h source.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)geometry.o: geometry.cpp geometry.h input.h  error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)mesh.o: mesh.cpp matrix_structure.hpp vector_structure.hpp linear_solvers_structure.hpp input.h  error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)matrix_structure.o: matrix_structure.cpp matrix_structure.hpp vector_structure.hpp linear_solvers_structure.hpp input.h  error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)vector_structure.o: vector_structure.cpp matrix_structure.hpp vector_structure.hpp linear_solvers_structure.hpp input.h  error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)linear_solvers_structure.o: linear_solvers_structure.cpp matrix_structure.hpp vector_structure.hpp linear_solvers_structure.hpp input.h  error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)solver.o: solver.cpp solver.h input.h  error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)output.o: output.cpp output.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)eles.o: eles.cpp eles.h arrayt.h error.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)eles_tris.o: eles_tris.cpp eles_tris.h eles.h funcs.h input.h arrayt.h arrayt.h cubature_1d.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)eles_quads.o: eles_quads.cpp eles_quads.h eles.h funcs.h input.h arrayt.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)eles_hexas.o: eles_hexas.cpp eles_hexas.h eles.h funcs.h input.h arrayt.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)eles_tets.o: eles_tets.cpp eles_tets.h eles.h funcs.h input.h arrayt.h error.h cubature_tri.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)eles_pris.o: eles_pris.cpp eles_pris.h eles.h funcs.h input.h arrayt.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)inters.o: inters.cpp inters.h flux.h funcs.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)int_inters.o: int_inters.cpp int_inters.h inters.h flux.h funcs.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)bdy_inters.o: bdy_inters.cpp bdy_inters.h inters.h flux.h funcs.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

ifeq ($(PARALLEL),MPI)
$(OBJ)mpi_inters.o: mpi_inters.cpp mpi_inters.h inters.h flux.h funcs.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
endif

$(OBJ)funcs.o: funcs.cpp funcs.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)cubature_1d.o: cubature_1d.cpp cubature_1d.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)cubature_tri.o: cubature_tri.cpp cubature_tri.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)cubature_quad.o: cubature_quad.cpp cubature_quad.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)cubature_hexa.o: cubature_hexa.cpp cubature_hexa.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)cubature_tet.o: cubature_tet.cpp cubature_tet.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)flux.o: flux.cpp flux.h arrayt.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)source.o: source.cpp source.h arrayt.h input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

$(OBJ)input.o: input.cpp input.h error.h
	$(CC) $(OPTS)  -c -o $@ $<
	
$(OBJ)global.o: global.cpp global.h error.h
	$(CC) $(OPTS)  -c -o $@ $<

ifeq ($(NODE),GPU)	
NVCC_FLAGS=
# For very verbose compilation:
#NVCC_FLAGS=--ptxas-options=-v
$(OBJ)cuda_kernels.o: cuda_kernels.cu cuda_kernels.h error.h util.h
	$(NVCC) -arch=sm_20 $(NVCC_FLAGS) -Xcudafe "--diag_suppress=subscript_out_of_range" $(OPTS) -c -o $@ $<
endif

clean: 
	rm $(BIN)HiFiLES $(OBJ)*.o
