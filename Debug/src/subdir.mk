################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/bicg.cpp \
../src/bicgstab.cpp \
../src/cg.cpp \
../src/cholesky.cpp \
../src/cudaParallel.cpp \
../src/jacobi.cpp \
../src/linearAlgebra.cpp \
../src/linearMatrix.cpp \
../src/linearSolver.cpp \
../src/main.cpp \
../src/nonParallel.cpp \
../src/pcg.cpp \
../src/testing.cpp 

CU_SRCS += \
../src/cudabicg.cu \
../src/cudabicgstab.cu \
../src/cudacg.cu \
../src/cudajacobi.cu \
../src/cudapcg.cu 

CU_DEPS += \
./src/cudabicg.d \
./src/cudabicgstab.d \
./src/cudacg.d \
./src/cudajacobi.d \
./src/cudapcg.d 

OBJS += \
./src/bicg.o \
./src/bicgstab.o \
./src/cg.o \
./src/cholesky.o \
./src/cudaParallel.o \
./src/cudabicg.o \
./src/cudabicgstab.o \
./src/cudacg.o \
./src/cudajacobi.o \
./src/cudapcg.o \
./src/jacobi.o \
./src/linearAlgebra.o \
./src/linearMatrix.o \
./src/linearSolver.o \
./src/main.o \
./src/nonParallel.o \
./src/pcg.o \
./src/testing.o 

CPP_DEPS += \
./src/bicg.d \
./src/bicgstab.d \
./src/cg.d \
./src/cholesky.d \
./src/cudaParallel.d \
./src/jacobi.d \
./src/linearAlgebra.d \
./src/linearMatrix.d \
./src/linearSolver.d \
./src/main.d \
./src/nonParallel.d \
./src/pcg.d \
./src/testing.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-5.5/bin/nvcc -I/home/facu/cuda-workspace/NSL/include -G -g -O0 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=sm_21 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-5.5/bin/nvcc -I/home/facu/cuda-workspace/NSL/include -G -g -O0 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-5.5/bin/nvcc -I/home/facu/cuda-workspace/NSL/include -G -g -O0 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=sm_21 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-5.5/bin/nvcc --compile -G -I/home/facu/cuda-workspace/NSL/include -O0 -g -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_21  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


