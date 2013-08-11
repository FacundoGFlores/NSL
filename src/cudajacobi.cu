/*****************************************************************************
Name		: 	cudajacobi.cu
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the Jacobi method implemented on CUDA.
License     :
    			Copyright (C) 2012 Flores, Facundo Gabriel

    			This program is free software: you can redistribute it and/or modify
    			it under the terms of the GNU General Public License as published by
    			the Free Software Foundation, either version 3 of the License, or
    			(at your option) any later version.

    			This program is distributed in the hope that it will be useful,
    			but WITHOUT ANY WARRANTY; without even the implied warranty of
    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    			GNU General Public License for more details.

    			You should have received a copy of the GNU General Public License
    			along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************** */

#include <iostream>

#include <cstdlib>

#include <cuda_runtime.h>

#include <cublas_v2.h>

#include <vector>

#include <sys/time.h>

#include <cudajacobi.cuh>


using namespace std;

/**
 * This is the "( b-Ax ) / A_ii" part of the
 * jacobi algorithm. We need a diagonal vector
 * from the Matrix A, and we have to take a
 * LU matrix from Matrix A too.
 * @param diag Diagonal vector
 * @param vec  b - Ax
 * @param dim  Width of Diagonal
 */
__global__ void parallel_jacobi_diagmult(float * diag, float * vec, int dim)
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     while ( tid < dim)
     {
          if (diag[tid] != 0)
          {
               vec[tid] /= diag[tid];
          }
          tid += gridDim.x * blockDim.x;
     }

}

/**
 * Set a null vector.
 * Remember you have to allocate memory previously
 * @param dst the null vector
 * @param N the vector Width
 */
void cudajacFillZeros(float *dst, const int N)
{
	for(int i = 0; i < N; i++)
		dst[i] = 0;
}

/**
 * Take the vector from the linear container
 * @param dst the C++ vector
 * @param vec the C array
 */
void cudajacGetArray(float *dst, const vector<float> &vec)
{
	for(int i = 0; i < vec.size(); i++)
		dst[i] = vec[i];
}

/**
 * Solve a linear system using the Jacobi method
 * @param Matrix_A
 * @param Vector_B
 * @param Vector_X
 * @param e
 * @param MaxIter
 */
long cudajacobi(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter)
{
	int N = Matrix_A.Get_Width();
	int M = Matrix_A.Get_Height();
	float myerror = 1;

	linearMatrix LU(N, M);

	LU = Matrix_A.Get_LU();

	/* ================ HOST VARIABLES ================ */

	float *mLU = (float *)malloc(sizeof(float) * N * N);
	if(!mLU)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	float *vB = (float *)malloc(sizeof(float) * N);
	if(!vB)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	float *vX = (float *)malloc(sizeof(float) * N);
	if(!vX)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	float *vDiagonal = (float *)malloc(sizeof(float) * N);
	if(!vDiagonal)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	float *nullvec = (float *)malloc(sizeof(float) * N);
	if(!nullvec)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}


	//Set vectors
	cudajacGetArray(mLU, LU.Get_Vector());
	cudajacGetArray(vB, Vector_B.Get_Vector());
	cudajacGetArray(vX, Vector_X.Get_Vector());
	cudajacGetArray(vDiagonal, Matrix_A.Get_Diagonal());

	cudajacFillZeros(nullvec, N);


	/* ================================================ */

	/* ================ DEVICE VARIABLES ================ */
	cudaError_t cudaStat ;
	cublasStatus_t stat ;
	cublasHandle_t handle ;

	float *dev_LU;
	float *dev_B;
	float *dev_X;
	float *dev_X_Old;
	float *dev_r;
	float *dev_Diagonal;
	float *dev_tmp;
	float *dev_null;

	/* **************** ALLOCATING ********************** */
	cudaStat = cudaMalloc( (void **)& dev_LU, N * N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_tmp, N * N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_B, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_X, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_X_Old, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_Diagonal, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_r, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_null, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	/* ************************************************** */

	/* **************** SETTINGS ********************** */
	stat = cublasCreate(&handle);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS Initialization failed!" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetMatrix(N, N, sizeof(float), mLU, N, dev_LU, N);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting matrix failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), vB, 1, dev_B, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorB failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_X, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorX failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_tmp, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorX failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_null, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorX failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_X_Old, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorX failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_r, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector r failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), vDiagonal, 1, dev_Diagonal, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorX failed" << endl;
		exit(EXIT_FAILURE);
	}
	/* ************************************************ */

	/* ================================================== */

    int dimBlocks, dimThreads;
    dimThreads = 32;
    dimBlocks = (N / dimThreads) + 1;

	float gemv_alpha;
	float gemv_beta;

    timeval start, end;
    long seconds, useconds, final;


    gettimeofday(&start, 0);

    int k = 0;

	while(myerror > e && k < MaxIter)
	{
		k++;

		cudaMemcpy(dev_tmp, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		gemv_alpha = -1;
		gemv_beta = 1;
		cublasSgemv(handle, CUBLAS_OP_N, N, N, &gemv_alpha, dev_LU, N, dev_X, 1, &gemv_beta, dev_tmp, 1);
		cudaMemcpy(dev_X, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		gemv_alpha = 1;
		cublasSaxpy(handle, N, &gemv_alpha, dev_B, 1, dev_X, 1);

		parallel_jacobi_diagmult<<<dimBlocks, dimThreads>>>(dev_Diagonal, dev_X, N);


		gemv_alpha = -1;
		cudaMemcpy(dev_tmp, dev_X, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSaxpy(handle, N, &gemv_alpha, dev_X_Old, 1, dev_X, 1);
		cudaMemcpy(dev_r, dev_X, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemcpy(dev_X, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		cublasSasum(handle, N, dev_r, 1, &myerror);

		cudaMemcpy(dev_X_Old, dev_X, N * sizeof(float), cudaMemcpyDeviceToDevice);
	}


	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	cudaMemcpy(vX, dev_X, N * sizeof(float), cudaMemcpyDeviceToHost);
	vector<float> v_tmp(vX, vX + N);
	Vector_X.Set_Matrix(v_tmp, N, 1);

	cudaFree(dev_LU);
	cudaFree(dev_B);
	cudaFree(dev_X);
	cudaFree(dev_X_Old);
	cudaFree(dev_r);
	cudaFree(dev_Diagonal);
	cudaFree(dev_tmp);
	cudaFree(dev_null);

	free(mLU);
	free(vB);
	free(vX);
	free(vDiagonal);
	free(nullvec);

	return final;
}
