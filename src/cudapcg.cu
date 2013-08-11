/*****************************************************************************
Name		: 	cudapcg.cu
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the imeplementation of the pre conditionated
				conjugate gradient method using CUDA and cuBlas.
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

#include <cudapcg.cuh>

using namespace std;

/**
 * Set a null vector.
 * Remember you have to allocate memory previously
 * @param dst the null vector
 * @param N the vector Width
 */
void cudapcgFillZeros(float *dst, const int N)
{
	for(int i = 0; i < N; i++)
		dst[i] = 0;
}

/**
 * Take the vector from the linear container
 * @param dst the C++ vector
 * @param vec the C array
 */
void cudapcgGetArray(float *dst, const vector<float> &vec)
{
	for(int i = 0; i < vec.size(); i++)
		dst[i] = vec[i];
}


/**
 * Solve a linear system of equations using the
 * conjugate gradient method.
 * @param Matrix_A Matrix A
 * @param Vector_B Vector B
 * @param Vector_X Vector X
 * @param e allowed error
 * @param MaxIter allowed iterations
 * @return the time spent
 */
long cuda_PCG(const linearMatrix &Matrix_A, const linearMatrix &Matrix_M,
				const linearMatrix &Vector_B, linearMatrix &Vector_X,
				const float e, const int MaxIter)
{
	int k = 0;

	int N = Vector_B.Get_Width();

	float *nullvec = (float *)malloc(sizeof(float) * N);
	if(!nullvec)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	cudapcgFillZeros(nullvec, N);

	/* ********** HOST VARIABLES ********** */
	float *mA = (float *)malloc(sizeof(float) * N * N);
	if(!mA)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	float *mM = (float *)malloc(sizeof(float) * N * N);
	if(!mM)
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


	//Set vectors
	cudapcgGetArray(mA, Matrix_A.Get_Vector());
	cudapcgGetArray(vB, Vector_B.Get_Vector());
	cudapcgGetArray(vX, Vector_X.Get_Vector());
	cudapcgGetArray(mM, Matrix_M.Get_Vector());

	/* ************************************ */

	/* ********** DEVICE VARIABLES ********** */

	cudaError_t cudaStat ;
	cublasStatus_t stat ;
	cublasHandle_t handle ;

	float *dev_A;
	float *dev_M;
	float *dev_z;
	float *dev_B;
	float *dev_X;
	float *dev_tmp;
	float *dev_r;
	float *dev_p;
	float *dev_w;
	float *dev_r_tmp;
	float *dev_null;

	/* ======================== Allocating ======================== */
	cudaStat = cudaMalloc( (void **)& dev_null, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_A, N * N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_M, N * N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_z, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_tmp, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_p, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_w, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_r_tmp, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	/* ================================================================ */

	/* ==================== Settings ==================== */
	stat = cublasCreate(&handle);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS Initialization failed!" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetMatrix(N, N, sizeof(float), mA, N, dev_A, N);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting matrix failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetMatrix(N, N, sizeof(float), mM, N, dev_M, N);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting matrix failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_null, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector NULL failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_z, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector Z failed" << endl;
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
		cout << "CUBLAS setting vector tmp failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_r, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector r failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_p, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector p failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_w, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector w failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_r_tmp, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector r tmp failed" << endl;
		exit(EXIT_FAILURE);
	}
	/* ======================================================== */
	/* ************************************** */

	/* ********** CG Method ********** */

	timeval start;
	timeval end;

	long seconds, useconds, final;

	gettimeofday(&start, 0);

	float delta;
	float delta_old;
	float myerror;
	float beta = 0;
	float alpha;
	float norm_b;
	float new_e;
	float dotprod;

	//dev_tmp = - dev_A * dev_X(matrix-vector multiplication
	float gemv_alpha = -1;
	float gemv_beta = 1;
	cublasSgemv(handle, CUBLAS_OP_N, N, N, &gemv_alpha, dev_A, N, dev_X, 1, &gemv_beta, dev_tmp, 1);

	//dev_tmp =  dev_B + dev_tmp
	gemv_alpha = 1;
	cublasSaxpy(handle, N, &gemv_alpha, dev_B, 1, dev_tmp, 1);

	//r = dev_tmp
	cudaMemcpy(dev_r, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

	//delta = <r,r>
	cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotprod);
	delta = dotprod;
	myerror = sqrt(delta);

	//||Vector_B||
	cublasSdot(handle, N, dev_B, 1, dev_B, 1, &dotprod);
	norm_b = sqrt(dotprod);

	new_e = e * norm_b;


	while(myerror > new_e  && k < MaxIter)
	{
		k++;

		gemv_alpha = 1;
		gemv_beta = 1;
		cudaMemcpy(dev_z, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSgemv(handle, CUBLAS_OP_T, N, N, &gemv_alpha, dev_M, N, dev_r, 1, &gemv_beta, dev_z, 1);

		cublasSdot(handle, N, dev_z, 1, dev_r, 1, &delta);

		if(k == 1)
		{
			beta = 0;
			cudaMemcpy(dev_p, dev_z, N * sizeof(float), cudaMemcpyDeviceToDevice);
		}
		else
		{

			beta = delta / delta_old;
			//helper
			cudaMemcpy(dev_r_tmp, dev_z, N * sizeof(float), cudaMemcpyDeviceToDevice);

			//p = z + beta * p
			cublasSaxpy(handle, N, &beta, dev_p, 1, dev_r_tmp, 1);
			cudaMemcpy(dev_p, dev_r_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);
		}


		cudaMemcpy(dev_w, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		//w = A*p
		gemv_alpha = 1; // tmp
		gemv_beta = 1;
		cublasSgemv(handle, CUBLAS_OP_T, N, N, &gemv_alpha, dev_A, N, dev_p, 1, &gemv_beta, dev_w, 1);

		//alpha = delta / <p,w>
		cublasSdot(handle, N, dev_p, 1, dev_w, 1, &dotprod);
		alpha = delta / dotprod;

		//x = alpha * p + x
		gemv_alpha = alpha;
		cublasSaxpy(handle, N, &gemv_alpha, dev_p, 1, dev_X, 1);

		//r = r - alpha * w
		gemv_alpha = alpha * (-1);
		cublasSaxpy(handle, N, &gemv_alpha, dev_w, 1, dev_r, 1);

		delta_old = delta;

		//delta = <r,r>
		cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotprod);
		delta = dotprod;
		myerror = sqrt(delta);
	}
	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	/* ****************************** */

	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	//Give the result to host
	cudaMemcpy(vX, dev_X, N * sizeof(float), cudaMemcpyDeviceToHost);
	vector<float> v_tmp(vX, vX + N);
	Vector_X.Set_Matrix(v_tmp, N, 1);

	cudaFree(dev_M);
	cudaFree(dev_r_tmp);
	cudaFree(dev_w);
	cudaFree(dev_p);
	cudaFree(dev_r);
	cudaFree(dev_z);
	cudaFree(dev_tmp);
	cudaFree(dev_X);
	cudaFree(dev_B);
	cudaFree(dev_A);
	cudaFree(dev_null);

	cublasDestroy(handle);

	free(mA);
	free(vB);
	free(vX);
	free(nullvec);
	free(mM);

	return final;
}

