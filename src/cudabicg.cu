/*****************************************************************************
Name		: 	cudabicg.cu
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the imeplementation of the
				biconjugate gradient method using CUDA and cuBlas.
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

#include <cudacg.cuh>


using namespace std;

/**
 * Set a null vector.
 * Remember you have to allocate memory previously
 * @param dst the null vector
 * @param N the vector Width
 */
void cudabicgFillZeros(float *dst, const int N)
{
	for(int i = 0; i < N; i++)
		dst[i] = 0;
}

/**
 * Take the vector from the linear container
 * @param dst the C++ vector
 * @param vec the C array
 */
void cudabicgGetArray(float *dst, const vector<float> &vec)
{
	for(int i = 0; i < vec.size(); i++)
		dst[i] = vec[i];
}


/**
 * Solve a linear system of equations using the
 * biconjugate gradient method.
 * @param Matrix_A Matrix A
 * @param Vector_B Vector B
 * @param Vector_X Vector X
 * @param e allowed error
 * @param MaxIter allowed iterations
 * @return the time spent
 */
long cudaBicg(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
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

	cudabicgFillZeros(nullvec, N);

	/* ********** HOST VARIABLES ********** */
	float *mA = (float *)malloc(sizeof(float) * N * N);
	if(!mA)
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
	cudabicgGetArray(mA, Matrix_A.Get_Vector());
	cudabicgGetArray(vB, Vector_B.Get_Vector());
	cudabicgGetArray(vX, Vector_X.Get_Vector());

	/* ************************************ */

	/* ********** DEVICE VARIABLES ********** */

	cudaError_t cudaStat ;
	cublasStatus_t stat ;
	cublasHandle_t handle ;

	float *dev_A;
	float *dev_B;
	float *dev_X;
	float *dev_tmp;
	float *dev_r;
	float *dev_rhat;
	float *dev_p;
	float *dev_phat;
	float *dev_null;
	float *dev_v;

	/* ======================== Allocating ======================== */

	cudaStat = cudaMalloc( (void **)& dev_v, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_A, N * N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_rhat, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_phat, N * sizeof(float) );
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

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_null, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorB failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_v, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vectorB failed" << endl;
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

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_rhat, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector rhat failed" << endl;
		exit(EXIT_FAILURE);
	}
	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_p, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector p failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_phat, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector phat failed" << endl;
		exit(EXIT_FAILURE);
	}

	/* ======================================================== */
	/* ************************************** */

	/* ********** BICG Method ********** */

    timeval start;
    timeval end;

    long seconds, useconds, final;

    gettimeofday(&start, 0);

	float gemv_alpha = -1;
	float gemv_betha = 1;

	//r = b - Ax
	cublasSgemv(handle, CUBLAS_OP_N, N, N, &gemv_alpha, dev_A, N, dev_X, 1, &gemv_betha, dev_tmp, 1);
	gemv_alpha = 1;
	cublasSaxpy(handle, N, &gemv_alpha, dev_B, 1, dev_tmp, 1);
	cudaMemcpy(dev_r, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);
	//rhat = r
	cudaMemcpy(dev_rhat, dev_r, N * sizeof(float), cudaMemcpyDeviceToDevice);


	float rho_old = 1;
	float dotproduct;
	cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotproduct);
	float rho = dotproduct;

	float alpha;
	float betha;

	cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotproduct);
	float lasterror = sqrt(dotproduct);

	cublasSdot(handle, N, dev_B, 1, dev_B, 1, &dotproduct);

	float norm_b = sqrt(dotproduct);

	float new_e = e * norm_b;


	while(lasterror > new_e  && k < MaxIter)
	{
		k++;

		betha = rho / rho_old;

		//p = r + betha * p
		gemv_alpha = betha;
		cudaMemcpy(dev_tmp, dev_r, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSaxpy(handle, N, &gemv_alpha, dev_p, 1, dev_tmp, 1);
		cudaMemcpy(dev_p, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);


		//phat = rhat + betha * phat
		gemv_alpha = betha;
		cudaMemcpy(dev_tmp, dev_rhat, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSaxpy(handle, N, &gemv_alpha, dev_phat, 1, dev_tmp, 1);
		cudaMemcpy(dev_phat, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//v = A*p
		gemv_alpha = 1;
		gemv_betha = 1;
		cudaMemcpy(dev_tmp, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSgemv(handle, CUBLAS_OP_N, N, N, &gemv_alpha, dev_A, N, dev_p, 1, &gemv_betha, dev_tmp, 1);
		cudaMemcpy(dev_v, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);


		//alpha = rhok / <phat,v>

		cublasSdot(handle, N, dev_phat, 1, dev_v, 1, &dotproduct);
		alpha = rho / dotproduct;


		//x = x + alpha * p
		gemv_alpha = alpha;
		cublasSaxpy(handle, N, &gemv_alpha, dev_p, 1, dev_X, 1);

		//r = r - alpha * v
		gemv_alpha = (-1) * alpha;
		cublasSaxpy(handle, N, &gemv_alpha, dev_v, 1, dev_r, 1);

		//rhat = rhat - alpha * A^T * phat

		//a)tmp = A^t * phat
		gemv_alpha = 1;
		gemv_betha = 1;
		cudaMemcpy(dev_tmp, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSgemv(handle, CUBLAS_OP_T, N, N, &gemv_alpha, dev_A, N, dev_phat, 1, &gemv_betha, dev_tmp, 1);
		//b)rhat = rhat - alpha * tmp
		gemv_alpha = (-1) * alpha;
		cublasSaxpy(handle, N, &gemv_alpha, dev_tmp, 1, dev_rhat, 1);

		rho_old = rho;
		cublasSdot(handle, N, dev_rhat, 1, dev_r, 1, &dotproduct);
		rho = dotproduct;

		cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotproduct);
		lasterror = sqrt(dotproduct);
	}

	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;

	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	//Give the result to host
	cudaMemcpy(vX, dev_X, N * sizeof(float), cudaMemcpyDeviceToHost);
	vector<float> v_tmp(vX, vX + N);
	Vector_X.Set_Matrix(v_tmp, N, 1);

	free(mA);
	free(vB);
	free(vX);
	free(nullvec);

	cudaFree(dev_A);
	cudaFree(dev_B);
	cudaFree(dev_X);
	cudaFree(dev_tmp);
	cudaFree(dev_r);
	cudaFree(dev_rhat);
	cudaFree(dev_p);
	cudaFree(dev_phat);
	cudaFree(dev_null);
	cudaFree(dev_v);


	return final;
}

