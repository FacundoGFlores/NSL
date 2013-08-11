/*****************************************************************************
Name		: 	bicgstab.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is implementation of the bi-conjugate gradient stable
				using cuBlas library
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

#include <cudabicgstab.cuh>

using namespace std;

void show(float *vec, int N)
{
	for(int i = 0; i < N; i++)
		cout << vec[i] << endl;
}

/**
 * Take the vector from the linear container
 * @param dst the C++ vector
 * @param vec the C array
 */
void cudabicgstabGetArray(float *dst, const vector<float> &vec)
{
	for(int i = 0; i < vec.size(); i++)
		dst[i] = vec[i];
}

void cudabicgstabFillZeros(float *dst, const int N)
{
	for(int i = 0; i < N; i++)
		dst[i] = 0;
}

/**
 * Solve a linear system of equations using the
 * bi-conjugate gradient stable method.
 * @param Matrix_A Matrix A
 * @param Vector_B Vector B
 * @param Vector_X Vector X
 * @param e allowed error
 * @param MaxIter allowed iterations
 * @return the time spent
 */
long cudaBicgstab(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter)
{

	/* ********** ALGORITHM HOST VARIABLES ********** */

	int N = Vector_B.Get_Width();

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

	float *nullvec = (float *)malloc(sizeof(float) * N);
	if(!nullvec)
	{
		cout << "Cannot allocate memory" << endl;
		exit(EXIT_FAILURE);
	}

	//Set vectors
	cudabicgstabGetArray(mA, Matrix_A.Get_Vector());
	cudabicgstabGetArray(vB, Vector_B.Get_Vector());
	cudabicgstabGetArray(vX, Vector_X.Get_Vector());

	cudabicgstabFillZeros(nullvec, N);
	/* ************************************ */

	/* ********** DEVICE VARIABLES ********** */

	cudaError_t cudaStat ;
	cublasStatus_t stat ;
	cublasHandle_t handle ;

	float *dev_A;
	float *dev_B;
	float *dev_X;
	float *dev_tmp;
	float *dev_tmp1;
	float *dev_tmp2;
	float *dev_r;
	float *dev_rhat0;
	float *dev_v;
	float *dev_p;
	float *dev_s;
	float *dev_t;
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

	cudaStat = cudaMalloc( (void **)& dev_tmp1, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_tmp2, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_rhat0, N * sizeof(float) );
	if(cudaStat != EXIT_SUCCESS)
	{
		cout << "device memory allocation failed" << endl;
		exit(EXIT_FAILURE);
	}

	cudaStat = cudaMalloc( (void **)& dev_v, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_s, N * sizeof(float) );
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

	cudaStat = cudaMalloc( (void **)& dev_t, N * sizeof(float) );
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

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_tmp1, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector tmp1 failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_tmp2, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector tmp2 failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_r, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector r failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_rhat0, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector r failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_v, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector v failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_p, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector p failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_s, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector s failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_r_tmp, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector r tmp failed" << endl;
		exit(EXIT_FAILURE);
	}

	stat = cublasSetVector(N, sizeof(float), nullvec, 1, dev_t, 1);
	if(stat != CUBLAS_STATUS_SUCCESS)
	{
		cout << "CUBLAS setting vector t failed" << endl;
		exit(EXIT_FAILURE);
	}

	/* ======================================================== */
	/* ************************************** */

	/* ********** BICGSTAB Method ********** */


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
	cudaMemcpy(dev_rhat0, dev_r, N * sizeof(float), cudaMemcpyDeviceToDevice);

	float betha;
	float alpha = 1;
	float w = 1;
	float rho_old = 1;
	float aux_norm;
	float dotprod;
	float lasterror;

	int k = 0;

	//rho = <r,r>
	cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotprod);
	float rho = dotprod;
	lasterror = sqrt(rho);

	//||Vector_B||
	cublasSdot(handle, N, dev_B, 1, dev_B, 1, &dotprod);
	float norm_b = sqrt(dotprod);
	float new_e = e * norm_b;



	while(lasterror > new_e  && k < MaxIter)
	{
		k++;

		betha = (rho / rho_old) * (alpha / w);

		//p = r + betha(p - w*v)
		gemv_betha = (-1) * w;
		cudaMemcpy(dev_tmp, dev_p, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSaxpy(handle, N, &gemv_betha, dev_v, 1, dev_tmp, 1);
		cudaMemcpy(dev_tmp1, dev_r, N * sizeof(float), cudaMemcpyDeviceToDevice);
		gemv_betha = betha;
		cublasSaxpy(handle, N, &gemv_betha, dev_tmp, 1, dev_tmp1, 1);
		cudaMemcpy(dev_p, dev_tmp1, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//v = A*p
		gemv_alpha = 1;
		gemv_betha = 1;
		cudaMemcpy(dev_tmp, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSgemv(handle, CUBLAS_OP_N, N, N, &gemv_alpha, dev_A, N, dev_p, 1, &gemv_betha, dev_tmp, 1);
		cudaMemcpy(dev_v, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//alpha = rho / <rhat0,v>
		cublasSdot(handle, N, dev_rhat0, 1, dev_v, 1, &dotprod);
		alpha = rho / dotprod;

		//s = r - alpha * v
		gemv_alpha = (-1) * alpha;
		cudaMemcpy(dev_tmp, dev_r, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSaxpy(handle, N, &gemv_alpha, dev_v, 1, dev_tmp, 1);
		cudaMemcpy(dev_s, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//t = A * s
		gemv_alpha = 1;
		gemv_betha = 1;
		cudaMemcpy(dev_tmp, dev_null, N * sizeof(float), cudaMemcpyDeviceToDevice);
		cublasSgemv(handle, CUBLAS_OP_N, N, N, &gemv_alpha, dev_A, N, dev_s, 1, &gemv_betha, dev_tmp, 1);
		cudaMemcpy(dev_t, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//w = <t,s>/sqrt(<t,t>)
		cublasSdot(handle, N, dev_t, 1, dev_t, 1, &dotprod);
		aux_norm = sqrt(dotprod);
		cublasSdot(handle, N, dev_t, 1, dev_s, 1, &dotprod);
		w = dotprod / aux_norm;

		//rho = -w*rhat0*t
		rho_old = rho;
		cublasSdot(handle, N, dev_rhat0, 1, dev_t, 1, &dotprod);
		rho = dotprod;
		rho = rho * (-1) * w;

		//x + alpha * p + w * s
		cudaMemcpy(dev_tmp1, dev_X, N * sizeof(float), cudaMemcpyDeviceToDevice);
		gemv_alpha = alpha;
		cublasSaxpy(handle, N, &gemv_alpha, dev_p, 1, dev_tmp1, 1);
		cudaMemcpy(dev_tmp2, dev_tmp1, N * sizeof(float), cudaMemcpyDeviceToDevice);
		gemv_alpha = w;
		cublasSaxpy(handle, N, &gemv_alpha, dev_s, 1, dev_tmp2, 1);
		cudaMemcpy(dev_X, dev_tmp2, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//r = s- w * t
		cudaMemcpy(dev_tmp, dev_s, N * sizeof(float), cudaMemcpyDeviceToDevice);
		gemv_alpha = (-1) * w;
		cublasSaxpy(handle, N, &gemv_alpha, dev_t, 1, dev_tmp, 1);
		cudaMemcpy(dev_r, dev_tmp, N * sizeof(float), cudaMemcpyDeviceToDevice);

		//lasterror = sqrt(<r,r>)
		cublasSdot(handle, N, dev_r, 1, dev_r, 1, &dotprod);
		lasterror = sqrt(dotprod);
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

	cudaFree(dev_A);
	cudaFree(dev_B);
	cudaFree(dev_X);
	cudaFree(dev_tmp);
	cudaFree(dev_tmp1);
	cudaFree(dev_tmp2);
	cudaFree(dev_r);
	cudaFree(dev_rhat0);
	cudaFree(dev_v);
	cudaFree(dev_p);
	cudaFree(dev_s);
	cudaFree(dev_t);
	cudaFree(dev_r_tmp);
	cudaFree(dev_null);



	free(mA);
	free(vB);
	free(vX);
	free(nullvec);

	return final;
}

