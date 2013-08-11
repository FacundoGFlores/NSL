/*****************************************************************************
Name		: 	cg.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is implementation of the Cholesky's factorization.
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

#include <cmath>

#include <vector>

#include <sys/time.h>

#include <linearAlgebra.hpp>

#include <cholesky.hpp>

using namespace std;

/**
 * Returns the L matrix from a square Matrix A using
 * Cholesky's factorization method.
 * @param Matrix_A the Matrix A
 * @return Matrix_L the Matrix L
 */
long npCholesky(const linearMatrix &Matrix_A, linearMatrix &Matrix_L)
{
	int N = Matrix_A.Get_Width();

	linearMatrix My_Matrix_L(N, N);
	My_Matrix_L.Generate_Null_Matrix();

	vector<float> A = Matrix_A.Get_Vector();
	vector<float> L = My_Matrix_L.Get_Vector();

	for(int k = 0; k < N; k++)
	{
		float sum = 0;
		for(int s = 0; s < k; s++)
			sum += L[k * N + s] * L[k * N + s];
		L[k * N + k] = sqrt(A[k * N + k] - sum);

		for(int i = k + 1; i < N; i++)
		{
			sum = 0;
			for(int s = 0; s < k; s++)
				sum += L[i * N + s] * L[k * N + s];
			L[i * N + k] = (A[i * N + k] - sum) / L[k * N + k];
		}
	}

	My_Matrix_L.Set_Matrix(L, N, N);
	Matrix_L = My_Matrix_L;
}
