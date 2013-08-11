/*****************************************************************************
Name		: 	linearAlgebra.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the implementation of the linearAlgebra library.
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
#include <vector>

#include <linearAlgebra.hpp>

using namespace std;

/**
 * Matrix - vector multiplication. Remember A and B
 * have to Height > 1. c = A*b
 * @param A the matrix
 * @param b the vector
 * @return a linearMatrix with Height = 1
 */
linearMatrix lAmatrixvecmult(const linearMatrix &A,
							 const linearMatrix &b)
{
	//Take sizes
	int W_A = A.Get_Width();
	int H_A = A.Get_Height();

	int H_B = b.Get_Height();
	int W_B = b.Get_Width();

	//Take vectors
	vector<float> vectorA = A.Get_Vector();
	vector<float> vectorB = b.Get_Vector();
	vector<float> vectorC;

	vectorC.resize(W_A * W_B);

	//Check sizes
	if(H_A == W_B && H_B == 1)
	{
		//Perform it
		linearMatrix C(W_A, 1);
		for(int i = 0; i < W_A; i++)
		{
			float sum = 0;
			for(int j = 0; j < W_B; j++)
				sum += vectorA[i * W_B + j] * vectorB[j];
			vectorC[i] = sum;
		}
		C.Set_Matrix(vectorC, W_A, 1);
		return C;
	}
	else
	{
		cout << "Different sizes! Cannot perform matrixmult" << endl;
		exit(1);
	}
}

/**
 * A simple matrix multiplication. C = A * B
 * @param A the matrix A
 * @param B the Matrix B
 * @return a matrix(a_width, b_height)
 */
linearMatrix lAmatrixmatrixmult(const linearMatrix &A,
							 const linearMatrix &B)
{
	//Take sizes
	int W_A = A.Get_Width();
	int H_A = A.Get_Height();

	int H_B = B.Get_Height();
	int W_B = B.Get_Width();

	//Take vectors
	vector<float> vectorA = A.Get_Vector();
	vector<float> vectorB = B.Get_Vector();
	vector<float> vectorC;

	vectorC.resize(W_A * H_B);

	//Check sizes
	if(H_A == W_B)
	{
		//Perform it
		linearMatrix C(W_A, H_B);
		for(int i = 0; i < W_A; i++)
		{
			for(int j = 0; j < H_B; j++)
			{
				float sum = 0;
				for(int k = 0; k < W_B; k++)
				{
					float a = vectorA[i * H_A + k];
					float b = vectorB[k * W_B + j];
					sum += a * b;
				}
				vectorC[i * W_A + j] = sum;
			}
		}
		C.Set_Matrix(vectorC, W_A, W_B);
		return C;
	}
	else
	{
		cout << "Different sizes! Cannot perform matrixmult" << endl;
		exit(1);
	}
}

/**
 * A simple dot product. alpha = <a,b>
 * @param a vector a
 * @param b vector b
 * @return a scalar
 */
float lAdotproduct(const linearMatrix &a,
						  const linearMatrix &b)
{
	//Take sizes
	int W_A = a.Get_Width();
	int H_A = a.Get_Height();

	int H_B = b.Get_Height();
	int W_B = b.Get_Width();

	//Check if they are 1D
	if(H_A != 1 || H_B != 1)
	{
		cout << "Non 1D vector" << endl;
		exit(1);
	}

	//Check if they are equals
	if(W_A != W_B)
	{
		cout << "They are not equals" << endl;
		exit(1);
	}

	//Assume we can perform it
	float sum = 0;
	vector<float> Vector_A = a.Get_Vector();
	vector<float> Vector_B = b.Get_Vector();
	for(int i = 0; i < W_A; i++)
		sum += Vector_A[i] * Vector_B[i];
	return sum;
}

/**
 * Scalar-vector multiplication. c = a * beta
 * @param a vector a
 * @param b beta(scalar)
 * @return a vector
 */
linearMatrix lAscalarvecmult(const linearMatrix &a,
							 const float b)
{
	//Check if it is a 1D vector
	if(a.Get_Height() != 1)
	{
		cout << "It isn't a 1D vector" << endl;
		exit(1);
	}

	vector<float> A = a.Get_Vector();
	int N = a.Get_Width();
	for(int i = 0; i < N; i++)
		A[i] = A[i] * b;

	linearMatrix myRet(N, 1);
	myRet.Set_Matrix(A, N, 1);

	return myRet;
}

/**
 * Vector addition. c = a + b
 * @param a vector a
 * @param b vector b
 * @return vector
 */
linearMatrix lAvectorvectorsum(const linearMatrix &a,
							   const linearMatrix &b)
{
	//Take sizes
	int W_A = a.Get_Width();
	int H_A = a.Get_Height();

	int H_B = b.Get_Height();
	int W_B = b.Get_Width();

	//Check if they are 1D
	if(H_A != 1 || H_B != 1 )
	{
		cout << "They are not 1D vectors" << endl;
		exit(1);
	}

	//Check if they are equals
	if(W_A != W_B)
	{
		cout << "They are not equals" << endl;
		exit(1);
	}

	vector<float> A = a.Get_Vector();
	vector<float> B = b.Get_Vector();
	vector<float> C;

	int N = a.Get_Width();

	C.resize(N);

	for(int i = 0; i < N; i++)
		C[i] = A[i] + B[i];

	//Prepare return
	linearMatrix myRet(W_A, 1);
	myRet.Set_Matrix(C, W_A, 1);
	return myRet;
}

/**
 * Matrix addition. C = A + B
 * @param A matrix A
 * @param B matrix B
 * @return matrix C
 */
linearMatrix lAmatrixmatrixsum(const linearMatrix &A,
							   const linearMatrix &B)
{
	//Take sizes;
	int W_A = A.Get_Width();
	int H_A = A.Get_Height();

	int H_B = B.Get_Height();
	int W_B = B.Get_Width();

	//Check sizes
	if(W_A != W_B || H_A != H_B)
	{
		cout << "Can't perform it. Bad sizes!" << endl;
		exit(1);
	}

	vector<float> Matrix_A = A.Get_Vector();
	vector<float> Matrix_B = B.Get_Vector();
	vector<float> Matrix_C;

	Matrix_C.resize(W_A * H_A);

	for(int i = 0; i < W_A; i++)
		for(int j = 0; j < H_A; j++)
			Matrix_C[i * H_A + j] = Matrix_A[i * H_A + j] + Matrix_A[i * H_A + j];

	//Prepare return
	linearMatrix myRet(W_A, H_A);
	myRet.Set_Matrix(Matrix_C, W_A, H_A);
	return myRet;
}

/**
 * Vector addition. c = a - b
 * @param a vector a
 * @param b vector b
 * @return vector
 */
linearMatrix lAvectorvectorsub(const linearMatrix &a,
							   const linearMatrix &b)
{
	//Take sizes
	int W_A = a.Get_Width();
	int H_A = a.Get_Height();

	int H_B = b.Get_Height();
	int W_B = b.Get_Width();

	//Check if they are 1D
	if(H_A != 1 || H_B != 1 )
	{
		cout << "They are not 1D vectors" << endl;
		exit(1);
	}

	//Check if they are equals
	if(W_A != W_B)
	{
		cout << "They are not equals" << endl;
		exit(1);
	}

	vector<float> A = a.Get_Vector();
	vector<float> B = b.Get_Vector();
	vector<float> C;

	int N = a.Get_Width();

	C.resize(N);

	for(int i = 0; i < N; i++)
		C[i] = A[i] - B[i];

	//Prepare return
	linearMatrix myRet(W_A, 1);
	myRet.Set_Matrix(C, W_A, 1);
	return myRet;
}

/**
 * Matrix addition. C = A - B
 * @param A matrix A
 * @param B matrix B
 * @return matrix C
 */
linearMatrix lAmatrixmatrixsub(const linearMatrix &A,
							   const linearMatrix &B)
{
	//Take sizes;
	int W_A = A.Get_Width();
	int H_A = A.Get_Height();

	int H_B = B.Get_Height();
	int W_B = B.Get_Width();

	//Check sizes
	if(W_A != W_B || H_A != H_B)
	{
		cout << "Can't perform it. Bad sizes!" << endl;
		exit(1);
	}

	vector<float> Matrix_A = A.Get_Vector();
	vector<float> Matrix_B = B.Get_Vector();
	vector<float> Matrix_C;

	Matrix_C.resize(W_A * H_A);

	for(int i = 0; i < W_A; i++)
		for(int j = 0; j < H_A; j++)
			Matrix_C[i * H_A + j] = Matrix_A[i * H_A + j] - Matrix_A[i * H_A + j];

	//Prepare return
	linearMatrix myRet(W_A, H_A);
	myRet.Set_Matrix(Matrix_C, W_A, H_A);

	return myRet;
}

/**
 * Returns the infinity norm of a vector
 * @param a the vector
 * @return the infinity norm
 */
float lAinfinitynorm(const linearMatrix &a)
{
	float max = 0;
	int N = a.Get_Width();

	vector<float> v = a.Get_Vector();

	for(int i = 0; i < N; i++)
	{
		if(max < fabs(v[i]))
			max = fabs(v[i]);
	}

	return max;
}
