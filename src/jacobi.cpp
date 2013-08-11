/*****************************************************************************
Name		: 	jacobi.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the Jacobi method implemented on C++.
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

#include <cmath>

#include <vector>

#include <sys/time.h>

#include <linearAlgebra.hpp>

using namespace std;

/**
 * This is the "( b-Ax ) / A_ii" part of the
 * jacobi algorithm. We need a diagonal vector
 * from the Matrix A, and we have to take a
 * LU matrix from Matrix A too.
 * @param X b - Ax
 * @param diagonal ( b-Ax ) / A_ii
 */
void jacobi_diagmult(linearMatrix &X, const vector<float> &diagonal)
{
	vector<float> X_Vec = X.Get_Vector();
	for(int i = 0; i < diagonal.size(); i++)
		if(diagonal[i] != 0)
			X_Vec[i] = X_Vec[i] / diagonal[i];
		else
		{
			cout << "Zero division" << endl;
			exit(1);
		}
	X.Set_Matrix(X_Vec, X.Get_Width(), X.Get_Height());
}

/**
 * Solve a linear system using the Jacobi method
 * @param Matrix_A
 * @param Vector_B
 * @param Vector_X
 * @param e
 * @param MaxIter
 */
long npjacobi(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter)
{
    timeval start;
    timeval end;

    long seconds, useconds, final;


	int N = Vector_B.Get_Width();
	float myerror = 1;

	vector<float> Diagonal = Matrix_A.Get_Diagonal();

	linearMatrix LU(Matrix_A.Get_Width(), Matrix_A.Get_Height());
	linearMatrix X_Old(N, 1);
	linearMatrix r(N, 1);

	//Take the LU matrix from Matrix A
	LU = Matrix_A.Get_LU();

	//Copy the Vector X
	X_Old = Vector_X;

	//Residual vector
	r.Generate_Null_Matrix();

    gettimeofday(&start, 0);

	int k = 0;


	while(myerror > e && k < MaxIter)
	{
		k++;
		Vector_X = lAmatrixvecmult(LU, X_Old);

		Vector_X = lAvectorvectorsub(Vector_B, Vector_X);

		jacobi_diagmult(Vector_X, Diagonal);

		r = lAvectorvectorsub(X_Old, Vector_X);

		myerror = lAinfinitynorm(r);
		X_Old = Vector_X;
	}

	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;

	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	return final;
}
