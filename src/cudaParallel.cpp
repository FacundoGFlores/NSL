/*****************************************************************************
Name		: 	cudaParallel.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the implementation of the Parallel class.
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
#include <sys/time.h>

#include <cudaParallel.hpp>
#include <cudacg.cuh>
#include <cudajacobi.cuh>
#include <cudapcg.cuh>
#include <cudabicgstab.cuh>
#include <cudabicg.cuh>

using namespace std;

/**
 * The constructor.
 * @param A Matrix A
 * @param B Vector B
 * @param X Vector X
 * @param myError allowed error
 * @param maxIter allowed max iterations
 */
Parallel::Parallel(const linearMatrix &A, const linearMatrix &B, const linearMatrix &X,
						const float myError, const int maxIter):
						linearSolver(A, B, X, myError, maxIter)
{
	mySolution = X;
}

/**
 * A copy constructor. The constructor copy
 * the attributes of the that object
 * insides itself.
 * @param that a Parallel object
 */
Parallel::Parallel(const Parallel &that): linearSolver(that.myMatrixA, that.myVectorB, that.myVectorX, that.myError, that.myMaximumIter)
{
	mySolution = that.myVectorX;
}


/**
 * This is part of the rule of three. We have
 * to allow the user using assignment operator.
 * @param that a Parallel object
 * @return a Parallel object
 */
Parallel &Parallel::operator =(const Parallel &that)
{
	if(this != &that)
	{
		linearSolver(that.myMatrixA, that.myVectorB,
						that.myVectorX, that.myError, that.myMaximumIter);
		mySolution = that.myVectorX;
	}
	return *this;
}

/**
 * We get the solution vector
 * @return a linearMatrix object
 */
linearMatrix Parallel::Get_Solution() const
{
	return mySolution;
}

/**
 * Solve the set linear system of equations
 * using the conjugate gradient method.
 */
void Parallel::With_Conjugate_Gradient()
{
	linearMatrix vector_copy = myVectorX;



	myTime = cudaConjugate_Gradient(myMatrixA, myVectorB,  myVectorX,
						myError, myMaximumIter);

	mySolution = myVectorX;
	myVectorX = vector_copy;
}

/**
 * Solve the set linear system of equations
 * using the Jacobi method
 */
void Parallel::With_Jacobi()
{
	linearMatrix vector_copy = myVectorX;



	myTime = cudajacobi(myMatrixA, myVectorB,  myVectorX,
						myError, myMaximumIter);

	mySolution = myVectorX;
	myVectorX = vector_copy;
}

/**
 * Solve the set linear system of equations
 * using the pre conditionated conjugate gradient method
 */
void Parallel::With_PCG(const linearMatrix &Conditioner)
{
	linearMatrix vector_copy = myVectorX;



	myTime = cuda_PCG(myMatrixA, Conditioner, myVectorB,  myVectorX,
						myError, myMaximumIter);

	mySolution = myVectorX;
	myVectorX = vector_copy;
}

/**
 * Solve the set linear system of equations
 * using the bi-conjugate gradient method
 */
void Parallel::With_BICG()
{
	linearMatrix vector_copy = myVectorX;

	myTime = cudaBicg(myMatrixA, myVectorB,  myVectorX,
						myError, myMaximumIter);

	mySolution = myVectorX;
	myVectorX = vector_copy;
}

/**
 * Solve the set linear system of equations
 * using the bi-conjugate gradient stable method
 */
void Parallel::With_BICGSTAB()
{
	linearMatrix vector_copy = myVectorX;

	myTime = cudaBicgstab(myMatrixA, myVectorB,  myVectorX,
						myError, myMaximumIter);

	mySolution = myVectorX;
	myVectorX = vector_copy;
}

/**
 * @return the time used by the executed method
 */
long Parallel::Get_Time_Used() const
{
	return myTime;
}
