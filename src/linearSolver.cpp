/*****************************************************************************
Name		: 	linearSolver.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the implementation of the linearSolver class.
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

#include <linearSolver.hpp>

/**
 * This constructor creates a linearSolver. It will be used with
 * the nonParallel and cudaParallel solvers, they are
 * linearSolvers.
 * @param Matrix_A the matrix A
 * @param Vector_B the vector B
 * @param Vector_X the vector X
 * @param Error the allowed error
 * @param maxIter the allowed iterations
 */
linearSolver::linearSolver(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, const linearMatrix &Vector_X,
						   const float Error, const int maxIter)
{
	myMatrixA = Matrix_A;
	myVectorB = Vector_B;
	myVectorX = Vector_X;
	myError = Error;
	myMaximumIter = maxIter;
}

/**
 * A copy constructor. It copies the attributes of
 * the object that inside its.
 * @param that a linearSolver object
 */
linearSolver::linearSolver(const linearSolver &that)
{
	myMatrixA = that.myMatrixA;
	myVectorB = that.myVectorB;
	myVectorX = that.myVectorX;
	myError = that.myError;
	myMaximumIter = that.myMaximumIter;
}

/**
 * This is part of the rule of three in C++.
 * We have to allow the user to use the assignment
 * operator. linearSolver A = linearSolver that
 * @param that a linearSolver object
 * @return a linearSolver object
 */
linearSolver &linearSolver::operator =(const linearSolver &that)
{
	if(this != &that)
	{
		myMatrixA = that.myMatrixA;
		myVectorB = that.myVectorB;
		myVectorX = that.myVectorX;
		myError = that.myError;
		myMaximumIter = that.myMaximumIter;
	}
	return *this;
}

/**
 * We set the allowed error in the numerical solver
 * @param Error the allowed error
 */
void linearSolver::Set_Error(const float Error)
{
	myError = Error;
}

/**
 * We set the maximum
 * @param MaxIter the maximum number of iterations
 */
void linearSolver::Set_MaximumIter(const int MaxIter)
{
	myMaximumIter = MaxIter;
}

/**
 * We get the allowed error
 * @return the allowed
 */
float linearSolver::Get_Error(void)
{
	return myError;
}

/**
 * We get the allowed number of iterations
 * @return allowed number of iterations
 */
int linearSolver::Get_MaximunIter(void)
{
	return myMaximumIter;
}

