/*****************************************************************************
Name		: 	linearSolver.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the linearSolver class. It allows us to solve linear
				system of equations( b = A*x ).
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

#ifndef LINEARSOLVER_HPP_
#define LINEARSOLVER_HPP_

#include "linearMatrix.hpp"

class linearSolver{
public:
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
	linearSolver(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, const linearMatrix &Vector_X,
				 const float Error, const int maxIter);

	/**
	 * A copy constructor. It copies the attributes of
	 * the object that inside its.
	 * @param that a linearSolver object
	 */
	linearSolver(const linearSolver &that);

	/**
	 * This is part of the rule of three in C++.
	 * We have to allow the user to use the assignment
	 * operator. linearSolver A = linearSolver that
	 * @param that a linearSolver object
	 * @return a linearSolver object
	 */
	linearSolver &operator =(const linearSolver &that);

	/**
	 * We set the allowed error in the numerical solver
	 * @param Error the allowed error
	 */
	void Set_Error(const float Error);

	/**
	 * We set the maximum
	 * @param MaxIter the maximum number of iterations
	 */
	void Set_MaximumIter(const int MaxIter);

	/**
	 * We get the allowed error
	 * @return the allowed
	 */
	float Get_Error(void);

	/**
	 * We get the allowed number of iterations
	 * @return allowed number of iterations
	 */
	int Get_MaximunIter(void);


protected:
	linearMatrix myMatrixA, myVectorB, myVectorX;
	float myError;
	int myMaximumIter;
};




#endif /* LINEARSOLVER_HPP_ */
