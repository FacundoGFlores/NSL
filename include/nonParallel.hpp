/*****************************************************************************
Name		: 	nonParallel.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the nonParallel class. It has got the different
				numerical methods for solving linear equations.
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


#ifndef NONPARALLEL_HPP_
#define NONPARALLEL_HPP_

#include "linearSolver.hpp"
#include <cg.hpp>
#include <jacobi.hpp>
#include <pcg.hpp>
#include <bicgstab.hpp>
#include <bicg.hpp>

class nonParallel : public linearSolver{
public:

	/**
	 * The constructor.
	 * @param A Matrix A
	 * @param B Vector B
	 * @param X Vector X
	 * @param myError allowed error
	 * @param maxIter allowed max iterations
	 */
	nonParallel(const linearMatrix &A, const linearMatrix &B, const linearMatrix &X,
				const float myError, const int maxIter);

	/**
	 * A copy constructor. The constructor copy
	 * the attributes of the that object
	 * insides itself.
	 * @param that a nonParallel object
	 */
	nonParallel(const nonParallel &that);

	/**
	 * This is part of the rule of three. We have
	 * to allow the user using assignment operator.
	 * @param that a nonParallel object
	 * @return a nonParallel object
	 */
	nonParallel &operator =(const nonParallel &that);

	/**
	 * We get the solution vector
	 * @return a linearMatrix object
	 */
	linearMatrix Get_Solution() const;

	/**
	 * Solve the set linear system of equations
	 * using the conjugate gradient method.
	 */
	void With_Conjugate_Gradient();

	/**
	 * Solve the set linear system of equations
	 * using the bi-conjugate gradient stable method.
	 */
	void With_Bicgstab();

	/**
	 * Solve the set linear system of equations
	 * using the bi-conjugate gradient method.
	 */
	void With_BICG();

	/**
	 * Solve the set linear system of equations
	 * using the Jacobi method.
	 */
	void With_Jacobi();

	/**
	 * Solve the set linear system of equations
	 * using the pre conditionated conjugate gradient method
	 */
	void With_PCG(const linearMatrix &Conditioner);


	/**
	 * @return the time used by the executed method
	 */
	long Get_Time_Used() const;

private:

	linearMatrix mySolution;
	long myTime;

};



#endif /* NONPARALLEL_HPP_ */
