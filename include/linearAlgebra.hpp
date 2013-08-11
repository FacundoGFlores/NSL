/*****************************************************************************
Name		: 	linearAlgebra.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This library has got many of the basic linear algebra operations.
				It uses the linearMatrix class.
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

#ifndef LINEARALGEBRA_HPP_
#define LINEARALGEBRA_HPP_

#include "linearMatrix.hpp"
#include <cmath>

/**
 * Matrix - vector multiplication. Remember A and B
 * have to Height > 1. c = A*b
 * @param A the matrix
 * @param b the vector
 * @return a linearMatrix with Height = 1
 */
linearMatrix lAmatrixvecmult(const linearMatrix &A,
							 const linearMatrix &b);

/**
 * A simple matrix multiplication. C = A * B
 * @param A the matrix A
 * @param B the Matrix B
 * @return a matrix(a_width, b_height)
 */
linearMatrix lAmatrixmatrixmult(const linearMatrix &A,
								const linearMatrix &B);

/**
 * A simple dot product. alpha = <a,b>
 * @param a vector a
 * @param b vector b
 * @return a scalar
 */
float lAdotproduct(const linearMatrix &a,
						  const linearMatrix &b);

/**
 * Scalar-vector multiplication. c = a * beta
 * @param a vector a
 * @param b beta(scalar)
 * @return a vector
 */
linearMatrix lAscalarvecmult(const linearMatrix &a,
							 const float b);

/**
 * Vector addition. c = a + b
 * @param a vector a
 * @param b vector b
 * @return vector
 */
linearMatrix lAvectorvectorsum(const linearMatrix &a,
							   const linearMatrix &b);

/**
 * Matrix addition. C = A + B
 * @param A matrix A
 * @param B matrix B
 * @return matrix C
 */
linearMatrix lAmatrixmatrixsum(const linearMatrix &A,
							   const linearMatrix &B);
/**
 * Vector addition. c = a - b
 * @param a vector a
 * @param b vector b
 * @return vector
 */
linearMatrix lAvectorvectorsub(const linearMatrix &a,
							   const linearMatrix &b);

/**
 * Matrix addition. C = A - B
 * @param A matrix A
 * @param B matrix B
 * @return matrix C
 */
linearMatrix lAmatrixmatrixsub(const linearMatrix &A,
							   const linearMatrix &B);

/**
 * Returns the infinity norm of a vector
 * @param a the vector
 * @return the infinity norm
 */
float lAinfinitynorm(const linearMatrix &a);

#endif /* LINEARALGEBRA_HPP_ */
