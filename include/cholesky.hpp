/*****************************************************************************
Name		: 	cg.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the Cholesky's factorization method.
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

#ifndef CHOLESKY_HPP_
#define CHOLESKY_HPP_

/**
 * Returns the L matrix from a square Matrix A using
 * Cholesky's factorization method.
 * @param Matrix_A the Matrix A
 * @return Matrix_L the Matrix L
 */
long npCholesky(const linearMatrix &Matrix_A, linearMatrix &Matrix_L);


#endif /* CHOLESKY_HPP_ */
