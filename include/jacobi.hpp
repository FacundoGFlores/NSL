/*****************************************************************************
Name		: 	jacobi.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the Jacobi method.
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

/**
 * Solve a linear system using the Jacobi method
 * @param Matrix_A
 * @param Vector_B
 * @param Vector_X
 * @param e
 * @param MaxIter
 */
long npjacobi(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter);
