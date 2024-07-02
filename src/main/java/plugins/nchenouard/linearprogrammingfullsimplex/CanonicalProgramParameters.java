/*
 * Copyright (c) 2010-2024. Institut Pasteur.
 *
 * This file is part of Icy.
 * Icy is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Icy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Icy. If not, see <https://www.gnu.org/licenses/>.
 */

package plugins.nchenouard.linearprogrammingfullsimplex;

/**
 * Parameters used for the simplex algorithm
 * <p>
 * By definition, the optimization problem is of the form:
 *
 * <pre>
 * min_x c'*x or max_x c'*x
 * A*x &lt;= b
 * with x &gt;=0
 * </pre>
 * <p>
 * where x, b, c are column vectors and c' is the transpose of c. * stands for
 * the matrix multiplication.
 * <p>
 * Here each constraint A[i]*x can be either "equal" or "lower than or equal"
 * constraints based on the boolean value in the equalityConstraints vector:
 * <p>
 * A[i]*x == b if equalityConstraints[i]
 * <p>
 * A[i]*x &lt;= b if equalityConstraints[i]
 * <p>
 * <p>
 * <p>
 * Part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org
 *
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */
public class CanonicalProgramParameters {
    /**
     * The objective function coefficient: c'*x or sum_i c[i]*x[i] for the variable x.
     */
    public double[] c;
    /**
     * Defines whether the problem is a maximization or a minimization with respect to x.
     */
    public boolean maximization;

    /**
     * The constraint matrix structure. Each line defines a linear combination of the input variable.
     */
    public double[][] A;
    /**
     * The value for each constraint.
     */
    public double[] b;

    /**
     * Defines the type of constraint
     * <p>
     * A[i]*x == b if equalityConstraints[i]
     * <p>
     * A[i]*x &lt;= b if equalityConstraints[i]
     */
    public boolean[] equalityConstraints;
}