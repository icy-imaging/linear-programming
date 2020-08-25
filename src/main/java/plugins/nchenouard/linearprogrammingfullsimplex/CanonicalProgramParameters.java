package plugins.nchenouard.linearprogrammingfullsimplex;

/**
 * Parameters used for the simplex algorithm
 * 
 * By definition, the optimization problem is of the form:
 * min_x c'*x or max_x c'*x
 * A*x <= b
 * with x >=0
 * where x, b, c are column vectors and c' is the transpose of c. * stands for the matrix multiplication.
 * Here each constraint A[i]*x can be either "equal" or "lower than or equal" constraints
 * based on the boolean value in the equalityConstraints vector:
 * A[i]*x == b if equalityConstraints[i]
 * A[i]*x <= b if equalityConstraints[i]
 * 
 * 
 * Part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org
 * 
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 * @version 1.0
 * @date 2014-04-21
 * @license gpl v3.0
 * 
 */

public class CanonicalProgramParameters
{
	/**
	 * The objective function coefficient: c'*x or sum_i c[i]*x[i] for the variable x.
	 * */
	public double[] c;
	/**
	 * Defines whether the problem is a maximization or a minimization with respect to x.
	 * */
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
	 * A[i]*x == b if equalityConstraints[i]
	 * A[i]*x <= b if equalityConstraints[i]
	 */
	public boolean[] equalityConstraints;
}