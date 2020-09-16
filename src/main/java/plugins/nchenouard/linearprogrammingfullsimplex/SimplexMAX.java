package plugins.nchenouard.linearprogrammingfullsimplex;

/**
 * Implements the maximum score rule for pivoting for the Simplex algorithm.
 * 
 * Warning: this rule shows high probability of cycling.
 * 
 * Part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org
 * 
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */
public class SimplexMAX extends CanonicalSimplexProgram
{
	/**
	 * Constructor specifying parameters of the linear programming problem.
	 * 
	 * @param A
	 *            constraint matrix
	 * @param b
	 *            constraint values
	 * @param c
	 *            linear objective function weights
	 * @param maximization
	 *            maximization problem if true, minimization else.
	 * @param equality
	 *            indicates whether each constraint is an equality. Else it is
	 *            &lt;=.
	 */
	public SimplexMAX(final double[][] A, final double[] b, final double[] c,
			final boolean maximization, final boolean[] equality) {
		super(A, b, c, maximization, equality);
	}
	/**
	 * Constructor specifying parameters of the linear programming problem.
	 * @param p Parameters of the linear programming problem.
	*/
	public SimplexMAX(final CanonicalProgramParameters p)
	{
		super(p);
	}

	/**
	 * Get the column index for pivoting using the maximum score rule
	 * @param tableau the tableau for Gauss-Jordan elimination
	 * @param max true if maximization problem, false if minization problem
	 * @return the column index in the tableau
	 * */
	@Override
	protected int getPivotColumn(final TableauWithSlackVariables tableau, final boolean max)
	{
		if (max)
		{
			int idx = -1;
			double score = 0;
			for (int i = 0; i <tableau.scoreRow.length; i++)
				if (tableau.scoreRow[i] < score && tableau.originalColOrder[i] >= 0)
				{
					score = tableau.scoreRow[i];
					idx = i;
				}
			if (idx >= 0)
				return idx;
			else
				return -1;
		}
		else
		{
			int idx = -1;
			double score = 0;
			for (int i = 0; i <tableau.scoreRow.length; i++)
				if (tableau.scoreRow[i] > score && tableau.originalColOrder[i] >= 0)
				{
					score = tableau.scoreRow[i];
					idx = i;
				}
			if (idx >= 0)
				return idx;
			else
				return -1;
		}
	}
	/**
	 * Get the row index for pivoting
	 * @param tableau the tableau for Gauss-Jordan elimination
	 * @param colIdx index of the pivoting column in the tableau
	 * @return the row index in the tableau
	 * */
	@Override
	protected int getRowidx(final TableauWithSlackVariables tableau, final int colIdx)
	{
		boolean found = false;
		double maxVal = -1;
		int rowIdx = -1;
		for (int i = 0; i < tableau.xCol.length; i++)
		{
			final double val = tableau.tableau[i][colIdx];
			if (val > tol)
			{
				if (!found)
				{
					found = true;
					maxVal = tableau.xCol[i]/val;
					rowIdx = i;
				}
				else
				{
					if (tableau.xCol[i]/val < maxVal)
					{
						maxVal = tableau.xCol[i]/val;
						rowIdx = i;
					}
				}
			}
		}	
		return rowIdx;
	}
	
	/**
	 * Duplicate the program
	 * 
	 * @param A
	 *            constraint matrix
	 * @param b
	 *            constraint values
	 * @param c
	 *            linear objective function weights
	 * @param maximization
	 *            maximization problem if true, minimization else.
	 * @param equality
	 *            indicates whether each constraint is an equality. Else it is
	 *            &lt;=.
	 */
	@Override
	protected CanonicalSimplexProgram createNewProgram(final double[][] A,
			final double[] b, final double[] c, final boolean maximization, final boolean[] equality) {
		return new SimplexMAX(A, b, c, maximization, equality);
	}
}