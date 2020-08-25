package plugins.nchenouard.linearprogrammingfullsimplex;

import java.util.ArrayList;

/**
 * Implements the lexicographic rule for pivoting for the Simplex algorithm.
 * 
 * Part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org
 * 
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 * @version 1.0
 * @date 2014-04-21
 * @license gpl v3.0
 * 
 */
public class SimplexLEXICO extends CanonicalSimplexProgram
{
	/**
	 * Constructor specifying parameters of the linear programming problem.
	 * @param A constraint matrix
	 * @param b constraint values
	 * @param c linear objective function weights
	 * @param maximization maximization problem if true, minimization else.
	 * @param equality indicates whether each constraint is an equality. Else it is <=.
	 * */
	public SimplexLEXICO(double[][] A, double[] b, double[] c,
			boolean maximization, boolean[] equality) {
		super(A, b, c, maximization, equality);
	}
	/**
	 * Constructor specifying parameters of the linear programming problem.
	 * @param p Parameters of the linear programming problem.
	 */
	public SimplexLEXICO(CanonicalProgramParameters p) {
		super(p);
	}

	/**
	 * Get the column index for pivoting using the lexicographic order rule
	 * @param tableau the tableau for Gauss-Jordan elimination
	 * @param max true if maximization problem, false if minization problem
	 * @return the column index in the tableau
	 * */
	@Override
	protected int getPivotColumn(TableauWithSlackVariables tableau, boolean max) {
		// max score choice strategy
		if (max)
		{
			int idx = -1;
			double score = 0;
			for (int i = 0; i <tableau.scoreRow.length; i++)
				if (tableau.scoreRow[i] < score)
				{
					if (tableau.originalColOrder[i] >= 0 || !parameters.equalityConstraints[i])
					{
						score = tableau.scoreRow[i];
						idx = i;
					}
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
				if (tableau.scoreRow[i] > score)
					if (tableau.originalColOrder[i] >= 0 || !parameters.equalityConstraints[i])
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
	protected int getRowidx(TableauWithSlackVariables tableau, int colIdx)
	{
		int rowIdx = -1;
		int numBasicRow = 0;
		ArrayList<Integer> basicCols = new ArrayList<Integer>();
		for (int i = 0; i < tableau.originalColOrder.length; i++)
			if (tableau.originalColOrder[i] < 0)
			{
				basicCols.add(i);
				numBasicRow ++;
			}
		double[] minLexico = new double[numBasicRow + 1];
		boolean found = false;
		for (int i = 0; i < tableau.xCol.length; i++)
		{
			double val = tableau.tableau[i][colIdx];
			if (val > tol)
			{
				int k = 0;
				boolean isMinLexico = false;
				if (!found)
				{
					found = true;
					isMinLexico = true;
				}
				else
				{
					if (tableau.xCol[i]/val == minLexico[k])
					{
						k++;
						for (Integer ii:basicCols)
						{
							double t = tableau.tableau[i][ii]/val;
							if (t > minLexico[k])
							{
								isMinLexico = false;
								break;
							}
							else if (t < minLexico[k])
							{
								isMinLexico = true;
								break;
							}
							k++;
						}
					}
					else if (tableau.xCol[i]/val < minLexico[k])
						isMinLexico = true;
				}
				if (isMinLexico)
				{
					if (k == 0)
					{
						minLexico[k] = tableau.xCol[i]/val;
						k++;
					}
					while (k <= basicCols.size() )
					{
						minLexico[k] = tableau.tableau[i][basicCols.get(k - 1)]/val;
						k++;
					}						
					rowIdx = i;
				}
			}
		}	
		return rowIdx;
	}
	/**
	 * Duplicate the program
	 * @param A constraint matrix
	 * @param b constraint values
	 * @param c linear objective function weights
	 * @param maximization maximization problem if true, minimization else.
	 * @param equality indicates whether each constraint is an equality. Else it is <=.
	 * */
	@Override
	protected CanonicalSimplexProgram createNewProgram(double[][] A,
			double[] b, double[] c, boolean maximization, boolean[] equality) {
		return new SimplexLEXICO(A, b, c, maximization, equality);
	}
}