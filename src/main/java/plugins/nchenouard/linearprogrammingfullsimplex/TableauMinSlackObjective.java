package plugins.nchenouard.linearprogrammingfullsimplex;

/**
 * Tableau used for Gauss-Jordan elimination in which slack variables are given a 1 score depending on whether constraints are equalities.
 * Slack variables are added to the set of column vectors in a leading position.
 * 
 * Part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org
 * 
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 * @version 1.0
 * @date 2014-04-21
 * @license gpl v3.0
 * 
 */
public class TableauMinSlackObjective extends TableauWithSlackVariables
{

	/**
	 * Create the tableau based on a linear program.
	 * The scores in the tableau will not be initialized to their standard value but set to handle equality constraints.
	 * 
	 * @param parameters parameters defining the linear optimization problem.
	 */
	public TableauMinSlackObjective(CanonicalSimplexProgram simplexProgram) {
		super(simplexProgram.parameters.A);
		// now make the slacks feasible
		xCol = simplexProgram.parameters.b.clone(); // equal to b since we start with unit vectors as initial left vectors
		// build the score row
		scoreRow = new double[numCol];
		// add a score for slacks corresponding to equality constraints
		double[] slackScores = new double[simplexProgram.parameters.equalityConstraints.length];
		this.scoreValue = 0;
		for (int i = 0; i < slackScores.length; i++)
			if (simplexProgram.parameters.equalityConstraints[i])
			{
				scoreRow[i] = 1;
				scoreValue += xCol[i];
			}
		for (int i = 0; i < numCol; i++)
		{
			if (originalColOrder[i] >= 0)
			{
				double pi = 0;
				for (int j = 0; j < numRows; j++)
				{
					if (simplexProgram.parameters.equalityConstraints[j])
						pi += tableau[j][i];
				}
				scoreRow[i] = pi;
			}
		}
	}
	
	/**
	 * Initialize scores for the tableau.
	 * The scores are not set to their standard value.
	 * Instead score 1 is given to slack variables corresponding to equality constraints and 0 to others.
	 * 
	 * @param c standard score for variables.
	 * @param slackScore score values for slack variables.
	 */
	public void initScores(double[] c, double[] slackScore)
	{
		double[] scoreValues = new double[numCol];
		for (int i = 0; i < slackScore.length; i++)
			scoreValues[i] = slackScore[i];
		for (int i = 0; i < numCol; i++)
		{
			double pi = 0;
			for (int j = 0; j < numRows; j++)
			{
				if (leftVectorsIdx[j] >= 0)
					pi += tableau[j][i]*scoreValues[leftVectorsIdx[j]];
				// else do nothing as unit vectors have no cost
			}
			if (originalColOrder[i] >= 0)
				scoreRow[i] = pi - c[originalColOrder[i]];
			else				
				scoreRow[i] = pi;
		}
		scoreValue = 0;
		for (int i = 0; i < numRows; i++)
		{
			if (leftVectorsIdx[i] >= 0 && originalColOrder[leftVectorsIdx[i]] >= 0)
			{
				scoreValue += c[originalColOrder[leftVectorsIdx[i]]]*xCol[i];
			}
		}
	}
}
