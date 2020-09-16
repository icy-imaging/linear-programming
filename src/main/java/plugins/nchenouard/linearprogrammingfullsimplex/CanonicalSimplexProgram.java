package plugins.nchenouard.linearprogrammingfullsimplex;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Linear program for the Simplex algorithm. Inheriting classes should implement
 * the pivoting rule.
 * 
 * By definition, the optimization problem is of the form: min_x c'*x or max_x
 * c'*x A*x &lt;= b with x &gt;=0 where x, b, c are column vectors and c' is the
 * transpose of c.
 * 
 * * stands for the matrix multiplication.
 * 
 * Here each constraint <code>A[i]*x</code> can be either "equal" or "lower than
 * or equal" constraints based on the boolean value in the equality constraints
 * vector:
 * <p>
 * A[i]*x == b if equalityConstraints[i]
 * <p>
 * A[i]*x &lt;= b if equalityConstraints[i]
 * <p>
 * 
 * Part of the Linear Programming plugin for ICY http://icy.bioimageanalysis.org
 * 
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */

public abstract class CanonicalSimplexProgram
{	
	/**
	 * Parameters defining the linear programming problem
	 * */
	CanonicalProgramParameters parameters;
	/**
	 * Number of constraints
	 * */
	int numConstraints;
	/**
	 * Number of variables
	 * */
	int numVariables;
	/**
	 * Solution obtained with the Simplex algorithm
	 * */
	double[] solution;
	/**
	 * Intermediary tableau used for some types of problems
	 * */
	TableauWithSlackVariables tableau0 = null;
	/**
	 * Numerical accuracy limit under which a number is considered to be 0
	 * */
	double tol = 1e-12;
	/**
	 * Display detail of steps
	 * */
	public boolean verbose = false;

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
	public CanonicalSimplexProgram(final double[][] A, final double[] b, final double[] c, final boolean maximization, final boolean[] equality)
	{
		this.parameters = new CanonicalProgramParameters();
		this.parameters.A = A;
		this.parameters.b = b;
		this.parameters.c = c;
		this.parameters.maximization = maximization;
		this.parameters.equalityConstraints = equality;
		this.numVariables = c.length;
		this.numConstraints = b.length;
	}
	/**
	 * Constructor specifying parameters of the linear programming problem.
	 * @param p Parameters of the linear programming problem.
	*/
	public CanonicalSimplexProgram(final CanonicalProgramParameters p)
	{
		this.parameters = p;
		this.numVariables = p.c.length;
		this.numConstraints = p.b.length;
	}

	/**
	 * Solve the linear programming problem using the Simplex algorithm for the primal problem
	 * @return true if optimization succeeded, false otherwise.
	 * */
	public boolean solvePrimalSimplex()
	{
		final double[][] A = parameters.A;
		final double[] b = parameters.b;
		final boolean[] eq = parameters.equalityConstraints;
		final double[] c = parameters.c;

		boolean onlyInequality = true;
		for (final boolean e:eq)
		{
			if (e)
			{
				onlyInequality = false;
				break;
			}
		}
		boolean hasNegativeRHS = false;
		for (int i = 0; i  < numConstraints; i++)
		{
			if (b[i] < 0)
			{
				if (eq[i]) // just change the signs for the constraint
				{
					b[i] = -b[i];
					for (int j = 0; j < A[i].length; j++)
						A[i][j] = -A[i][j];
				}
				else
					hasNegativeRHS = true;
			}
		}
		if (onlyInequality)
		{
			if(hasNegativeRHS)
			{
				// add slack variables and multiply by -1 when corresponding to a negative right-hand-side
				final double[][] A2 = new double[numConstraints][numConstraints + numVariables];
				final double[] c2 = new double[numConstraints + numVariables];
				System.arraycopy(c, 0, c2, 0, numVariables); // no cost for slack variables
				final boolean[] equality = new boolean[numConstraints];
				for (int i = 0; i < numConstraints; i++)
				{
					equality[i] = true;
					if (b[i] >= 0)
					{
						System.arraycopy(A[i], 0, A2[i], 0, numVariables);
						A2[i][numVariables + i] = 1;
					}
					else
					{
						for (int j = 0; j < numVariables; j++)
							A2[i][j] = -A[i][j];
						A2[i][numVariables + i] = -1;
						b[i] = -b[i];
					}
				}
				parameters.A = A2;
				parameters.c = c2;
				parameters.equalityConstraints = equality;
				// we have now a problem with only equality constraints and non-negative right-hand-side
				if(solveWithEqualities(false))
				{
					final double[] tmpSolution = this.solution;
					this.solution = new double[numVariables];
					System.arraycopy(tmpSolution, 0, solution, 0, numVariables);
				}
				else
					return false;
			}
			else
				solveWithInequalities(false);
		}
		else
		{
			if(hasNegativeRHS)
			{	
				// we need to solve first a problem to find the initial tableau
				// we will add slack variables and minimize the cost of the sum of the slack variables for the equality constraints
				final boolean maximization = false;
				final boolean[] eq2 = new boolean[numConstraints];
				final boolean[] isNegativeConstraint = new boolean[numConstraints];
				for (int i = 0; i < numConstraints; i++)
				{
					if (b[i] < 0)
					{
						isNegativeConstraint[i] = true;
						b[i] = -b[i];
						for (int j = 0; j < numVariables; j++)
							A[i][j] = -A[i][j];
					}
				}
				final CanonicalSimplexProgram modifiedProgram = createNewProgram(A, b, c, maximization, eq2);
				if (modifiedProgram.solvePhaseICanonicalSimplex())
				{
					modifiedProgram.tableau0.setUnitScoreToSlackVars(eq2);
					if (modifiedProgram.solvePhaseIICanonicalSimplex(modifiedProgram.tableau0, maximization))
					{
						if (Math.abs(modifiedProgram.tableau0.scoreValue) < tol) // optimal score should be 0
						{
							// we have found an initial tableau for the main problem
							// we should modify the scores now for the tableau
							this.tableau0 = modifiedProgram.tableau0;
							tableau0.modifyScores(c);
							return true;
						}
						else
							return false;
					}
					else
						return false;
				}
			}
			else
			{
				return solveWithEqualities(false);
			}
		}
		return true;
	}

	/**
	 * Solve a problem that contains only inequality constraints
	 * @return true if optimization succeeded, false otherwise.
	 * */
	private boolean solveWithInequalities(final boolean compareWithLPSolve)
	{
		// Canonical form: initial tableau is built by introducing unit vectors
		solvePhaseICanonicalSimplex();
		if (verbose)
		{
			System.out.println("Tableau after phase 1 :");
			tableau0.printTableau();
		}
		// Phase II: optimisation by Gauss-Jordan replacement of non-basic columns
		return solvePhaseIICanonicalSimplex(tableau0, parameters.maximization);
	}

	/**
	 * Solve a problem that contains some equality constraints
	 * @return true if optimization succeeded, false otherwise.
	 * */
	private boolean solveWithEqualities(final boolean compareWithLPSolve)
	{
		tableau0 = new TableauMinSlackObjective(this);
		if (verbose)
			tableau0.printTableau();
		// check for some duplicates in tableau0 columns
//		for (int i = 0; i < tableau0.columnVectors.length; i++)
//		{
//			for (int ii = 0; ii < tableau0.columnVectors.length; ii++)
//			{
//				if (i != ii)
//				{
//					int cnt = 0;
//					for (int j = 0; j < tableau0.columnVectors[i].length; j ++)
//					{
//						if (tableau0.columnVectors[i][j] != tableau0.columnVectors[ii][j])
//							cnt++;
//					}
//					if (cnt == 0)
//						System.out.println("Duplicate!!");
//				}
//			}
//		}		
		final boolean success = solvePhaseIICanonicalSimplex(tableau0, false); // minimization of the init tableau
		//		{
		//			double c[] = new double[tableau0.numCol];
		//			for (int i = 0; i < numConstraints; i ++)
		//				if (parameters.equalityConstraints[i])
		//					c[i] = 1;//put weights to slack variables which correspond to equalities
		////			tableau0.solverWithLPSolve(c, b0, false);		
		//		}
		if (!success)
			return false;
		if (verbose)
		{
			System.out.println("Phase 1 ");
			tableau0.printTableau();
		}
		if (Math.abs(tableau0.scoreValue) < tol) // optimal score should be 0
		{
			// check that the initial solution is valid
			//			for (int i = 0; i < A0.length; i++)
			//			{
			//				double v = 0;
			//				double[] s = tableau0.getSolution();
			//				// compute the constraint
			//				for (int j = 0; j < A0[i].length; j++)
			//					v += A0[i][j]*s[j];
			//				if (equality[i])
			//					System.out.println(v+" == "+b0[i]);
			//				else
			//					System.out.println(v+" <= "+b0[i]);
			//			}
			// modifies the score for the objective function
			tableau0.modifyScoresBis(parameters.c);
		}
		else
			return false;
		// phase II: optimization by Gauss-Jordan replacement of non-basic columns
		return solvePhaseIICanonicalSimplex(tableau0, parameters.maximization);
	}

	/**
	 * Solve the phase I of the Simplex algorithm for canonical problem: no equality constraints nor negative constraints
	 * @return true if optimization succeeded, false otherwise.
	 * */
	public boolean solvePhaseICanonicalSimplex()
	{
		// phase I: initial solution with slack variables
		tableau0 = new TableauWithSlackVariables(this.parameters);
		return true;
	}

	/**
	 * Solve the phase II of the Simplex algorithm for canonical problem: no equality constraints nor negative constraints
	 * @return true if optimization succeeded, false otherwise.
	 * */
	public boolean solvePhaseIICanonicalSimplex(final TableauWithSlackVariables tableau, final boolean max)
	{
		if (tableau == null)
			return false;

		// Gauss-Jordan simplex method
		boolean feasible = true;
		boolean reachedOptimum = false;

		while(feasible && !reachedOptimum)
		{
			// start optimization
			final int colIdx = getPivotColumn(tableau, max);
			if (colIdx < 0)// optimal solution achieved
			{
				reachedOptimum = true;
				break;
			}
			else
			{
				final int rowIdx = getRowidx(tableau, colIdx);
				if (rowIdx < 0) // not feasible
				{
					feasible = false;
					break;
				}
				else
				{
					if (verbose)
						System.out.println("replacement: row = "+ rowIdx + ", column  = " + colIdx);
					tableau.pivot(colIdx, rowIdx);
				}
			}
			if (verbose)
			{
				System.out.println("Tableau at iteration:");
				tableau.printTableau();
			}
		}
		if (verbose)
		{
			System.out.println("Final tableau:");
			tableau.printTableau();
		}		if (!feasible)
		{
			if (verbose)
				System.out.println("UNFEASIBLE problem");
			return false;
		}
		else
		{
			solution = tableau.getSolution();
		}
		return true;
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
	protected abstract CanonicalSimplexProgram createNewProgram(double[][] A, double[] b, double[] c, boolean maximization, boolean[] equality);

	/**
	 * Get the row index for pivoting
	 * @param tableau the tableau for Gauss-Jordan elimination
	 * @param colIdx index of the pivoting column in the tableau
	 * @return the row index in the tableau
	 * */
	protected abstract int getRowidx(TableauWithSlackVariables tableau, int colIdx);

	/**
	 * Get the column index for pivoting
	 * @param tableau the tableau for Gauss-Jordan elimination
	 * @param max true if maximization problem, false if minization problem
	 * @return the column index in the tableau
	 * */
	protected abstract int getPivotColumn(TableauWithSlackVariables tableau, boolean max);

	/**
	 * Display whether the computed solution satisfies the constraints
	 * */
	public void checkSolution()
	{
		final double[][] A0 = parameters.A;
		final double[] b0 = parameters.b;
		if (solution != null)
		{
			for (int i = 0; i < numConstraints; i++)
			{
				double v = 0;
				for (int j = 0; j < solution.length; j++)
				{
					v += A0[i][j]*solution[j];
				}
				System.out.println(v + "<?>" + b0[i]);
			}
		}
	}

	/**
	 * Returns the solution of the program, once already computed.
	 * @return computed solution for the linear program
	 * */
	public double[] getSolution()
	{
		return solution;
	}
	/**
	 * Returns the score of the solution of the program, once already computed.
	 * @return score for the computed solution for the linear program
	 * */
	public double getScore()
	{
		if (solution != null)
			return tableau0.scoreValue;
		else
			return 0;
	}
	/**
	 * Returns the number of variable for the optimization problem
	 * @return number of variables
	 * */
	public int getNumVariables() {
		return numVariables;
	}
	/**
	 * Returns the number of constraints for the optimization problem
	 * @return number of constraints
	 * */
	public int getNumConstraints()
	{
		return numConstraints;
	}


	/**
	 * @return Parameters of the linear programming problem
	 * */
	public CanonicalProgramParameters getParameters() {
		return parameters;
	}

	/**
	 * Set the parameters of the linear programming problem
	 * */
	public void setParameters(final CanonicalProgramParameters parameters) {
		this.parameters = parameters;
	}


	public void displayProblem()
	{
		System.out.println("Linear Programming problem max c'x st. Ax <= b, x >= 0 and some equality contraints defined by eq");
		if (parameters.maximization)
			System.out.println("Maximization problem");
		else
			System.out.println("Minimization problem");
		System.out.print("c = [");
		for (int i = 0; i < parameters.c.length; i++)
			System.out.print(parameters.c[i]+", ");
		System.out.println("]");
		System.out.print("b = [");
		for (int i = 0; i < parameters.b.length; i++)
			System.out.print(parameters.b[i]+", ");
		System.out.println("]");
		System.out.print("eq = [");
		for (int i = 0; i < parameters.equalityConstraints.length; i++)
			System.out.print(parameters.equalityConstraints[i]+", ");
		System.out.println("]");
		System.out.print("A = [");
		for (int i = 0; i < parameters.A.length; i++)
		{
			System.out.print("[");
			for (int j = 0; j < parameters.A[i].length; j++)			
			System.out.print(parameters.A[i][j]+", ");
			System.out.println("],");
		}
		System.out.println("]");
	}

	public void displayProblemTutorial()
	{
		System.out.println("Linear Programming problem:");
		System.out.println();
		String line = "";
		if (parameters.maximization)
			line = line + "max_x";
		else
			line = line + "min_x";
		line = line + " ";
		for (int i = 0; i < parameters.c.length; i++)
		{
			if (parameters.c[i] >= 0 && i > 0)
				line = line + " + " + parameters.c[i] + "*x[" + i+"]";
			else
				line = line + " - " + Math.abs(parameters.c[i]) + "*x[" + i+"]";
		}		
		System.out.println(line);
		System.out.println();
		System.out.println("Such that x >= 0 and:");
		System.out.println();
		for (int i = 0; i < parameters.b.length; i++)
		{
			line = "";
			for (int j = 0; j < parameters.A[i].length; j++)
			{
				if (parameters.A[i][j] >= 0 && j > 0)
					line = line + " + " + parameters.A[i][j] + "*x[" + i+"]";
				else
					line = line + " - "+ Math.abs(parameters.A[i][j]) + "*x[" + i+"]";
			}
			if (parameters.equalityConstraints[i])
				line = line + " = ";
			else
				line = line + " <= ";
			line = line + parameters.b[i];
			System.out.println(line);
		}
	}
	
	/**
	 * @return Manual of the library for stand-alone use through command line
	 * */
	public static String getHelp()
	{
		final String helpString = "Run the solver for example problems or for a user-defined system\n" 
				+ "Empty arguments to use the default example.\n"
				+ "Enter an integer from 0 to 3 for different example scenarios.\n"
				+ "Enter -1 for a user-defined scenario through text files.\n "
				+ "For custom scenarios, enter filenames (text files) for different parameters of the problem preceded by the appropriate prefix:\n"
				+ "-c for the objective function file.\n"
				+ "-A for constraint matrix file.\n"
				+ "-b for the constraint value file.\n"
				+ "-e for the equality constraint file.\n"
				+ "-max or -min to indicate a maximization or minimization problem, respectively. Default in minimization.\n"
				+ "Example arguments: java -jar linearProgrammingICY.jar -c c.txt -A A.txt -b b.txt -e eq.txt -max -o solution.txt"
				+ "Each text file must contain a series of double values separated by ',' in a single line, except for the constraint file which contains one line per constraint.\n"
				+ "For the equality file '0' stands for 'false' and '1' for true.\n\n"
				+ "Version 1.0. April 2014. Author: Nicolas Chenouard. nicolas.chenouard.dev@gmail.com. Licence GPL V3.0";
		return helpString;
	}

	/**
	 * Display the manual of the library for stand-alone use through command line
	 * */
	public static void displayHelp()
	{
		System.out.println(getHelp());
	}
	
	
	/**
	 * Run the solver for example problems or for a user-defined system
	 * 
	 * @param args empty to used the default example.
	 * Enter an integer from 0 to 3 for different example scenarios.
	 * Enter -1 for a user-defined scenario through text files
	 * 
	 * For custom scenarios, enter filenames (text files) for different parameters of the problem preceded by the appropriate prefix:
	 * -c for the objective function file
	 * -A for constraint matrix file
	 * -b for the constraint value file
	 * -e for the equality constraint file
	 * -max or -min to indicate a maximization or minimization problem, respectively. Default in minimization.
	 * -v for verbose mode
	 * -o for an output file containing the solution
	 * Example arguments: -c c.txt -A A.txt -b b.txt -e eq.txt -max -o solution.txt
	 * Each text file must contain a series of double values separated by ',' in a single line, except for the constraint file which contains one line per constraint.
	 * For the equality file '0' stands for 'false' and '1' for true.
	 * 
	 * */
	public static void main(final String [] args)
	{
		int scenario = 0;
		if (args.length > 0)
		{
			if(args[0].compareTo("-h") == 0 || args[0].compareTo("-help") == 0 || args[0].compareTo("--h") == 0 || args[0].compareTo("--help") == 0)
			{
				displayHelp();
				return;
			}
			else
			{
				try {
					scenario = Integer.parseInt(args[0]);
				}
				catch (final NumberFormatException e)
				{
					System.out.println("Invalid input argument");
					System.out.println("Use input option -help to display the manual of the software");
					e.printStackTrace();
					return;
				}
			}
		}
		else
		{
			System.out.println("Solving a default simple linear problem with the Simplex algorithm.");
			System.out.println("Enter -help for the manual explicating the customization options.");
		}
		if (scenario < -1 || scenario > 3)
		{
			System.out.println("Invalid input argument");
			System.out.println("Use input option -help to display the manual of the software");
			return;
		}
		boolean maximization = false;
		double[] b = null;
		double[] c = null;
		boolean[] equality = null;
		double[][] A = null;
		String outputFile = null;
		boolean verbose = true;
		switch (scenario) {
		case -1:
		{
			verbose = false;
			for (int i = 0; i < args.length; i++)
			{
				if (args[i].compareTo("-v") == 0)
				{
					verbose = true;
					break;
				}
			}
			for (int i = 0; i < args.length; i++)
			{
				if (args[i].compareTo("-max") == 0)
				{
					maximization = true;
					break;
				}
			}
			String fileA = null;
			for (int i = 0; i < args.length -1; i++)
			{
				if (args[i].compareTo("-A") == 0 || args[i].compareTo("-a") == 0)
				{
					fileA = args[i + 1];
					break;
				}
			}
			if (fileA == null)
			{
				System.out.println("Missing -A input argument");
				System.out.println("Use input option -help to display the manual of the software");
				return;
			}
			String fileB = null;
			for (int i = 0; i < args.length -1; i++)
			{
				if (args[i].compareTo("-b") == 0)
				{
					fileB = args[i + 1];
					break;
				}
			}
			if (fileB == null)
			{
				System.out.println("Missing -b input argument");
				System.out.println("Use input option -help to display the manual of the software");
				return;
			}
			
			String fileC = null;
			for (int i = 0; i < args.length -1; i++)
			{
				if (args[i].compareTo("-c") == 0)
				{
					fileC = args[i + 1];
					break;
				}
			}
			if (fileC == null)
			{
				System.out.println("Missing -c input argument");
				System.out.println("Use input option -help to display the manual of the software");
				return;
			}
			for (int i = 0; i < args.length -1; i++)
			{
				if (args[i].compareTo("-o") == 0)
				{
					outputFile = args[i + 1];
					break;
				}
			}
			String fileE = null;
			for (int i = 0; i < args.length -1; i++)
			{
				if (args[i].compareTo("-e") == 0)
				{
					fileE = args[i + 1];
					break;
				}
			}
			if (fileE == null)
			{
				System.out.println("Missing -e input argument");
				System.out.println("Use input option -help to display the manual of the software");
				return;
			}			
			try {
				File file = new File(fileB);
				FileInputStream fis;
				fis = new FileInputStream(file);
				BufferedReader reader =  new BufferedReader(new InputStreamReader(fis));
				String line = reader.readLine();
				String[] tab = line.split(",");
				b = new double[tab.length-1];
				for (int i = 0; i < tab.length-1; i++)
					b[i] = Double.parseDouble(tab[i]);	
				reader.close();
				fis.close();

				file = new File(fileC);
				fis = new FileInputStream(file);
				reader =  new BufferedReader(new InputStreamReader(fis));
				line = reader.readLine();
				tab = line.split(",");
				c = new double[tab.length-1];
				for (int i = 0; i < tab.length-1; i++)
					c[i] = Double.parseDouble(tab[i]);	
				reader.close();
				fis.close();

				file = new File(fileE);
				fis = new FileInputStream(file);
				reader =  new BufferedReader(new InputStreamReader(fis));
				line = reader.readLine();
				tab = line.split(",");
				equality = new boolean[tab.length-1];
				for (int i = 0; i < tab.length-1; i++)
				{
					final double d = Double.parseDouble(tab[i]);
					equality[i] = d > 0;
				}
				reader.close();
				fis.close();

				final int numRow = b.length;
				final int numCol = c.length;

				file = new File(fileA);
				fis = new FileInputStream(file);
				reader =  new BufferedReader(new InputStreamReader(fis));
				A = new double[numRow][numCol];
				line = reader.readLine();
				int j = 0;
				while (line!=null)
				{
					tab = line.split(",");
					for (int i = 0; i < tab.length-1; i++)
					{
						A[j][i] = Double.parseDouble(tab[i]);
					}
					line = reader.readLine();
					j++;
				}
				reader.close();
				fis.close();
			} catch (final FileNotFoundException e) {
				e.printStackTrace();
				return;
			} catch (final IOException e) {
				e.printStackTrace();
				return;
			}
			break;
		}
		case 0:
		{
			System.out.println("Basic example with 3 pivots");
			System.out.println("Optimal is 6.5");
			System.out.println("Solution is [1.0, 1.0, 0.5, 0.0]");
			A = new double[][]{{2, 1, 0, 0},
					{0, 1, 4, 1},
					{1, 3, 0, 1}
							};
			b = new double[]{3, 3, 4};
			c = new double[]{2, 4, 1, 1};
			equality = new boolean[b.length];
			maximization = true;
			break;
		}
		case 1:
		{
			System.out.println("Example that would cycle with Bland rule.");
			System.out.println("Optimal score is 1/20.");
			System.out.println("Solution is [1/25, 0, 1, 0].");
			A = new double[][]{{1d/4, -60, -1d/25, 9},
									{1d/2, -90, -1d/50, 3},
									{0, 0, 1, 0}
									};
			b = new double[]{0, 0, 1};
			c = new double[]{3d/4, -150, 1d/50, -6};
			equality = new boolean[b.length];
			maximization = true;
			break;
		}
		case 2:
		{
			System.out.println("An example with equality constraints");
			System.out.println("Solution is [0, 4/7, 1 + 5/7, 0, 0].");
			A = new double[][]
							{{5, -4, 13, -2, 1},
							{1, -1, 5, -1, 1},
							};
			b = new double[]{20, 8};
			c = new double[]{1, 6, -7, 1, 5};
			equality = new boolean[]{true, true};
			maximization = false;
			break;
		}
		case 3:
		{
			System.out.println("An example with equality constraints and negative right-hand-side");
			System.out.println("Optimal score is -3 - 1/16.");
			System.out.println("Solution is [3/16, 1 + 1/4, 0, 5/16].");
			A = new double[][]
					{{1, 2, 1, 1},
					{1, -2, 2, 1},
					{3, -1, 0, -1}
					};
			b = new double[]{3, -2, -1};
			c = new double[]{2, -3, 1, 1};
			equality = new boolean[]{true, true, true};
			maximization = true;
		}
		}
		if (A == null || b == null || c == null || equality == null)
		{
			System.out.println("Missing arguments");
			System.out.println("Use input option -help to display the manual of the software");
			return;
		}
		final CanonicalSimplexProgram program = new SimplexLEXICO(A, b, c, maximization, equality);
		if (verbose)
		{
			program.displayProblem();
		}
		if(program.solvePrimalSimplex())
		{
			final double[] solution = program.solution;
			if (verbose)
			{
				System.out.print("Computed solution = [");
				for (int i = 0; i < solution.length; i++)
					System.out.print(solution[i] + ", ");
				System.out.println("]");
				System.out.println("Computed score = " + program.tableau0.scoreValue);
				System.out.println("Computed constraint values:");
				for (int i = 0; i < A.length; i++)
				{
					double v = 0;
					// compute the constraint
					for (int j = 0; j < A[i].length; j++)
						v += A[i][j]*solution[j];
					if (equality[i])
						System.out.println(v+" == "+b[i]);
					else
						System.out.println(v+" <= "+b[i]);
				}
			}
			if (outputFile != null)
			{
	            BufferedWriter out;
				try {
					out = new BufferedWriter(new FileWriter(outputFile));
					out.write("solution\n");
		            for (int i = 0; i < solution.length; i++)
		            {
		            	out.write(solution[i]+", ");
		            }
					out.write("\nscore\n");
					out.write(Double.toString(program.tableau0.scoreValue));
					out.close();					
				} catch (final IOException e) {
					e.printStackTrace();
					return;
				}
			}
			return;
		}
		else
		{
			System.out.println("Solve simplex failed");
			return;
		}
	}
	
	public static void runExampleScenario(final int scenario)
	{
		if (scenario < 0 || scenario > 3)
		{
			System.out.println("Invalid input argument");
			System.out.println("Scenario indices are: 0, 1, 2 and 3");
			return;
		}
		boolean maximization = false;
		double[] b = null;
		double[] c = null;
		boolean[] equality = null;
		double[][] A = null;
		final boolean verbose = true;
		switch (scenario) {
		case 0:
		{
			System.out.println("Basic example with 3 pivots");
			System.out.println("Optimal score is 6.5");
			System.out.println("Solution is [1.0, 1.0, 0.5, 0.0]");
			A = new double[][]{{2, 1, 0, 0},
					{0, 1, 4, 1},
					{1, 3, 0, 1}
							};
			b = new double[]{3, 3, 4};
			c = new double[]{2, 4, 1, 1};
			equality = new boolean[b.length];
			maximization = true;
			break;
		}
		case 1:
		{
			System.out.println("Example that would cycle with Bland rule.");
			System.out.println("Optimal score is 1/20.");
			System.out.println("Solution is [1/25, 0, 1, 0].");
			A = new double[][]{{1d/4, -60, -1d/25, 9},
									{1d/2, -90, -1d/50, 3},
									{0, 0, 1, 0}
									};
			b = new double[]{0, 0, 1};
			c = new double[]{3d/4, -150, 1d/50, -6};
			equality = new boolean[b.length];
			maximization = true;
			break;
		}
		case 2:
		{
			System.out.println("An example with equality constraints");
			System.out.println("Solution is [0, 4/7, 1 + 5/7, 0, 0].");
			A = new double[][]
							{{5, -4, 13, -2, 1},
							{1, -1, 5, -1, 1},
							};
			b = new double[]{20, 8};
			c = new double[]{1, 6, -7, 1, 5};
			equality = new boolean[]{true, true};
			maximization = false;
			break;
		}
		case 3:
		{
			System.out.println("An example with equality constraints and negative right-hand-side");
			System.out.println("Optimal score is -3 - 1/16.");
			System.out.println("Solution is [3/16, 1 + 1/4, 0, 5/16].");
			A = new double[][]
					{{1, 2, 1, 1},
					{1, -2, 2, 1},
					{3, -1, 0, -1}
					};
			b = new double[]{3, -2, -1};
			c = new double[]{2, -3, 1, 1};
			equality = new boolean[]{true, true, true};
			maximization = true;
			break;
		}
		}
		final CanonicalSimplexProgram program = new SimplexLEXICO(A, b, c, maximization, equality);
		if (verbose)
		{
			System.out.println();
			program.displayProblemTutorial();
		}
		if(program.solvePrimalSimplex())
		{
			final double[] solution = program.solution;
			if (verbose)
			{
				System.out.println();
				System.out.print("Computed solution = [");
				for (int i = 0; i < solution.length; i++)
					System.out.print(solution[i] + ", ");
				System.out.println("]");
				System.out.println();
				System.out.println("Computed score = " + program.tableau0.scoreValue);
				System.out.println();
				System.out.println("Computed constraint values:");
				for (int i = 0; i < A.length; i++)
				{
					double v = 0;
					// compute the constraint
					for (int j = 0; j < A[i].length; j++)
						v += A[i][j]*solution[j];
					if (equality[i])
						System.out.println(v+" == "+b[i]);
					else
						System.out.println(v+" <= "+b[i]);
				}
			}
			return;
		}
		else
		{
			System.out.println("Solve simplex failed");
			return;
		}
	}
}