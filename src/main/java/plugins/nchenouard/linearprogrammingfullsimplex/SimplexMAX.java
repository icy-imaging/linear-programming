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

import org.jetbrains.annotations.NotNull;

/**
 * Implements the maximum score rule for pivoting for the Simplex algorithm.
 * <p>
 * Warning: this rule shows high probability of cycling.
 * <p>
 * Part of the Linear Programming plugin for ICY <a href="https://icy.bioimageanalysis.org">https://icy.bioimageanalysis.org</a>
 *
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */
public class SimplexMAX extends CanonicalSimplexProgram {
    /**
     * Constructor specifying parameters of the linear programming problem.
     *
     * @param A            constraint matrix
     * @param b            constraint values
     * @param c            linear objective function weights
     * @param maximization maximization problem if true, minimization else.
     * @param equality     indicates whether each constraint is an equality. Else it is
     *                     &lt;=.
     */
    public SimplexMAX(final double[][] A, final double[] b, final double[] c, final boolean maximization, final boolean[] equality) {
        super(A, b, c, maximization, equality);
    }

    /**
     * Constructor specifying parameters of the linear programming problem.
     *
     * @param p Parameters of the linear programming problem.
     */
    public SimplexMAX(final CanonicalProgramParameters p) {
        super(p);
    }

    /**
     * Get the column index for pivoting using the maximum score rule
     *
     * @param tableau the tableau for Gauss-Jordan elimination
     * @param max     true if maximization problem, false if minization problem
     * @return the column index in the tableau
     */
    @Override
    protected int getPivotColumn(final TableauWithSlackVariables tableau, final boolean max) {
        int idx = -1;
        double score = 0;

        if (max) {
            for (int i = 0; i < tableau.scoreRow.length; i++)
                if (tableau.scoreRow[i] < score && tableau.originalColOrder[i] >= 0) {
                    score = tableau.scoreRow[i];
                    idx = i;
                }
        }
        else {
            for (int i = 0; i < tableau.scoreRow.length; i++)
                if (tableau.scoreRow[i] > score && tableau.originalColOrder[i] >= 0) {
                    score = tableau.scoreRow[i];
                    idx = i;
                }
        }

        if (idx >= 0)
            return idx;
        else
            return -1;
    }

    /**
     * Get the row index for pivoting
     *
     * @param tableau the tableau for Gauss-Jordan elimination
     * @param colIdx  index of the pivoting column in the tableau
     * @return the row index in the tableau
     */
    @Override
    protected int getRowidx(final @NotNull TableauWithSlackVariables tableau, final int colIdx) {
        boolean found = false;
        double maxVal = -1;
        int rowIdx = -1;
        for (int i = 0; i < tableau.xCol.length; i++) {
            final double val = tableau.tableau[i][colIdx];
            if (val > tol) {
                if (!found) {
                    found = true;
                    maxVal = tableau.xCol[i] / val;
                    rowIdx = i;
                }
                else {
                    if (tableau.xCol[i] / val < maxVal) {
                        maxVal = tableau.xCol[i] / val;
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
     * @param A            constraint matrix
     * @param b            constraint values
     * @param c            linear objective function weights
     * @param maximization maximization problem if true, minimization else.
     * @param equality     indicates whether each constraint is an equality. Else it is
     *                     &lt;=.
     */
    @Override
    protected CanonicalSimplexProgram createNewProgram(final double[][] A, final double[] b, final double[] c, final boolean maximization, final boolean[] equality) {
        return new SimplexMAX(A, b, c, maximization, equality);
    }
}