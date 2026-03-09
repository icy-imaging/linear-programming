/*
 * Copyright (c) 2010-2026. Institut Pasteur.
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

import java.util.ArrayList;

/**
 * Implements the lexicographic rule for pivoting for the Simplex algorithm.
 * <p>
 * Part of the Linear Programming plugin for ICY <a href="https://icy.bioimageanalysis.org">https://icy.bioimageanalysis.org</a>
 *
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */
public class SimplexLEXICO extends CanonicalSimplexProgram {
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
    public SimplexLEXICO(final double[][] A, final double[] b, final double[] c, final boolean maximization, final boolean[] equality) {
        super(A, b, c, maximization, equality);
    }

    /**
     * Constructor specifying parameters of the linear programming problem.
     *
     * @param p Parameters of the linear programming problem.
     */
    public SimplexLEXICO(final CanonicalProgramParameters p) {
        super(p);
    }

    /**
     * Get the column index for pivoting using the lexicographic order rule
     *
     * @param tableau the tableau for Gauss-Jordan elimination
     * @param max     true if maximization problem, false if minization problem
     * @return the column index in the tableau
     */
    @Override
    protected int getPivotColumn(final TableauWithSlackVariables tableau, final boolean max) {
        // max score choice strategy
        int idx = -1;
        double score = 0;

        if (max) {
            for (int i = 0; i < tableau.scoreRow.length; i++)
                if (tableau.scoreRow[i] < score) {
                    if (tableau.originalColOrder[i] >= 0 || !parameters.equalityConstraints[i]) {
                        score = tableau.scoreRow[i];
                        idx = i;
                    }
                }
        }
        else {
            for (int i = 0; i < tableau.scoreRow.length; i++)
                if (tableau.scoreRow[i] > score)
                    if (tableau.originalColOrder[i] >= 0 || !parameters.equalityConstraints[i]) {
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
        int rowIdx = -1;
        int numBasicRow = 0;
        final ArrayList<Integer> basicCols = new ArrayList<>();
        for (int i = 0; i < tableau.originalColOrder.length; i++)
            if (tableau.originalColOrder[i] < 0) {
                basicCols.add(i);
                numBasicRow++;
            }
        final double[] minLexico = new double[numBasicRow + 1];
        boolean found = false;
        for (int i = 0; i < tableau.xCol.length; i++) {
            final double val = tableau.tableau[i][colIdx];
            if (val > tol) {
                int k = 0;
                boolean isMinLexico = false;
                if (!found) {
                    found = true;
                    isMinLexico = true;
                }
                else {
                    if (tableau.xCol[i] / val == minLexico[k]) {
                        k++;
                        for (final Integer ii : basicCols) {
                            final double t = tableau.tableau[i][ii] / val;
                            if (t > minLexico[k]) {
                                isMinLexico = false;
                                break;
                            }
                            else if (t < minLexico[k]) {
                                isMinLexico = true;
                                break;
                            }
                            k++;
                        }
                    }
                    else if (tableau.xCol[i] / val < minLexico[k])
                        isMinLexico = true;
                }
                if (isMinLexico) {
                    if (k == 0) {
                        minLexico[k] = tableau.xCol[i] / val;
                        k++;
                    }
                    while (k <= basicCols.size()) {
                        minLexico[k] = tableau.tableau[i][basicCols.get(k - 1)] / val;
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
        return new SimplexLEXICO(A, b, c, maximization, equality);
    }
}