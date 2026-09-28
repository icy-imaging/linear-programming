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

import fr.icy.shared.logging.CustomLevel;
import org.jetbrains.annotations.NotNull;

import java.util.logging.Logger;

/**
 * Tableau used for Gauss-Jordan elimination.
 * Initial basis is built by using slack variables.
 * Slack variables are added to the set of column vectors in a leading position.
 * <p>
 * Part of the Linear Programming plugin for ICY <a href="https://icy.bioimageanalysis.org">https://icy.bioimageanalysis.org</a>
 *
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */
public class TableauWithSlackVariables {
    private static final Logger LOGGER = Logger.getLogger(TableauWithSlackVariables.class.getName());

    /**
     * Set of column vectors for the tableau. The first numVariables ones are unit vectors corresponding to slack variables.
     */
    double[][] columnVectors;
    /**
     * Basis vectors for the tableau. In the initial state they are unit vectors corresponding to slack variables.
     */
    double[][] basisVectors;
    /**
     * Indices of basis vectors in the set of column vectors columnVectors.
     */
    int[] leftVectorsIdx;
    /**
     * Index of column vectors in the original constraint matrix. Unit vectors added for slack variables have index -1.
     */
    int[] originalColOrder;
    /**
     * Row of the tableau corresponding to scores
     */
    double[] scoreRow;
    /**
     * Number of columns of the tableau
     */
    int numCol;
    /**
     * Number of rows of the tableau
     */
    int numRows;
    /**
     * Tableau for Gauss-Jordan elimination.
     */
    double[][] tableau;
    /**
     * Solution column. First rows correspond to slack variables.
     */
    double[] xCol;
    /**
     * Number of variables
     */
    int numVariables;
    /**
     * Score value for the current solution
     */
    double scoreValue;

    /**
     * Create the tableau based on a constraint matrix.
     * Each element A[i] defines a linear combination of the variable corresponding to a constraint.
     *
     * @param A table of linear combinations of the variable, each defining a constraint.
     */
    public TableauWithSlackVariables(final double[][] A) {
        createVectors(A);
    }

    /**
     * Create the tableau based on a linear program.
     *
     * @param parameters parameters defining the linear optimization problem.
     */
    public TableauWithSlackVariables(final @NotNull CanonicalProgramParameters parameters) {
        createVectors(parameters.A); // constraint value column
        xCol = parameters.b.clone(); // equal to b since we start with unit vectors as initial left vectors
        // build the score row
        scoreRow = new double[numCol];
        modifyScores(parameters.c);
    }

    /**
     * Change the score of slack variables to 1
     *
     * @param changeSlack indicate whether the score of each slack variable should be changed
     */
    protected void setUnitScoreToSlackVars(final boolean[] changeSlack) {
        final double[] scoreValues = new double[numCol];
        // recompute scoreValues based on c
        int k = 0;
        for (int i = 0; i < originalColOrder.length; i++) {
            if (originalColOrder[i] < 0) {
                if (changeSlack[k])
                    scoreValues[i] = 1d;
                k++;
            }
        }
        scoreRow = new double[numCol];
        for (int i = 0; i < numCol; i++) {
            double pi = 0;
            for (int j = 0; j < numRows; j++) {
                if (leftVectorsIdx[j] >= 0)
                    pi += tableau[j][i] * scoreValues[leftVectorsIdx[j]];
            }
            scoreRow[i] = pi - scoreValues[i];
        }
        scoreValue = 0;
        for (int i = 0; i < leftVectorsIdx.length; i++) {
            scoreValue += scoreValues[leftVectorsIdx[i]] * xCol[i];
        }

        k = 0;
        for (int i = 0; i < originalColOrder.length; i++) {
            if (originalColOrder[i] < 0) {
                if (changeSlack[k])
                    scoreRow[i] = 1;
                k++;
            }
        }
    }

    /**
     * Create column and basis vectors based on the constraint matrix
     *
     * @param A Constraint matrix. Each element A[i] defines a linear combination of the variable corresponding to a constraint.
     */
    protected void createVectors(final double @NotNull [] @NotNull [] A) {
        final int numConstraints = A.length;
        numVariables = A[0].length;
        numCol = numConstraints + numVariables;
        numRows = numConstraints;
        tableau = new double[numRows][numCol];

        columnVectors = new double[numConstraints + numVariables][numConstraints];
        originalColOrder = new int[numConstraints + numVariables];
        basisVectors = new double[numConstraints][];
        leftVectorsIdx = new int[numConstraints];
        // set unit vectors as left vectors of the tableau, and add them to the column vector list
        for (int i = 0; i < numConstraints; i++) {
            columnVectors[i][i] = 1d;
            tableau[i][i] = 1d;
            basisVectors[i] = columnVectors[i];
            originalColOrder[i] = -1;
            leftVectorsIdx[i] = i;
        }
        // build the set of columns corresponding to initial variables
        for (int i = 0; i < numVariables; i++) {
            for (int j = 0; j < numConstraints; j++)
                columnVectors[i + numConstraints][j] = A[j][i];
            originalColOrder[i + numConstraints] = i;
        }
        for (int j = 0; j < numConstraints; j++)
            System.arraycopy(A[j], 0, tableau[j], numConstraints, numVariables);
    }

    /**
     * Print the tableau
     */
    public void printTableau() {
        if (!LOGGER.isLoggable(CustomLevel.TRACE))
            return;
        LOGGER.log(CustomLevel.TRACE, "Column vectors");
        for (int i = 0; i < columnVectors.length; i++) {
            final double[] vector = columnVectors[i];
            LOGGER.log(CustomLevel.TRACE, "a" + i + " = [");
            for (final double v : vector)
                LOGGER.log(CustomLevel.TRACE, v + ", ");
            LOGGER.log(CustomLevel.TRACE, "]'");
        }
        LOGGER.log(CustomLevel.TRACE, "Basis vectors");
        for (int i = 0; i < leftVectorsIdx.length; i++) {
            final double[] vector = columnVectors[leftVectorsIdx[i]];
            LOGGER.log(CustomLevel.TRACE, "v" + i + " = [");
            for (final double v : vector)
                LOGGER.log(CustomLevel.TRACE, v + ", ");
            LOGGER.log(CustomLevel.TRACE, "]'");
        }
        LOGGER.log(CustomLevel.TRACE, "Tableau");
        for (final double[] doubles : tableau) {
            for (final double aDouble : doubles)
                LOGGER.log(CustomLevel.TRACE, aDouble + " ");
        }

        LOGGER.log(CustomLevel.TRACE, "x = [");
        for (final double v : xCol)
            LOGGER.log(CustomLevel.TRACE, v + ", ");
        LOGGER.log(CustomLevel.TRACE, "]'");

        LOGGER.log(CustomLevel.TRACE, "score row = [");
        for (final double v : scoreRow)
            LOGGER.log(CustomLevel.TRACE, v + ", ");
        LOGGER.log(CustomLevel.TRACE, "]'");

        LOGGER.log(CustomLevel.TRACE, "Score = " + scoreValue);
    }

    /**
     * Pivot at a given position in the tableau
     *
     * @param colIdx index of the pivoting column in the tableau
     * @param rowIdx index of the pivoting row in the tableau
     */
    public void pivot(final int colIdx, final int rowIdx) {
        leftVectorsIdx[rowIdx] = colIdx;

        final double pivotVal = tableau[rowIdx][colIdx];
        final double[] pivotRow = tableau[rowIdx];
        final double[] pivotCol = new double[numRows];
        for (int i = 0; i < numRows; i++)
            pivotCol[i] = tableau[i][colIdx];
        final double scorePivot = scoreRow[colIdx];

        // update the score row
        for (int j = 0; j < colIdx; j++)
            scoreRow[j] = scoreRow[j] - pivotRow[j] * scorePivot / pivotVal;
        for (int j = colIdx + 1; j < numCol; j++)
            scoreRow[j] = scoreRow[j] - pivotRow[j] * scorePivot / pivotVal;
        scoreRow[colIdx] = 0;

        scoreValue = scoreValue - xCol[rowIdx] * scorePivot / pivotVal;

        // update the variable column
        for (int i = 0; i < rowIdx; i++)
            xCol[i] = xCol[i] - xCol[rowIdx] * pivotCol[i] / pivotVal;
        for (int i = rowIdx + 1; i < numRows; i++)
            xCol[i] = xCol[i] - xCol[rowIdx] * pivotCol[i] / pivotVal;
        xCol[rowIdx] = xCol[rowIdx] / pivotVal;

        // update the tableau
        for (int i = 0; i < rowIdx; i++) {
            for (int j = 0; j < colIdx; j++)
                tableau[i][j] = tableau[i][j] - pivotCol[i] * pivotRow[j] / pivotVal; // tij - til*tkj/tkl for col l and row k
            for (int j = colIdx + 1; j < numCol; j++)
                tableau[i][j] = tableau[i][j] - pivotCol[i] * pivotRow[j] / pivotVal;
        }
        for (int i = rowIdx + 1; i < numRows; i++) {
            for (int j = 0; j < colIdx; j++)
                tableau[i][j] = tableau[i][j] - pivotCol[i] * pivotRow[j] / pivotVal;
            for (int j = colIdx + 1; j < numCol; j++)
                tableau[i][j] = tableau[i][j] - pivotCol[i] * pivotRow[j] / pivotVal;
        }
        for (int j = 0; j < numCol; j++)
            tableau[rowIdx][j] = pivotRow[j] / pivotVal;
        for (int i = 0; i < numRows; i++)
            tableau[i][colIdx] = 0;
        tableau[rowIdx][colIdx] = 1;
    }

    /**
     * Get the solution. Slack variables are removed.
     *
     * @return a table of variable values defining the current solution
     */
    public double[] getSolution() {
        final double[] solution = new double[numVariables];
        for (int i = 0; i < leftVectorsIdx.length; i++) {
            final int originalIdx = originalColOrder[leftVectorsIdx[i]];
            if (originalIdx >= 0)
                solution[originalIdx] = xCol[i];

        }
        return solution;
    }

    /**
     * Modify scores in the tableau according to a set of scores for the constraints
     */
    public void modifyScores(final double[] c) {
        final double[] scoreValues = new double[numCol];
        for (int i = 0; i < numCol; i++) {
            if (originalColOrder[i] >= 0)
                scoreValues[i] = c[originalColOrder[i]];
            else
                scoreValues[i] = 0;
        }
        for (int i = 0; i < numCol; i++) {
            double pi = 0;
            for (int j = 0; j < numRows; j++) {
                if (leftVectorsIdx[j] >= 0)
                    pi += tableau[j][i] * scoreValues[leftVectorsIdx[j]];
                // else do nothing as unit vectors have no cost
            }
            if (originalColOrder[i] >= 0)
                scoreRow[i] = pi - c[originalColOrder[i]];
            else
                scoreRow[i] = pi;
        }
        scoreValue = 0;
        for (int i = 0; i < numRows; i++) {
            if (leftVectorsIdx[i] >= 0 && originalColOrder[leftVectorsIdx[i]] >= 0) {
                scoreValue += c[originalColOrder[leftVectorsIdx[i]]] * xCol[i];
            }
        }
    }

    /**
     * Modify scores in the tableau according to a set of scores for the constraints
     * Same as modifyScoresBis. Just a backup function. Still working on it.
     */
    public void modifyScoresBis(final double[] c) {
        final double[] scoreValues = new double[numCol];
        for (int i = 0; i < numCol; i++) {
            if (originalColOrder[i] >= 0)
                scoreValues[i] = c[originalColOrder[i]];
            else
                scoreValues[i] = 0;
        }
        for (int i = 0; i < numCol; i++) {
            double pi = 0;
            for (int j = 0; j < numRows; j++)
                pi += tableau[j][i] * scoreValues[leftVectorsIdx[j]];
            scoreRow[i] = pi - scoreValues[i];
        }
        scoreValue = 0;
        for (int i = 0; i < numRows; i++)
            scoreValue += scoreValues[leftVectorsIdx[i]] * xCol[i];
    }
}