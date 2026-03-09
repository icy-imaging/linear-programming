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

/**
 * Tableau used for Gauss-Jordan elimination in which slack variables are given a 1 score depending on whether constraints are equalities.
 * Slack variables are added to the set of column vectors in a leading position.
 * <p>
 * Part of the Linear Programming plugin for ICY <a href="https://icy.bioimageanalysis.org">https://icy.bioimageanalysis.org</a>
 *
 * @author Nicolas Chenouard (nicolas.chenouard.dev@gmail.com)
 */
public class TableauMinSlackObjective extends TableauWithSlackVariables {
    /**
     * Create the tableau based on a linear program. The scores in the tableau
     * will not be initialized to their standard value but set to handle
     * equality constraints.
     *
     * @param simplexProgram parameters defining the linear optimization problem.
     */
    public TableauMinSlackObjective(final @NotNull CanonicalSimplexProgram simplexProgram) {
        super(simplexProgram.parameters.A);
        // now make the slacks feasible
        xCol = simplexProgram.parameters.b.clone(); // equal to b since we start with unit vectors as initial left vectors
        // build the score row
        scoreRow = new double[numCol];
        // add a score for slacks corresponding to equality constraints
        final double[] slackScores = new double[simplexProgram.parameters.equalityConstraints.length];
        this.scoreValue = 0;
        for (int i = 0; i < slackScores.length; i++)
            if (simplexProgram.parameters.equalityConstraints[i]) {
                scoreRow[i] = 1;
                scoreValue += xCol[i];
            }
        for (int i = 0; i < numCol; i++) {
            if (originalColOrder[i] >= 0) {
                double pi = 0;
                for (int j = 0; j < numRows; j++) {
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
     * @param c          standard score for variables.
     * @param slackScore score values for slack variables.
     */
    public void initScores(final double[] c, final double[] slackScore) {
        final double[] scoreValues = new double[numCol];
        System.arraycopy(slackScore, 0, scoreValues, 0, slackScore.length);
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
}
