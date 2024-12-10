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

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class TestCanonicalSimplexProgram {
    @Test
    void scenarioZero() {
        final double[][] A = new double[][]{
                {2, 1, 0, 0},
                {0, 1, 4, 1},
                {1, 3, 0, 1}
        };
        final double[] b = new double[]{3, 3, 4};
        final double[] c = new double[]{2, 4, 1, 1};
        final boolean[] equality = new boolean[b.length];
        final boolean maximization = true;

        final CanonicalSimplexProgram program = new SimplexLEXICO(A, b, c, maximization, equality);

        assertTrue(program.solvePrimalSimplex());
        final double score = program.getScore();
        final double[] solution = program.solution;

        assertEquals(6.5d, score);
        assertArrayEquals(new double[]{1.d, 1.d, .5d, .0d}, solution);
    }

    @Test
    void scenarioOne() {
        final double[][] A = new double[][]{
                {1d / 4, -60, -1d / 25, 9},
                {1d / 2, -90, -1d / 50, 3},
                {0, 0, 1, 0}
        };
        final double[] b = new double[]{0, 0, 1};
        final double[] c = new double[]{3d / 4, -150, 1d / 50, -6};
        final boolean[] equality = new boolean[b.length];
        final boolean maximization = true;

        final CanonicalSimplexProgram program = new SimplexLEXICO(A, b, c, maximization, equality);

        assertTrue(program.solvePrimalSimplex());
        final double score = program.getScore();
        final double[] solution = program.solution;

        assertEquals(1.d / 20.d, score);
        assertArrayEquals(new double[]{1.d / 25.d, .0d, 1.d, .0d}, solution);
    }

    @Test
    void scenarioTwo() {
        final double[][] A = new double[][]{
                {5, -4, 13, -2, 1},
                {1, -1, 5, -1, 1},
        };
        final double[] b = new double[]{20, 8};
        final double[] c = new double[]{1, 6, -7, 1, 5};
        final boolean[] equality = new boolean[]{true, true};
        final boolean maximization = false;

        final CanonicalSimplexProgram program = new SimplexLEXICO(A, b, c, maximization, equality);

        assertTrue(program.solvePrimalSimplex());
        final double[] solution = program.solution;

        assertArrayEquals(new double[]{.0d, 4.d / 7.d, 1.d + (5.d / 7.d), .0d, .0d}, solution, .0001d);
    }

    @Test
    void scenarioThree() {
        final double[][] A = new double[][]{
                {1, 2, 1, 1},
                {1, -2, 2, 1},
                {3, -1, 0, -1}
        };
        final double[] b = new double[]{3, -2, -1};
        final double[] c = new double[]{2, -3, 1, 1};
        final boolean[] equality = new boolean[]{true, true, true};
        final boolean maximization = true;

        final CanonicalSimplexProgram program = new SimplexLEXICO(A, b, c, maximization, equality);

        assertTrue(program.solvePrimalSimplex());
        final double score = program.getScore();
        final double[] solution = program.solution;

        assertEquals(-3.d - (1.d / 16.d), score, .0001d);
        assertArrayEquals(new double[]{3.d / 16.d, 1.d + 1.d / 4.d, .0d, 5.d / 16.d}, solution, .0001d);
    }
}
