/*
 * Copyright (c) 2010-2023. Institut Pasteur.
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

import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.PluginActionable;

public class LinearProgrammingICYPlugin extends PluginActionable {
    @Override
    public void run() {
        new AnnounceFrame("It is now running a series of example optimization problems. See the console for the output.");
        new AnnounceFrame("This is a utility plugin, intended to be used by other ICY plugins.");
        for (int i = 0; i < 3; i++) {
            System.out.println();
            System.out.println("=== Linear Programming test scenario " + i + " ===");
            System.out.println();
            CanonicalSimplexProgram.runExampleScenario(i);
        }
    }

}
