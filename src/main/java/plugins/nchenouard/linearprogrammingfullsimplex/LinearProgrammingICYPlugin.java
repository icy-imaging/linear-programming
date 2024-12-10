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

import org.bioimageanalysis.icy.extension.plugin.abstract_.PluginActionable;
import org.bioimageanalysis.icy.extension.plugin.annotation_.IcyPluginIcon;
import org.bioimageanalysis.icy.extension.plugin.annotation_.IcyPluginName;
import org.bioimageanalysis.icy.gui.frame.progress.AnnounceFrame;
import org.bioimageanalysis.icy.system.logging.IcyLogger;

@IcyPluginName("Linear Programming Example")
@IcyPluginIcon(path = "/plugins/nchenouard/linearprogramming/linear-programming.png")
public class LinearProgrammingICYPlugin extends PluginActionable {
    @Override
    public void run() {
        new AnnounceFrame("It is now running a series of example optimization problems. See the console for the output.");
        new AnnounceFrame("This is an utility plugin, intended to be used by other ICY plugins.");
        for (int i = 0; i < 3; i++) {
            IcyLogger.info(this.getClass(), "=== Linear Programming test scenario " + i + " ===");
            CanonicalSimplexProgram.runExampleScenario(i);
        }
    }
}
