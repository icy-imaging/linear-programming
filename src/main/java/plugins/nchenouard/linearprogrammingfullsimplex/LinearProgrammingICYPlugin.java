package plugins.nchenouard.linearprogrammingfullsimplex;


import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.PluginActionable;

public class LinearProgrammingICYPlugin extends PluginActionable
{
	@Override
	public void run()
	{
		new AnnounceFrame("It is now running a series of example optimization problems. See the console for the output.");
		new AnnounceFrame("This is a utility plugin, intended to be used by other ICY plugins.");
		for (int i = 0; i < 3; i++)
		{
			System.out.println();
			System.out.println("=== Linear Programming test scenario "+i+" ===");
			System.out.println();
			CanonicalSimplexProgram.runExampleScenario(i);
		}
	}
	
}
