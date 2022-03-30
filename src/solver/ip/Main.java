package solver.ip;

import java.nio.file.Path;
import java.nio.file.Paths;

import ilog.concert.IloException;

public class Main
{  
  public static void main(String[] args) throws Exception
  {
		if(args.length == 0)
		{
			System.out.println("Usage: java Main <file>");
			return;
		}
		
		String input = args[0];
		Path path = Paths.get(input);
		String filename = path.getFileName().toString();

		Timer watch = new Timer();
		watch.start();
		
		IPInstance instance = DataParser.parseIPFile(input);
		System.out.println(instance);
		double result = instance.solve();
		String result_str;
		if (result==-1) {
			result_str= "\"--\"";
		} else {
			result_str= String.format("\"%.2f\"",result);
		}
		watch.stop();
		System.out.println("{\"Instance\": \"" + filename +
				"\", \"Time\": " + String.format("%.2f",watch.getTime()) +
				", \"Result\": " + result_str + 
				", \"Solution\": \"OPT?\"}");
  }
}
