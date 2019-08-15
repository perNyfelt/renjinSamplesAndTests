package se.alipsa;

import org.junit.Test;
import org.renjin.eval.EvalException;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.script.ScriptException;

import static org.junit.Assert.fail;

public class SimpleRenjinWithDBITest {

   Logger log = LoggerFactory.getLogger(SimpleRenjinWithDBITest.class);

   /**
    * here we see a nasty stacktrace being printed due to
    * java.lang.NoClassDefFoundError: org/renjin/sexp/Environment$EmptyEnv
    * but the code still works
    */
   @Test
   public void testSimpleDBI() {
      RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = factory.getScriptEngine();
      StringBuilder str = new StringBuilder()
         .append("library('org.renjin.cran:DBI')\n")
         .append("library('org.renjin:hamcrest')\n")
         .append("a <- 24; print(a)");
      try {
         engine.eval(str.toString());
      } catch (ScriptException | EvalException e) {
         System.err.println("Failed to execute script");
         e.printStackTrace();
         fail(e.toString());
      }
   }


   /**
    * ...but as seen here, as soon as we try to
    * do something involving graphics it fails.
    */
   @Test
   public void testPlotToFile() {
      String script = "library(grDevices)\n" +
         "library(graphics)\n" +
         "library('org.renjin.cran:DBI')\n" +
         "\n" +
         "sysname <- tolower(Sys.info()['sysname'])\n" +
         "if(startsWith(sysname, 'windows')) {\n" +
         "    fileName <- 'c:/tmp/svgplot.svg'\n" +
         "} else {\n" +
         "    fileName <- '/tmp/svgplot.svg'\n" +
         "}\n" +
         "\n" +
         "# First plot\n" +
         "svg(fileName)\n" +
         "plot(sin, -pi, 2*pi)\n" +
         "\n" +
         "dev.off()\n" +
         "stopifnot(file.exists(fileName))";


      RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = factory.getScriptEngine();
      System.out.println("Running script: \n" + script);
      try {
         engine.eval(script);
      } catch (ScriptException | EvalException e) {
         System.err.println("Failed to execute script");
         e.printStackTrace();
         fail(e.toString());
      }
   }
}
