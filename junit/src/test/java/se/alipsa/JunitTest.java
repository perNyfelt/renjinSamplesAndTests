package se.alipsa;

import org.junit.Test;
import org.renjin.eval.EvalException;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;

import javax.script.ScriptException;

import static org.junit.Assert.fail;

public class JunitTest {

   @Test
   public void TestJunit() {
      RenjinScriptEngineFactory renjinScriptEngineFactory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = renjinScriptEngineFactory.getScriptEngine();
      try {
         engine.eval("a <- 24; print(a)");
      } catch (ScriptException | EvalException e) {
         System.err.println("Failed to execute script");
         e.printStackTrace();
         fail(e.toString());
      }
   }

}
