package se.alipsa;

import org.junit.Test;
import org.renjin.eval.Context;
import org.renjin.eval.EvalException;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.Environment;
import org.renjin.sexp.StringVector;

import javax.script.ScriptException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

public class JunitTest {

   /**
    * here we see a nasty stacktrace being printed due to
    * java.lang.NoClassDefFoundError: org/renjin/sexp/Environment$EmptyEnv
    * but the code still works
    */
   @Test
   public void testJunit() {
      RenjinScriptEngineFactory renjinScriptEngineFactory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = renjinScriptEngineFactory.getScriptEngine();
      try {
         engine.eval("a <- 'testJunit'; print(a)");
         Environment global = engine.getSession().getGlobalEnvironment();
         Context topContext = engine.getSession().getTopLevelContext();
         StringVector strVec = (StringVector)global.getVariable(topContext, "a");
         assertEquals("testJunit", strVec.asString());
      } catch (ScriptException | EvalException e) {
         System.err.println("Failed to execute script");
         e.printStackTrace();
         fail(e.toString());
      }
   }

}
