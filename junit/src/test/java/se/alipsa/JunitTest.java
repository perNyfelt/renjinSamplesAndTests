package se.alipsa;

import org.junit.Test;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;

import javax.script.ScriptException;

public class JunitTest {

   @Test
   public void TestJunit() throws ScriptException {
      RenjinScriptEngineFactory renjinScriptEngineFactory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = renjinScriptEngineFactory.getScriptEngine();
      engine.eval("a <- 24; print(a)");
   }

}
