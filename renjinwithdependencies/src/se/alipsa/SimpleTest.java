package se.alipsa;

import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.SEXP;

import javax.script.ScriptException;

public class SimpleTest {

   public static void main(String[] args) throws ScriptException {
      RenjinScriptEngineFactory renjinScriptEngineFactory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = renjinScriptEngineFactory.getScriptEngine();

      SEXP result = (SEXP)engine.eval("p <- 48; foo <- function(a, b) { return( a * b) }; bar <- foo(p, 2)");
      System.out.println("Result is " + result.asInt());
   }
}
