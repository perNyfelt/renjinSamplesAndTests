package se.alipsa.renjin.examples.helloworld;

import org.junit.jupiter.api.Test;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.SEXP;

import javax.script.ScriptException;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HelloWorldTest {

  @Test
  public void testHelloWorld() throws ScriptException {
    RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
    RenjinScriptEngine engine = factory.getScriptEngine();

    // See http://docs.renjin.org/en/latest/library/capture.html for how to capture results from the script
    SEXP result = (SEXP)engine.eval("paste('Hello', 'World')");
    assertEquals("Hello World", result.asString());
  }
}
