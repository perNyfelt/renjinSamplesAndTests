package se.alipsa.refrenceclasses;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.CoreMatchers.*;

import org.junit.jupiter.api.Test;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.SEXP;

import java.io.IOException;
import java.io.InputStreamReader;
import javax.script.ScriptException;

public class ReferenceClassesTest {

  @Test
  public void testReferenceClasses() throws IOException, ScriptException {
    try (InputStreamReader in = new InputStreamReader(ScriptLoader.getScript().openStream())) {

      RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = factory.getScriptEngine();

      // See http://docs.renjin.org/en/latest/library/capture.html for how to capture results from the script
      SEXP result = (SEXP)engine.eval(in);
      // Do some assertions to ensure the script is working correctly
      assertThat(result, is(notNullValue()));

      SEXP s1 = (SEXP)engine.eval("s1 <- StandardAssignment$new()\n" +
          "s1$setAttribute('name', 'foo')");

      // The session is still the same so we can reference variables created in the previous eval
      SEXP nameVal = (SEXP)engine.eval("s1$getAttribute('name')");

      assertThat(nameVal.asString(), equalTo("foo"));

    }
  }

}