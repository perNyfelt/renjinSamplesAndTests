package se.alipsa.json;

import org.junit.jupiter.api.Test;
import org.renjin.eval.EvalException;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.SEXP;

import java.io.StringWriter;
import javax.script.ScriptException;

/**
 * Related to issue <a href="https://github.com/bedatadriven/renjin/issues/518"/>
 */
public class JsonTest {

  static RenjinScriptEngine engine;

  @Test
  public void testJsonLite() throws ScriptException {
      RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
      engine = factory.getScriptEngine();
      engine.eval("saveRDS(c(1,2,3), 'foo.rds'); robj <- readRDS('foo.rds')");

      SEXP jobj = (SEXP)engine.eval("library('org.renjin.cran:jsonlite'); toJSON(robj)");
      System.out.println(jobj.asString());
  }

  @Test
  public void testJsonLiteAltOutput() throws ScriptException {
    RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
    engine = factory.getScriptEngine();
    engine.eval("foo <- c(1,2,3); saveRDS(foo, 'foo.rds'); robj <- readRDS('foo.rds')");

    engine.eval("library('org.renjin.cran:jsonlite'); jobj <- toJSON(robj)");

    StringWriter outputWriter = new StringWriter();
    engine.getContext().setWriter(outputWriter);
    String pcmd = "print(jobj)";
    engine.eval(pcmd);
    String result = outputWriter.toString();
    System.out.println(result);
  }

  @Test
  public void testRjson() throws ScriptException {
    RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
    engine = factory.getScriptEngine();
    engine.eval("saveRDS(c(1,2,3), 'foo.rds'); robj <- readRDS('foo.rds')");

    SEXP jobj = (SEXP)engine.eval("library('org.renjin.cran:rjson'); toJSON(robj)");
    System.out.println(jobj.asString());
  }

  @Test
  public void testJsonify() throws ScriptException {
    RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
    engine = factory.getScriptEngine();
    engine.eval("saveRDS(c(1,2,3), 'foo.rds'); robj <- readRDS('foo.rds')");

    SEXP jobj = null;
    try {
      jobj = (SEXP)engine.eval("library('org.renjin.cran:jsonify'); to_json(to_json)");
      System.out.println(jobj.asString());
    } catch (EvalException e) {
      e.printStackTrace();
    }
  }

}