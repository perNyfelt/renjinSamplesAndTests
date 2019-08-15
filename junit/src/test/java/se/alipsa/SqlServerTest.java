package se.alipsa;

import org.junit.Rule;
import org.junit.Test;
import org.renjin.eval.EvalException;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.testcontainers.containers.MSSQLServerContainer;

import javax.script.ScriptException;

import static org.junit.Assert.fail;

public class SqlServerTest {

   Logger log = LoggerFactory.getLogger(SqlServerTest.class);

   @Rule
   public MSSQLServerContainer mssqlserver = new MSSQLServerContainer();


   @Test
   public void testSqlServer() {
      RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = factory.getScriptEngine();
      StringBuilder str = new StringBuilder()
         .append("library('org.renjin.cran:DBI')\n")
         .append("library('org.renjin:hamcrest')\n")
         .append("a <- 24; print(a)");

      engine.put("mssqlserver", mssqlserver);
      try {
         engine.eval(str.toString());
      } catch (ScriptException | EvalException e) {
         System.err.println("Failed to execute script");
         e.printStackTrace();
         fail(e.toString());
      }
   }
}
