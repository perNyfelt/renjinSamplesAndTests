package se.alipsa;

import org.junit.Rule;
import org.junit.Test;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.testcontainers.containers.MSSQLServerContainer;

import javax.script.ScriptEngine;
import javax.script.ScriptException;

public class SqlServerTest {

   Logger log = LoggerFactory.getLogger(SqlServerTest.class);

   //@Rule
   //public MSSQLServerContainer mssqlserver = new MSSQLServerContainer();


   @Test
   public void TestSqlServer() throws ScriptException {
      RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = factory.getScriptEngine();
      StringBuilder str = new StringBuilder()
      .append("library('org.renjin.cran:DBI')\n")
      //.append("library('se.alipsa:R2JDBC')\n")
      //  .append("library('org.renjin:hamcrest')\n");
      //engine.put("mssqlserver", mssqlserver)
         .append("a <- 24; print(a)");
      ;
      engine.eval(str.toString());
   }
}
