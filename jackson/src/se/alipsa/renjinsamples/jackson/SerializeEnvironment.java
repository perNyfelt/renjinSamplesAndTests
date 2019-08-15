package se.alipsa.renjinsamples.jackson;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.renjin.eval.Session;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.Environment;

import javax.script.ScriptException;
import java.util.Iterator;

public class SerializeEnvironment {

   public static void main(String[] args) throws JsonProcessingException, ScriptException {
      RenjinScriptEngineFactory renjinScriptEngineFactory = new RenjinScriptEngineFactory();
      RenjinScriptEngine engine = renjinScriptEngineFactory.getScriptEngine();

      engine.eval("p <- 48; foo <- function(a, b) { return( a * b) }; bar <- foo(12, 2)");

      Session session = engine.getSession();
      Environment environment = session.getGlobalEnvironment();

      Iterator it = environment.getNames().iterator();
      while (it.hasNext()) {
         Object obj = it.next();
         System.out.println(obj);
      }
      ObjectMapper objectMapper = new ObjectMapper();
      //System.out.println(objectMapper.writeValueAsString(environment));
   }
}
