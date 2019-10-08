package se.alipsa.singleton;

import java.net.URL;

public class ScriptLoader {

  public static URL getScript() {
    return ScriptLoader.class.getResource("/Singleton.R");
  }
}