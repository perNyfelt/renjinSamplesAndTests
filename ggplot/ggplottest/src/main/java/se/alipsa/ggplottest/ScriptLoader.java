package se.alipsa.ggplottest;

import java.net.URL;

public class ScriptLoader {

  public static URL getScript() {
    return ScriptLoader.class.getResource("/Ggplottest.R");
  }
}