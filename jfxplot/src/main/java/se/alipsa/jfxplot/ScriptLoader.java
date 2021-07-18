package se.alipsa.jfxplot;

import java.net.URL;

public class ScriptLoader {

  public static URL getScript() {
    return ScriptLoader.class.getResource("/Jfxplot.R");
  }
}