package se.alipsa.charsequence;

import java.net.URL;

public class ScriptLoader {

  public static URL getScript() {
    return ScriptLoader.class.getResource("/Charsequence.R");
  }
}