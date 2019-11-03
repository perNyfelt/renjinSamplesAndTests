package se.alipsa.ctsmrexample;

import java.net.URL;

public class ScriptLoader {

  public static URL getScript() {
    return ScriptLoader.class.getResource("/CtsmrExample.R");
  }
}