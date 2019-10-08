package se.alipsa.singleton;

public class Validator {

  final static Validator INSTANCE = new Validator();
  
  public static Validator getInstance() {
    return INSTANCE;
  }
  
  public boolean validate(String str) {
    return str == null ? false : true;
  }
}