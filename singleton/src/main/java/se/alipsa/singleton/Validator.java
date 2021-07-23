package se.alipsa.singleton;

public class Validator {

  final static Validator INSTANCE = new Validator();
  
  public static Validator getInstance() {
    return INSTANCE;
  }
  
  public boolean validate(String str) {
    System.out.println("Validating an object of class " + (str == null ? null : str.getClass()));
    return str != null && !str.isEmpty();
  }

  // If we pass in NULL from R Renjin will not be able to find the right method so
  // we must overload the method to handle null values
  public boolean validate(Object str) {
    return validate(str == null ? null : String.valueOf(str));
  }
}