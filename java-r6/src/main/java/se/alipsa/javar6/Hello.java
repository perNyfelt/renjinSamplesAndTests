package se.alipsa.javar6;

public class Hello {

  private String name;
  
  public String sayHello() {
    return "Hello " + name;
  }
  
  public Hello setName(String name) {
    this.name = name;
    return this;
  }

}