package se.alipsa.javarc;

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