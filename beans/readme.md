# Handling java beans in Renjin

As descriped in the [Renjin documentation](http://docs.renjin.org/en/latest/importing-java-classes-in-r-code.html#bean-classes),
Renjin provides som "secret sauce" to make access 
to java beans from R more natural i.e.

acessors (get methods) can be accessed "directly". 
A Java Bean class "Student" that has a private member variable called name
and a `public setName(String name)` as well as a `public String getName()`
defined can be constructed and accessed as follows:

```R
import(com.acme.Student)
bob <- Student$new(name = "Bob")

print(paste("The students name is", bob$name))
```

This example uses beans in a package similar to the above and exposes 
usage of this in an R API. 