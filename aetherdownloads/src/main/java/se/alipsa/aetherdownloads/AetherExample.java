package se.alipsa.aetherdownloads;

import org.eclipse.aether.repository.RemoteRepository;
import org.renjin.aether.AetherFactory;
import org.renjin.aether.AetherPackageLoader;
import org.renjin.aether.ConsoleRepositoryListener;
import org.renjin.aether.ConsoleTransferListener;
import org.renjin.eval.Session;
import org.renjin.eval.SessionBuilder;
import org.renjin.script.RenjinScriptEngine;
import org.renjin.script.RenjinScriptEngineFactory;
import org.renjin.sexp.SEXP;

import java.util.ArrayList;
import java.util.List;
import javax.script.ScriptException;

/** This example was created to study the details of AetherPackageLoader downloads
 * The transfer listener output is somewhat confusing as it prints out a lot of TransferException stacktraces for
 * Each repo that does not have the artefact in question which leads you to think that something is wrong
 * but this is just how it is. Once resolved it just continues on to the next dependency and in the end everything works.
  */
public class AetherExample {

  public static void main(String[] args) throws ScriptException {

    List<RemoteRepository> repoList = new ArrayList<>();
    // THis will download maven-metadata.xml as maven-metadata-be-data-driven.xml
    repoList.add(new RemoteRepository
        .Builder("be-data-driven","default", "https://nexus.bedatadriven.com/content/groups/public/")
        .build());

    repoList.add(AetherFactory.mavenCentral());

    AetherPackageLoader loader = new AetherPackageLoader(AetherExample.class.getClassLoader(), repoList);
    loader.setRepositoryListener(new ConsoleRepositoryListener(System.out));
    loader.setTransferListener(new ConsoleTransferListener(System.out));

    ClassLoader classLoader = loader.getClassLoader();

    SessionBuilder builder = new SessionBuilder();
    Session session = builder
        .withDefaultPackages()
        .setPackageLoader(loader) // allows library to work without having to include in the pom
        .setClassLoader(classLoader) //allows imports in r code to work
        .build();

    RenjinScriptEngineFactory factory = new RenjinScriptEngineFactory();
    RenjinScriptEngine engine = factory.getScriptEngine(session);

    String code = "library('digest'); digest(paste0('secretSalt','9901011234'), algo='sha256', serialize=FALSE)";

    SEXP sexp = (SEXP)engine.eval(code);
    System.out.println(sexp);
  }
}
