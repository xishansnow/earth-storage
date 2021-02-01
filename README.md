


#Prerequisite:

  Beast is a package for spatial big data distribution storage & parallel processing.  For its version on maven repository isn't the latest. It may good to compile it yourself.
  (1) Compile Beast packages        
        You can find beast on https://bitbucket.org/eldawy/beast/src/master/. Clone it and compile with instruction in README.md.
  (2) Install Beast jars to local maven repository.
        Go to the directory which contain beast-*-*.jar, and use 'mvn install' command like below.
        mvn install:install-file -Dfile=beast-uber-spark-0.9.0-SNAPSHOT.jar -DgroupId=edu.ucr.cs.bdlab -DgroupId=edu.ucr.cs.bdlab -DartifactId=beast -Dversion=0.9.0 -Dpackaging=jar
    (3) Test Beast 
        In your project, edit POM.xml, add dependency like that:
'''
            <dependencies>
                <dependency>
                    <groupId>edu.ucr.cs.bdlab</groupId>
                    <artifactId>beast</artifactId>
                    <version>0.9.0</version>
                </dependency>
            </dependencies>
''''
        In your code file, add codes to test if beast can be used: 
        import edu.ucr.cs.bdlab.beast.BeastOptions;