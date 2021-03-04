


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
#Source Code Structure
```
earth-storage
    ├── README.md
    ├── analysis-spark          空间分析的库
    ├── cli-wrapper             有关命令行的库
    ├── common                  数据格式读写、地理信息模型、计算几何、数据统计、命令行、杂项工具等基础类库
    ├── compute-spark           空间计算的库
    ├── dggs                    地球离散网格的库        
    ├── indexing                有关空间索引的基础库 （OLAP-Data Indexing to Optimize Data Access.
    ├── pom.xml                 
    ├── preparation-spark       有关数据预处理的库（OLAP-Data Preparation & Integretation）
    └── query-spark             有关数据查询的库  （OLAP-Data Query ）
```
##Module 'common'

##Model 'dggs'

##'MIn your code file, add codes to test if beast can be used: 
        import edu.ucr.cs.bdlab.beast.AppOptions;