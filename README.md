# MAYGEN

Copyright 2021 Mehmet Aziz Yirik

## Introduction

MAYGEN is an open source chemical structure generator based on the orderly graph generation method. The principles of this method were outlined in the works by Grund et al. [1]. The theoretical basis and the outlines of the functions can be found in [1,2].  
MAYGEN takes a molecular formula (such as C10H16O) as input and generates all constitutional isomers of this formula, i. e. all non-isomorphic molecules that can be constructed with the set of atoms in the input formula. For the case of C10H16O, for example, there are 452,458 non-identical molecules.  
[[resources/C10H16O.png|Twelve out of 400k+ isomers of C10H16O]]

## Download jar File

The main class of the current jar file is MAYGEN class. The file can be downloaded from [here](https://github.com/MehmetAzizYirik/MAYGEN/releases/tag/V1.0)

## Download Source Code

You can download the source code as a ZIP file from the landing page of this repository. 
Alternatively, you can clone the repository using GIT. For more information [set-up-git](https://help.github.com/articles/set-up-git/ )

To download MAYGEN source code:

```
$ git clone https://github.com/MehmetAzizYirik/MAYGEN.git
```
## Compiling

To compile MAYGEN, Apache Maven and Java 1.8 (or later) are required.
```
MAYGEN/$ mvn package
```
This command will create jar file named as "MAYGEN-jar-with-dependencies" under the target folder.

## Usage

MAYGEN.jar can be run from command line with the specified arguments. An example command is given below.

The definitions of the arguments are given below:

```
usage: java -jar MAYGEN.jar -f <arg> [-v] -d <arg>

Generates molecular structures for a given molecular formula. The input is 
a molecular formula string. For example 'C2OH4'. Besides molecular formula, the
directory is needed to be specified for the output file.

 -f,--molecularFormula <arg>   Molecular formula as a string (required)
 
 -v,--verbose                  Print messages about structure generation
 
 -d,--filedir <arg>            Creates and store the output txt file 
                               in the directory 

Please report issues at https://github.com/MehmetAzizYirik/MAYGEN
```

```
java -jar MAYGEN.jar -f C2OH4 -v -d C:\Users\UserName\Desktop\
```

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/MehmetAzizYirik/MAYGEN/blob/master/LICENSE) file for details

## Authors

 - Mehmet Aziz Yirik - [MehmetAzizYirik](https://github.com/MehmetAzizYirik)
 
## Acknowledgements
![YourKit](https://camo.githubusercontent.com/97fa03cac759a772255b93c64ab1c9f76a103681/68747470733a2f2f7777772e796f75726b69742e636f6d2f696d616765732f796b6c6f676f2e706e67)

The developer uses YourKit to profile and optimise code.

YourKit supports open source projects with its full-featured Java Profiler. YourKit, LLC is the creator of YourKit Java Profiler and YourKit .NET Profiler, innovative and intelligent tools for profiling Java and .NET applications.

![cdk](https://github.com/MehmetAzizYirik/HMD/blob/master/cdk.png)

This project relies on the Chemistry Development Project (CDK), hosted under [CDK GitHub](http://cdk.github.io/). Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is available though our [Github site](http://cdk.github.io/cdk/).

## References

1- Grund, R. and Müller, R., 1995. Konstruktion molekularer Graphen mit gegebenen Hybridisierungen und überlappungsfreien Fragmenten. Lehrstuhl II für Mathematik.

2- Kerber, A., Laue, R., Meringer, M., Rücker, C. and Schymanski, E., 2013. Mathematical chemistry and chemoinformatics: structure generation, elucidation and quantitative structure-property relationships. Walter de Gruyter.


