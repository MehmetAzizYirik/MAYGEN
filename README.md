<img src="docs/logo.png" alt="drawing" width="100" align = "right"/> [![DOI](https://zenodo.org/badge/359922963.svg)](https://zenodo.org/badge/latestdoi/359922963)
[![Maintainability Rating](https://sonarcloud.io/api/project_badges/measure?project=MehmetAzizYirik_MAYGEN&metric=sqale_rating)](https://sonarcloud.io/summary/overall?id=MehmetAzizYirik_MAYGEN)
[![Reliability Rating](https://sonarcloud.io/api/project_badges/measure?project=MehmetAzizYirik_MAYGEN&metric=reliability_rating)](https://sonarcloud.io/summary/overall?id=MehmetAzizYirik_MAYGEN)
[![Security Rating](https://sonarcloud.io/api/project_badges/measure?project=MehmetAzizYirik_MAYGEN&metric=security_rating)](https://sonarcloud.io/summary/overall?id=MehmetAzizYirik_MAYGEN)
[![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=MehmetAzizYirik_MAYGEN&metric=ncloc)](https://sonarcloud.io/summary/overall?id=MehmetAzizYirik_MAYGEN)
[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=MehmetAzizYirik_MAYGEN&metric=bugs)](https://sonarcloud.io/summary/overall?id=MehmetAzizYirik_MAYGEN)
[![Build](https://github.com/MehmetAzizYirik/MAYGEN/actions/workflows/maven.yml/badge.svg)](https://github.com/MehmetAzizYirik/MAYGEN/actions/workflows/maven.yml)
[![CodeQL](https://github.com/MehmetAzizYirik/MAYGEN/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/MehmetAzizYirik/MAYGEN/actions/workflows/codeql-analysis.yml)
# MAYGEN - A chemical structure generator for constitutional isomers based on the orderly generation principle

Copyright 2021 Mehmet Aziz Yirik

## Introduction

MAYGEN is an open source chemical structure generator based on the orderly graph generation method. The principles of this method were outlined in the [MAYGEN article](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00529-9)[1].
MAYGEN takes a molecular formula (such as C<sub>10</sub>H<sub>16</sub>O) as input and generates all constitutional isomers of this formula, i. e. all non-isomorphic molecules that can be constructed with the set of atoms in the input formula. For the case of C<sub>10</sub>H<sub>16</sub>O, for example, there are 452,458 non-identical molecules. Here are 12 out of those.

<p align="center">
  <img src=/resources/C10H16O.png />
</p>

As can be seen from these examples, MAYGEN makes no assumptions on chemical stability. In particular in small ring systems, this may lead to unlikely structures, such as C=1C=C1. 

We benchmarked MAYGEN V.1.4 against the current state-of-the-art, but [closed-source structure generator MOLGEN 5.0](http://www.molgen.de) from the University of Bayreuth as well as against the [Parallel Molecule Generator (PMG)](https://sourceforge.net/projects/pmgcoordination/)[2], the fastest available open source structure generator. Since PMG can be run in multi-threaded mode, the benchmark was performed in single-threaded mode for algorithmic comparability. For randomly selected 50 formulae, MAYGEN was in average 3 times slower than MOLGEN but 47 times faster than PMG. For some formulae, PMG could not generate isomers. These are shown by gaps on the its plot.

<p align="center">
  <img src=/resources/Benchmarking.png />
</p>

## Download jar File

Executable JAR files can be downloaded from [the release page](https://github.com/MehmetAzizYirik/MAYGEN/releases)

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
This command will create jar file named as "MAYGEN-1.8" under the target folder.

## Usage

MAYGEN-1.8.jar can be run from command line with the specified arguments. An example command is given below.

The definitions of the arguments are given below:

```
usage: java -jar MAYGEN-1.8.jar [-f <arg>] [-fuzzy <arg>] [-setElements]
       [-v] [-t] [-o <arg>] [-b] [-m] [-smi] [-sdf] [-sdfCoord] [-h]

Generates molecular structures for a given molecular formula.
The input is a molecular formula string.

For example 'C2OH4'.

If user wants to store output file in a specific directory, that is needed
to be specified. It is also possible to generate SMILES instead of an SDF
file, but it slows down the generation time. For this, use the '-smi'
option.

 -f,--formula <arg>               formula

 -fuzzy,--fuzzyFormula <arg>      fuzzy formula

 -setElements,--settingElements   User defined valences

 -v,--verbose                     print message

 -t,--tsvoutput                   Output formula, number of structures and
                                  execution time in CSV format. In
                                  multithread, the 4th column in the
                                  output is the number of threads.

 -o,--outputFile <arg>            Store output file

 -b,--boundaryConditions          Setting the boundary conditions option

 -m,--multithread                 Use multi thread

 -smi,--SMILES                    Output in SMILES format

 -sdf,--SDF                       Output in SDF format

 -sdfCoord,--coordinates          Output in SDF format with atom
                                  coordinates

 -h,--help                        Displays help message

Please report issues at https://github.com/MehmetAzizYirik/MAYGEN
```

```
java -jar MAYGEN-1.8.jar -f C2OH4 -v -t -o C:\Users\UserName\Desktop\

java -jar MAYGEN-1.8.jar -fuzzy C[2-5]O2H[4-8] -v -t -o C:\Users\UserName\Desktop\

java -jar MAYGEN-1.8.jar -f N(val=4)6H6 -setElements -v -t -o C:\Users\UserName\Desktop\
```
## Webservice

A [webservice](https://www.cheminfo.org/?viewURL=https%3A%2F%2Fcouch.cheminfo.org%2Fcheminfo-public%2F7cead2ae3da71090cb17baa3856ea38b%2Fview.json) is also developed for MAYGEN software for easy usage and educational purposes. The documentation for the webservice is given [here](https://maygen.cheminfo.org/documentation/static/index.html). 

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

1- Yirik, M.A., Sorokina, M. & Steinbeck, C. MAYGEN: an open-source chemical structure generator for constitutional isomers based on the orderly generation principle. J Cheminform 13, 48 (2021). https://doi.org/10.1186/s13321-021-00529-9

2- Jaghoori MM, Jongmans SS, De Boer F, Peironcely J, Faulon JL, Reijmers T, Hankemeier T. PMG: multi-core metabolite identification. Electronic Notes in Theoretical Computer Science. 2013 Dec 25;299:53-60.
  
  


