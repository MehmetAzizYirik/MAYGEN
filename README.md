# MAYGEN - A chemical structure generator for constitutional isomers based on the orderly generation principle

Copyright 2021 Mehmet Aziz Yirik

## Introduction

MAYGEN is an open source chemical structure generator based on the orderly graph generation method. The principles of this method were outlined in the works by Grund et al. [1]. The theoretical basis and the outlines of the functions can be found in [1,2].  
MAYGEN takes a molecular formula (such as C<sub>10</sub>H<sub>16</sub>O) as input and generates all constitutional isomers of this formula, i. e. all non-isomorphic molecules that can be constructed with the set of atoms in the input formula. For the case of C<sub>10</sub>H<sub>16</sub>O, for example, there are 452,458 non-identical molecules. Here are 12 out of those.

![Twelve out of 400k+ isomers of C10H16O](/resources/C10H16O.png)

As can be seen from these examples, MAYGEN makes no assumptions on chemical stability. In particular in small ring systems, this may lead to unlikely structures, such as C=1C=C1. 

We benchmarked MAYGEN against the current state-of-the-art, but [closed-source structure generator Molgen](http://www.molgen.de) from the University of Bayreuth as well as against the [Open Molecule Generator (OMG)](https://sourceforge.net/projects/openmg/)[3], the only available open source structure generator. We did not test against the parallelised version of OMG - called PMG - to maintain the focus on algorithmic rather than technological speed. The following plot compares the timings. The generated numbers  of structures are identical for MAYGEN and Molgen 3.5.

![A speed comparison of maygen against molgen. Maygen is consistently slower than molgen](/resources/maygen-molgen.png)

MAYGEN is consistently slower than Molgen, but gets close to its performance for pure carbohydrates (Factor 1.4 slower). For example, MAYGEN generates the 400 mio isomers of C<sub>13</sub>H<sub>8</sub> in 14h on a current unix machine, whereas Molgen takes 11h. 

![Results for the carbohydrates with oxygens](/resources/CarbohydratesWithOxygens.png)

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
usage: java -jar MAYGEN.jar -f <arg> [-v] [-t] -o <arg>

Generates molecular structures for a given molecular formula. The input is 
a molecular formula string. For example 'C2OH4'. Besides molecular formula, the
directory is needed to be specified for the output file.

 -f,--molecularFormula <arg>   Molecular formula as a string (required)
 
 -v,--verbose                  Print messages about structure generation
 
 -t,--tsvoutput                Output formula, number of structures and execution time in CSV format
 
 -o,--filename <arg>           Store output in given file 

Please report issues at https://github.com/MehmetAzizYirik/MAYGEN
```

```
java -jar MAYGEN.jar -f C2OH4 -v -d C:\Users\UserName\Desktop\
```

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/MehmetAzizYirik/MAYGEN/blob/main/LICENSE) file for details

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

3- Peironcely, J. E., Rojas-Chertó, M., Fichera, D., Reijmers, T., Coulier, L., Faulon, J.-L. & Hankemeier, T. OMG: Open Molecule Generator. J Cheminformatics 4, 21 (2012).
  
  


