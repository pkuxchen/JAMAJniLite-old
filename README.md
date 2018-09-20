# JAMAJniLite
JAMAJniLite is a sibling package of JAMAJni. It is a JAVA package providing a java interface for lapack and blas library and using the classes defined by JAMA Package. It's a tool for numerical linear algebra computations. JAMAJniLite calls lapack and blas libraries which require previous installation (Detailed installation instructions are listed below). You can also use openblas library(which is an optimized blas and lapack library by changing the Makefile).


Build Instructions
------------------

* JAMAJniLite requires the installation of lapack and blas. (The command in Ubuntu is "sudo apt-get install libblas-dev liblapack-dev"). 

* To compile the package, enter src directory and execute "make". Notice that you may have to change the extension of generated libraries in the Makefile based on your operating system. On OS X you have to change all the extensions of dynamic library to .dylib while on Linux the corresponding extensions are .so or .a. 

* To clean generated file, type “make clean” on the command line.  

* If you would like to use openblas library, you can install openblas by command "sudo apt-get install libopenblas-dev". Then you need to replace all the "-llapack" and "-lblas" by the "-lopenblas" in Makefile.

Running the tests
-----------------
* For testing, enter test directory and execute “make” . If you want to clean testing results and all class files, type "make clean". 

* There are four test files. The "JAMAJniLiteTest.java" will test all the methods in JAMAJniLite and report the errors. The "JAMAJniLiteExamples.java" will provide specific examples for basic linear algebra operations. It can clearly show you how to use methods defined in JAMAJniLite. If you are interested in how to use functions in blas and lapack libraries to do matrix operations, the "JAMAJniLiteExamplesBLAS.java" and "JAMAJniLiteExamplesLAPACK.java" will give you specific examples. However, it is not necessary to go into blas and lapack if you just want to be a user of JAMAJniLite.

Notes
---------
This package is intended for some basic problems we encounter when solving linear algebra problems. So we only include several most basic and widely used routines of the blas and lapack library.


Source Repository
-----------------
JAMAJniLite's source-code repository is hosted here on GitHub.


Authors
---------

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Xiang Chen (maintainer)| pkuchenxiang@pku.edu.cn   | Visiting student, Department of Biostatistics  UCLA|
| Lu Zhang | lu.zhang@ucla.edu    | PhD student, Department of Biostatistics UCLA  |                            
| Sudipto Banerjee | sudipto@ucla.edu   | Professor, Department of Biostatistics  UCLA |
<!--- --->
                             


Licensing
---------
JAMAJniLite is licensed under the Creative Commons Attribution License. 



