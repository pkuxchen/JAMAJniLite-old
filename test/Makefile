JAVA_HOME = ..

JAMAJNILITE_CLASS = JAMAJniLiteTest.class
JAMAJNILITE_RUN = JAMAJniLite_test_result.txt
BLAS_CLASS = JAMAJniLiteExamplesBLAS.class
BLAS_RUN = BLAS_examples_result.txt
LAPACK_CLASS = JAMAJniLiteExamplesLAPACK.class
LAPACK_RUN = LAPACK_examples_result.txt
EXAMPLE_CLASS = JAMAJniLiteExamples.class
EXAMPLE_RUN = JAMAJniLite_examples_result.txt

all:  $(JAMAJNILITE_CLASS) $(JAMAJNILITE_RUN) $(BLAS_CLASS) $(BLAS_RUN) $(LAPACK_CLASS) $(LAPACK_RUN) $(EXAMPLE_CLASS) $(EXAMPLE_RUN)

$(JAMAJNILITE_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniLiteTest.java 

$(JAMAJNILITE_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniLiteTest > JAMAJniLite_test_result.txt

$(BLAS_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniLiteExamplesBLAS.java 

$(BLAS_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniLiteExamplesBLAS > BLAS_examples_result.txt

$(LAPACK_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniLiteExamplesLAPACK.java 

$(LAPACK_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniLiteExamplesLAPACK > LAPACK_examples_result.txt

$(EXAMPLE_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniLiteExamples.java 

$(EXAMPLE_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniLiteExamples > JAMAJniLite_examples_result.txt

clean: 
	rm -f *.txt *.class *.out *.serial
