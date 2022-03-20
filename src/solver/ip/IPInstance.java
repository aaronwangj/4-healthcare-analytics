package solver.ip;

import ilog.cplex.*;
import ilog.concert.*;

import java.util.Arrays;

public class IPInstance
{
  // IBM Ilog Cplex Solver 
  IloCplex cplex;
	
  int numTests;			// number of tests
  int numDiseases;		// number of diseases
  double[] costOfTest;  // [numTests] the cost of each test
  int[][] A;            // [numTests][numDiseases] 0/1 matrix if test is positive for disease
  
  public IPInstance()
  {
	super();
  }
  
  void init(int numTests, int numDiseases, double[] costOfTest, int[][] A)
  {
    assert(numTests >= 0) : "Init error: numtests should be non-negative " + numTests;
    assert(numDiseases >= 0) : "Init error: numtests should be non-negative " + numTests;
    assert(costOfTest != null) : "Init error: costOfTest cannot be null";
    assert(costOfTest.length == numTests) : "Init error: costOfTest length differ from numTests" + costOfTest.length + " vs. " + numTests;
    assert(A != null) : "Init error: A cannot be null";
    assert(A.length == numTests) : "Init error: Number of rows in A differ from numTests" + A.length + " vs. " + numTests;
    assert(A[0].length == numDiseases) : "Init error: Number of columns in A differ from numDiseases" + A[0].length + " vs. " + numDiseases;
    
    this.numTests = numTests;
    this.numDiseases = numDiseases;
    this.costOfTest = new double[numTests];
    for(int i=0; i < numTests; i++)
      this.costOfTest[i] = costOfTest[i];
    this.A = new int[numTests][numDiseases];
    for(int i=0; i < numTests; i++)
      for(int j=0; j < numDiseases; j++)
        this.A[i][j] = A[i][j];
  }
  
  public String toString()
  {
	StringBuffer buf = new StringBuffer();
	buf.append("Number of tests: " + numTests + "\n");
	buf.append("Number of diseases: " + numDiseases + "\n");
	buf.append("Cost of tests: " + Arrays.toString(costOfTest) + "\n");
	buf.append("A:\n");
	for(int i=0; i < numTests; i++)
		buf.append(Arrays.toString(A[i]) + "\n");
	return buf.toString();
  }
}