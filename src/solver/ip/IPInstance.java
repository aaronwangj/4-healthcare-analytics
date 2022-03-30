package solver.ip;

import ilog.cplex.*;
import solver.lp.IloNumExpr;
import solver.lp.IloNumVar;
import ilog.concert.*;

import java.lang.Math;
import java.util.Arrays;

public class IPInstance
{
  // IBM Ilog Cplex Solver 
  IloCplex cplex;
	
  int numTests;			// number of tests
  int numDiseases;		// number of diseases
  double[] costOfTest;  // [numTests] the cost of each test
  int[][] A;            // [numTests][numDiseases] 0/1 matrix if test is positive for disease
  int[][][] absDiffTDD;   // absDiffTDD[numTests][numDiseases][numDiseases] 0/1 matrix of the absolute difference of the 2 disease ([i][j][k]=[i][k][j])
  
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
    
    for(int t=0; t < numTests; t++) { // preprocessing (upper triangle excluding diagnol used)
        for(int d1 =0; d1 < numDiseases; d1++) {
        	for(int d2 =0; d2 < numDiseases; d2++) {
        		this.absDiffTDD[t][d1][d2] = Math.abs(A[t][d1]- A[t][d2]);
        		this.absDiffTDD[t][d2][d1] = Math.abs(A[t][d1]- A[t][d2]);
        	}
        }
       
    }
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
  
  public void solve() {
	  // I. Variables: binary variables of whether to include a test
	  IloNumVar[] testVars = cplex.numVarArray(numTests, 0, 1, IloNumVarType.Float);
	  // II. Constraints: for any pair of diseases, have at least one test that can differentiate them
	  for(int d1=0; d1 < numDiseases; d1++) {
		  for(int d2=d1+1; d2<numDiseases; d2++) { // d1<d2
			  IloNumExpr testDiff = cplex.constant(0);
			  for(int t=0; t<numTests; t++) {
				  testDiff = cplex.sum(testDiff, cplex.prod(testVars[t],absDiffTDD[t][d1][d2]));
			  }
			  cplex.add(cplex.ge(testResDiff, 0)); // absolute difference is nonnegative
		  }
	  }
	  // III. Objective: minimizes total cost
	  IloNumExpr testCost = cplex.constant(0);
	  for(int t=0; t<numTests; t++) {
		  testCost = cplex.sum(testCost, cplex.prod(testVars[t],costOfTest[t]));
	  }
	  cplex.addMinimize(testCost);

  }
}

















