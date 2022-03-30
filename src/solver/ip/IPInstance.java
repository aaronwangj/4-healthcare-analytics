package solver.ip;

import ilog.cplex.*;
//import solver.lp.IloNumExpr;
//import solver.lp.IloNumVar;
import ilog.concert.*;

import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

public class IPInstance
{
  // IBM Ilog Cplex Solver 
  IloCplex cplex;
	
  int numTests;			// number of tests
  int numDiseases;		// number of diseases
  double[] costOfTest;  // [numTests] the cost of each test
  int[][] A;            // [numTests][numDiseases] 0/1 matrix if test is positive for disease
  int[][][] absDiffTDD;   // absDiffTDD[numTests][numDiseases][numDiseases] 0/1 matrix of the absolute difference of the 2 disease ([i][j][k]=[i][k][j])
  
  double bestObjective = -1;
  int[] incumbentSol = new int[numTests];
  IloNumVar[] testVars;  // for access in branchNBound
  
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
    
    this.absDiffTDD = new int[numTests][numDiseases][numDiseases];
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
  
  public double solve() throws IloException{
	 try{
		 cplex = new IloCplex();
		  
		  // I. Variables: binary variables of whether to include a test
		  testVars = cplex.numVarArray(numTests, 0, 1, IloNumVarType.Float);
		  // II. Constraints: for any pair of diseases, have at least one test that can differentiate them
		  for(int d1=0; d1 < numDiseases; d1++) {
			  for(int d2=d1+1; d2<numDiseases; d2++) { // d1<d2
				  IloNumExpr testDiff = cplex.constant(0);
				  for(int t=0; t<numTests; t++) {
					  testDiff = cplex.sum(testDiff, cplex.prod(testVars[t],absDiffTDD[t][d1][d2]));
				  }
				  cplex.add(cplex.ge(testDiff, 0)); // absolute difference is nonnegative
			  }
		  }
		  // III. Objective: minimizes total cost
		  IloNumExpr testCost = cplex.constant(0);
		  for(int t=0; t<numTests; t++) {
			  testCost = cplex.sum(testCost, cplex.prod(testVars[t],costOfTest[t]));
		  }
		  cplex.addMinimize(testCost);
		  
		  // IV. Solving
		Set<Integer> unconsTestVarIdx = new HashSet<Integer>();
		for (int t = 0; t<numTests; t++) {
			unconsTestVarIdx.add(t);
		}
		// optimization TODO: List<Integer> testVarBranchVal = new ArrayList<Integer>(numTests); // idx: 0 or 1
		  branchNBound(cplex, testVarIdx);
		  if( bestObjective == -1) {
			  System.out.println("No solution found");
		  } else {
			  System.out.println("Solution:");
		      for (int t = 0; t<numTests; t++) {
		         System.out.print(incumbentSol[t]);
		      }
		      System.out.println();
		  }
		  return bestObjective;
	 } catch(IloException e){
		  System.out.println("Error " + e);
	 }
  }
  
  public void branchNBound(IloCplex cplex, Set<Integer> unconsTestVarIdx) {
  /* Recursive implementation of the branch & bound algorithm for solving IP problem with successive LP solving. bestObjective and incumbentSol altered upon return.
   * cplex: LP solver instance whose constraints are modified through passes and whose results are loaded for checking.
   * unconsTestVarIdx: indices of test variables that are not yet constrained. Vars not on this are assigned vals that can be checked through cplex.
   */
	  // NOTE: termination condition of exhausting tree falls into 1. 2. or 3a. (guaranteed to be integral)
	  if (!cplex.solve()) return;  // 1.  Infeasible: prune this branch (has no solution) and backtrack
	  if (cplex.getObjValue() >= bestObjective && bestObjective != -1) return; // 2. Feasible + not dominant: prune and backtrack
	  // 3. Feasible + dominant
	  double[] testVals = new double[numTests]; testVals = cplex.getValues(testVars);
	  int[] testValsInt = new int[numTests]; Boolean isIntSol = isIntegralSol(testVals, testValsInt);
	  if(isIntSol) { // 3a. Integral: optimal for branch (may be entire tree if at root) -> update + backtrack
		  bestObjective = cplex.getObjValue();
		  incumbentSol = testValsInt;
		  return;
	  }  else { // 3b. Not integral: keep searching (guaranteed unconsTestVarIdx.size() > 0)->select variable to branch on (enforcing integral value)
		  int branchTestVarIdx = unconsTestVarIdx.iterator().next(); // TODO: optimization - which variable to pick (fractional in exisitng solution / best first)
		  unconsTestVarIdx.remove(branchTestVarIdx);
		// TODO: optimization: which value to pick (0 or 1)
		  IloRange branchingCon = myCplex.addEq(testVars[branchTestVarIdx],cplex.constant(1)); cplex.add(branchingCon);
		  branchNBound(cplex, unconsTestVarIdx);
		  cplex.remove(branchingCon); // delete
		  branchingCon = myCplex.addEq(testVars[branchTestVarIdx], cplex.constant(0)); cplex.add(branchingCon);
		  branchNBound(cplex, unconsTestVarIdx);
		  cplex.remove(branchingCon); // delete
		  unconsTestVarIdx.add(branchTestVarIdx); // May be unnecessary
	  }
  }
	  
  
  public boolean isIntegralSol(double[] testVals, int[] testValsInt) {
	  /* Returns true if all values in testVals are integers. testValsInt populated (passed by reference) */
	  for (int t = 0; t<numTests; t++) {
		  if ( (testVals[t]!= 1.0) && (testVals[t]!=0.0) ) {
			  return false; // has fractional value
		  }
		  testValsInt[t] = (int)(testVals[t]);
	  }
	  return true;
  }
}

















