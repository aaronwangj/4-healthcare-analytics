package solver.ip;

import java.io.File;
import java.io.FileNotFoundException;

import java.util.Scanner;

/**
 * File Format
 * #Tests (i.e., n)
 * #Diseases (i.e., m)
 * Cost_1 Cost_2 . . . Cost_n
 * A(1,1) A(1,2) . . . A(1, m)
 * A(2,1) A(2,2) . . . A(2, m)
 * . . . . . . . . . . . . . .
 * A(n,1) A(n,2) . . . A(n, m)
 */
public class DataParser
{
  public static IPInstance parseIPFile(String fileName)
  {
    IPInstance instance = new IPInstance();
    try
    {
      Scanner read = new Scanner(new File(fileName));
      
      int numTests = read.nextInt();  // n
      int numDiseases = read.nextInt();  // m
      
      double[] costOfTest = new double[numTests];
      for(int i=0; i < costOfTest.length; i++)
        costOfTest[i] = read.nextDouble();
      
      int[][] A = new int[numTests][numDiseases];
      for(int i=0; i < numTests; i++)
        for(int j=0; j < numDiseases; j++)
          A[i][j] = read.nextInt();
      
      // Initialize the instance
      instance.init(numTests, numDiseases, costOfTest, A);
    }
    catch (FileNotFoundException e)
    {
      System.out.println("Error: file not found " + fileName);
    }
    return instance;
  }
}
