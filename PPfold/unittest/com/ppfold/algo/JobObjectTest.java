package com.ppfold.algo;

public class JobObjectTest {

    public int[][] matrix1 = new int[8][8];
    public int[][] matrix2= new int[8][8];

    
    public JobObjectTest(int n){
    	 for(int i = 1; i<matrix1.length; i++){
    		 for(int j = 1; j<matrix2.length; j++){
    			 matrix1[i][j] = n;
    			 matrix2[i][j] = n; 
    		 }
    	 }
    }
    public int addNumbers(int n, int m){
    	return matrix1[n][m] + matrix2[m][n];
    }
    
    public void multiplyMatrices(){
    	int[][] resultmatrix = new int[matrix1.length][matrix2[0].length];
    	for(int i = 0; i<matrix1.length; i++){
    		for(int j = 0; j<matrix2[0].length; j++){
    			for(int k = 0; k<matrix1[0].length ; k++){
    				resultmatrix[i][j] += matrix1[i][k]*matrix2[k][j];
    				                              
    			}
    		}
    	}
    }
    
    
}
