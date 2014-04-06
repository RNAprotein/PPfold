package com.ppfold.algo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import com.ppfold.algo.extradata.ExtraDataBars;

import junit.framework.TestCase;

public class MLETest extends TestCase {
	public void setUp() {
	}

	public void testTrees() throws Exception {
		String paramfilename = "matrices.in";
		BufferedReader paramFileReader = null;
        Parameters param; 
		try {
			paramFileReader = new BufferedReader(new InputStreamReader(
					Thread.currentThread().getContextClassLoader().
					getResourceAsStream(paramfilename)));
			param = Parameters.readParam(paramFileReader);
		} catch (Exception e) {
			throw new Exception("Error reading parameter file! " +
					"Check that the file name and format are OK.");
		}
		finally {
			paramFileReader.close();
        }
		double [][] D = param.getrD();
		double [][] V = param.getrV();
		double [][] V1 = param.getrV1();
		double [] Pr = param.getPr();
		
		double[] vectorA = new double[]{1d, 0d, 0d, 0d};
		double [] vectorB = new double[]{0d, 0d, 1d, 0d};
		double [] vectorNA = new double[]{1d, 1d, 1d, 1d};
		double t = 2.3363; //distance between A and B
		double [][] expRt = MatrixTools.expRT(D, t,  V, V1);
		double [][] expR0 = MatrixTools.expRT(D, 0d, V, V1);
		double [] tmp = new double[4];
		double [] tmp2 = new double[4];
		double [] vecAcopy = vectorA.clone();
		double [] vecBcopy = vectorB.clone();
		
		System.out.println("Calculated prob of 'extra rooted' tree:");

		MatrixTools.multiplyMatrixVector(expR0,vecAcopy, tmp);
		MatrixTools.multiplyMatrixVector(expRt,vecBcopy, tmp2);
		System.out.println("Down-top A:");
		MatrixTools.prints(vecAcopy);
		System.out.println("Down-top B: ");
		MatrixTools.prints(vecBcopy);
		MatrixTools.multiplySeries(tmp, tmp2);
		System.out.println("Down-top A <*> Down-top B:");
		MatrixTools.prints(tmp);
		double result1 = MatrixTools.scalarProduct(tmp, Pr);
		System.out.println("Probability: ");
		System.out.println(result1);
		
		System.out.println();
		System.out.println("The same, calculated on node A: ");
		tmp = new double[4];
		tmp2 = new double[4];	
		vecAcopy = vectorA.clone();
		vecBcopy = vectorB.clone();
		MatrixTools.multiplyMatrixVector(expRt,vecBcopy,tmp);
		System.out.println("Downtop B:");
		MatrixTools.prints(vecBcopy);
		MatrixTools.multiplyVectorMatrix(vecBcopy, expR0, tmp2);
		System.out.println("Upbottom A = uptopA expR0:");
		MatrixTools.prints(vecBcopy);
		MatrixTools.multiplySeries(tmp, vectorA);
		System.out.println("(Uptop A expR0) <*> downbottom A: ");		
		MatrixTools.prints(tmp);
		double result1a = MatrixTools.scalarProduct(tmp, Pr);
		System.out.println("Probability: ");
		System.out.println(result1a);
		
		System.out.println();
		System.out.println("Calculated prob of 'leaf in root' tree, in A:");
		tmp = new double[4];
		tmp2 = new double[4];
		vecAcopy = vectorA.clone();
		vecBcopy = vectorB.clone();
		System.out.println("Down-top B:");
		MatrixTools.multiplyMatrixVector(
				expRt, vecBcopy, tmp);
		MatrixTools.prints(vecBcopy);
		System.out.println("expR(0) vectorA:");
		MatrixTools.multiplyMatrixVector(
				expR0, vecAcopy, tmp2);
		System.out.println("Down-bottom A:");
		MatrixTools.multiplySeries(vecBcopy, vecAcopy);
		//MatrixTools.multiplySeries(vecBcopy, vectorA);
		MatrixTools.prints(vecBcopy);
		double result2 = MatrixTools.scalarProduct(vecBcopy, Pr);
		System.out.println("(Uptop A is [1 1 1 1])");
		System.out.println("Probability: ");
		System.out.println(result2);

		System.out.println();
		System.out.println("Calculated prob of 'leaf in root' tree, in B:");
		tmp = new double[4];
		tmp2 = new double[4];
		vecAcopy = vectorA.clone();
		vecBcopy = vectorB.clone();
		double [] vecBcopy2 = vectorB.clone();
		double[] vecNA = new double[]{1d, 1d, 1d, 1d}; 
		System.out.println("Down-bottom B:");
		MatrixTools.prints(vecBcopy);
		System.out.println("Down-top B:");
		MatrixTools.multiplyMatrixVector(
				expRt, vecBcopy, tmp);
		MatrixTools.prints(vecBcopy);
		System.out.println("Uptop A = [1 1 1 1]");
		System.out.println("Upbottom A: ");
		MatrixTools.multiplyMatrixVector(
				expR0, vecNA, tmp2);
		MatrixTools.prints(vecNA);
		System.out.println("Uptop B = upbottom A");
		System.out.println("(Uptop B expRt) <*> downbottom B: ");
		double [] savedved = vecNA.clone();
		MatrixTools.multiplySeries(vecNA, vecBcopy2);
		MatrixTools.prints(vecNA);
		//MatrixTools.multiplySeries(vecBcopy, vectorA);
		double result2a = MatrixTools.scalarProduct(vecNA, Pr);
		System.out.println("Probability:");
		System.out.println(result2a);
		
		
		System.out.println("(Uptop B expRt) <*> downbottom B: ");
		
		
		
	}
}
