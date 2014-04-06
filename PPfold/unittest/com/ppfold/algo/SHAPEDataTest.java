package com.ppfold.algo;

import com.ppfold.algo.extradata.ExtraDataBars;

import junit.framework.TestCase;

public class SHAPEDataTest extends TestCase {
	
	public void setUp() {
	}

	public void testReading() throws Exception {
		String filename = "res/SHAPE_ribosome_all.dat";
		ExtraDataBars readdata = ExtraDataBars.readDistTable(filename);
		filename = "testfile/HIV_SHAPE.txt";
		readdata.importData(filename, 9173);
		
	}
	
}
