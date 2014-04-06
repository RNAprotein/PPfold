package com.ppfold.main;

import junit.framework.TestCase;

public class NewickReaderTest extends TestCase {

	public void setUp() {
	}

	public void testNewick() throws Exception {
	 
		System.out.println("GCA TEST:");
		NewickReader.parse("((gca_rat:0.047,gca_something:0.0767):0.3077,gca_chicken:0.3981,gca_bovine:0.1589);").print();
		
		System.out.println();
		System.out.println("Complicated name and number test:");
		NewickReader.parse("(A:0.1,(B:0.2,H:0.1):0.5,(C:0.3,D:0.4)E:0.5);").print();
		
		try{
			System.out.println();
			System.out.println("Failure handling test:");
			NewickReader.parse("(:0.1,:0.2,:0.3,:0.4):0.5):0.0;");
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		try{
			System.out.println();
			System.out.println("Failure handling test 2:");
			NewickReader.parse("(A:0.1,(B:0.2,H:0.1:0.5,(C:0.3,D:0.4)E:0.5);");
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		try{
			System.out.println();
			System.out.println("Failure handling test 3:");
			NewickReader.parse(":0.1,:0.2,:0.3,:0.4):0.5):0.0);");
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		System.out.println();
		System.out.println("Longer test with root removal:");
		NewickReader.parse("(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0");
		
		System.out.println();
		System.out.println("Another test:");
		NewickReader.parse("((a:1,b:2):3,c:4):5").print();
		
		
		System.out.println();
		System.out.println("Just number test:");
		NewickReader.parse("(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;").print();
		
		System.out.println();
		System.out.println("FINISHED TESTING");

	}
	
	
}
