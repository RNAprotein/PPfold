package com.ppfold.algo;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import junit.framework.TestCase;

public class BenchmarkingSCFGTest extends TestCase {

	public void setUp() {
	}

	
	public void testTiming() throws Exception {
//		Tests with one leaf in tree

		File file  = new File("../PPfold/res/output.log");
		PrintStream printStream = new PrintStream(new FileOutputStream(file));
	//	System.setOut(printStream);

		
		System.out.println("Running timing test: ");

		String my_path = "../PPfold/res/matrices.in";
		Parameters param = Parameters.readParam(new BufferedReader(new FileReader(new File(my_path))));

		List<String> names = new ArrayList<String>();
		List<char[]> columns = new ArrayList<char[]>();

		String sequence = "GGTCTCTCTTGCTAGACCAGATTTGAGCCTGGGAGCTCTCTGACTAGCAGGGGAACCCACT" + 
		"GCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTATAGTAGTGTGTGCCCGTCTGTTGTGT" + 
		"GACTCTGGTAACTAGAGATCCCTCAGACCACTCTAGACAGTGTAAAAATCTCTAGCAGTGG" + 
		"CGCCCGAACAGGGACTTGAAAGCGAAAGTTAATAGGGACTCGAAAGCGAAAGTTCCAGAGA" +
		"AGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTGCACACAGCAAGAGGCGAGAGCGGCGA" +
		"CTGGTGAGTACGCCATTTTTGACTAGCAGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCG" +
		"TCAGTATTAAGTGGGGGAAAATTAGATGAATGGGAAAAAATTCGGTTACGGCCAGGGGGAA" +
		"AGAAAAAATATAAAATGAAACACCTAGTATGGGCAAGCAGGGAGCTGGAAAGATTCGCAAT" + 
		"TAACCCTGGCCTTTTAGAAACAGCAGAAGGATGTCAACAGATAATAGAGCAGTTACAATCA" +
 		"ACTCTCAAGACAGGATCAGAAGAACTTAAATCATTATTTAATACAGTAGCAACCCTCTATT" + 
		"GTGTACATCAAAGGATAAAAGTAACAGACACCAAGGAAGCCATAGATAAAATAGAGGAAAT" + 
		"ACAAAAGAAGAGCAAGCAAAAGGCACATCAGGCAGCTGCCACAGGAAACAGCAGTAATCTC" + 
		"AGTCAAAATTACCCCATAGTGCAAAATGCACAAGGGCAAATGGTACATCAGGCCATGTCAC" + 
		"CACGAACTTTAAATGCATGGGTGAAAGTAATAGAAGATAAGGCTTTCAGCCCAGAAGTAAT" + 
		"ACCCATGTTCACAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACATGATGCTGAAC" + 
		"GTGTACATCAAAGGATAAAAGTAACAGACACCAAGGAAGCCATAGATAAAATAGAGGAAAT" + 
		"ACAAAAGAAGAGCAAGCAAAAGGCACATCAGGCAGCTGCCACAGGAAACAGCAGTAATCTC" + 
		"AGTCAAAATTACCCCATAGTGCAAAATGCACAAGGGCAAATGGTACATCAGGCCATGTCAC" + 
		"CACGAACTTTAAATGCATGGGTGAAAGTAATAGAAGATAAGGCTTTCAGCCCAGAAGTAAT" + 
		"ACCCATGT" + 
		"ATAGTGGGGGGACACCAGGCAGCTATGCAAATG";
		
		System.out.println("Length of sequence: " + sequence.length());

		Node root = new Node();
		root.setName("Root");
		root.setId(0);
		Tree tree = new Tree(root);

		Node node1 = new Node();
		node1.setId(1);
		node1.setName("sequence");
		root.addChild(node1);

		node1.setDistanceFromParent(0.5);

		names.add("sequence");

		char[] thiscolumn = new char[1]; 


		for(int i = 0; i < sequence.length(); i++){
			//iterates column number
			thiscolumn[0]=(char) sequence.charAt(i);
			columns.add(thiscolumn.clone());
		}


		//AsynchronousJobExecutor executor = new AsynchronousJobExecutorBlocking();
		 AsynchronousJobExecutor executor;
		
		//a random shuffling of distances 
		int[] distances = new int[35];
	
		int maxpool = 2; 
		
		for(int i = maxpool; i>0; i--){
			executor = new AsynchronousJobExecutorThreadPool(i);
			distances  = shuffle(distances); 
			System.out.println("*****************");
			System.out.println("POOL("+i+") EXECUTOR");
			System.out.println("*****************");
			
			System.out.println("SET ONE");
			for(int d : distances){
				System.out.println("------------------------------");
				System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
				System.out.println();
				FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
			
			}

			System.out.println("SET TWO");
			distances  = shuffle(distances); 
			
			for(int d : distances){
				System.out.println("------------------------------");
				System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
				FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
			}
			
			System.out.println("SET THREE");
			distances  = shuffle(distances); 
		
			for(int d : distances){
				System.out.println("------------------------------");
				System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
				FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
			}
			
			System.out.println("SET FOUR");
			distances  = shuffle(distances); 
		
			for(int d : distances){
				System.out.println("------------------------------");
				System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
				FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
			}
			
			System.out.println("SET FIVE");
			distances  = shuffle(distances); 
		
			for(int d : distances){
				System.out.println("------------------------------");
				System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
				FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);

			}
			
			System.out.println("SET LARGE");
			int[] largedivisions = {40, 50, 60, 70, 80, 90, 100, 250, 500, 100};  
			for(int d : largedivisions){
				System.out.println("------------------------------");
				System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
				FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
			}
			
		}
		
		
		System.out.println("*****************");
		System.out.println("BLOCKING EXECUTOR");
		System.out.println("*****************");
		
		executor = new AsynchronousJobExecutorBlocking();
		
		distances  = shuffle(distances); 
		
		System.out.println("SET ONE - BLOCKING");
		for(int d : distances){
			System.out.println("------------------------------");
			System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
			System.out.println();
			FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		
		}

		System.out.println("SET TWO  - BLOCKING");
		distances  = shuffle(distances); 
		
		for(int d : distances){
			System.out.println("------------------------------");
			System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
			FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		}
		
		System.out.println("SET THREE  - BLOCKING");
		distances  = shuffle(distances); 
	
		for(int d : distances){
			System.out.println("------------------------------");
			System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
			FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		}
		
		System.out.println("SET FOUR  - BLOCKING");
		distances  = shuffle(distances); 
	
		for(int d : distances){
			System.out.println("------------------------------");
			System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
			FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		}
		
		System.out.println("SET FIVE  - BLOCKING");
		distances  = shuffle(distances); 
	
		for(int d : distances){
			System.out.println("------------------------------");
			System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
			FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		}
		
		System.out.println("SET LARGE  - BLOCKING");
		int[] largedivisions = {40, 50, 60, 70, 80, 90, 100, 250, 500, 100};  
		for(int d : largedivisions){
			System.out.println("------------------------------");
			System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
			FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		}
		
		System.out.println("FINISHED SUCCESSFULLY");
		assertTrue(true);

	}
	
	
	private static int[] shuffle(int[] distances){
		Random rgen = new Random();  
		
		for (int i=0; i<distances.length; i++) {
		    distances[i] = i+1;
		}

		//--- Shuffle by exchanging each element randomly
		for (int i=0; i<distances.length; i++) {
		    int randomPosition = rgen.nextInt(distances.length);
		    int temp = distances[i];
		    distances[i] = distances[randomPosition];
		    distances[randomPosition] = temp;
		}
		
		return distances;
	}
	
}
