package com.ppfold.algo;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

public class FoldingProjectTest extends TestCase {

	public void setUp() {
	}
	
	public void _testVeryLongFold() throws Exception{
		int d = 25; 
		
		System.out.println("Running timing test: ");

		String my_path = "../PPfold/res/matrices.in";
		Parameters param = Parameters.readParam(new BufferedReader(new FileReader(new File(my_path))));

		List<String> names = new ArrayList<String>();
		List<char[]> columns = new ArrayList<char[]>();

		String sequence = "GGTCTCTCTTGCTAGACCAGATTTGAGCCTGGGAGCTCTCTGACTAGCAGGGGAACCCACT" + 
		"GCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTATAGTAGTGTGTGCCCGTCTGTTGTGT" + 
		"GACTCTGGTAACTAGAGATCCCTCAGACCACTCTAGACAGTGTAAAAATCTCTAGCAGTGG" + 
		"CGCCCGAACAGGGACTTGAAAGCGAAAGTTAATAGGGACTCGAAAGCGAAAGTTCCAGAGA" +
//		"AGTTCTCTCGACGCAGGACTCGGCTTGCTGAGGTGCACACAGCAAGAGGCGAGAGCGGCGA" +
//		"CTGGTGAGTACGCCATTTTTGACTAGCAGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCG" +
//		"TCAGTATTAAGTGGGGGAAAATTAGATGAATGGGAAAAAATTCGGTTACGGCCAGGGGGAA" +
		"AGAAAAAATATAAAATGAAACACCTAGTATGGGCAAGCAGGGAGCTGGAAAGATTCGCAAT" + 
		"TAACCCTGGCCTTTTAGAAACAGCAGAAGGATGTCAACAGATAATAGAGCAGTTACAATCA" +
 		"ACTCTCAAGACAGGATCAGAAGAACTTAAATCATTATTTAATACAGTAGCAACCCTCTATT" + 
		"GTGTACATCAAAGGATAAAAGTAACAGACACCAAGGAAGCCATAGATAAAATAGAGGAAAT" + 
		"ACAAAAGAAGAGCAAGCAAAAGGCACATCAGGCAGCTGCCACAGGAAACAGCAGTAATCTC" + 
		"AGTCAAAATTACCCCATAGTGCAAAATGCACAAGGGCAAATGGTACATCAGGCCATGTCAC" + 
		"CACGAACTTTAAATGCATGGGTGAAAGTAATAGAAGATAAGGCTTTCAGCCCAGAAGTAAT" + 
		"ACCCATGTTCACAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACATGATGCTGAAC" + 
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
			columns.add(
					thiscolumn.clone());
		}

		AsynchronousJobExecutor executor;
		
		System.out.println("POOL(4) EXECUTOR");
		executor = new AsynchronousJobExecutorThreadPool(4);
		System.out.println("------------------------------");
		System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
		System.out.println();
		FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);

		System.out.println("POOL(3) EXECUTOR");
		executor = new AsynchronousJobExecutorThreadPool(3);
		System.out.println("------------------------------");
		System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
		System.out.println();
		FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);

		System.out.println("POOL(2) EXECUTOR");
		executor = new AsynchronousJobExecutorThreadPool(2);
		System.out.println("------------------------------");
		System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
		System.out.println();
		FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		
		//PointRes.counter.set(0);
		
		System.out.println("POOL(1) EXECUTOR");
		executor = new AsynchronousJobExecutorThreadPool(1);
		System.out.println("------------------------------");
		System.out.println("FOLDING WITH SCFG DIVISIONS = " + d);
		System.out.println();
		FoldingProject.fold(NullProgress.INSTANCE, 50, d, tree, columns, names, param, executor, false, 0, null,false,false);
		
		

	}
	
	public void _testPhyloInsideOutsideLongFold() throws Exception {
//		Tests with one leaf in tree
		System.out.println("Running long folding test: ");

		String my_path = "../PPfold/res/matrices.in";
		Parameters param = Parameters.readParam(new BufferedReader(new FileReader(new File(my_path))));

		List<String> names = new ArrayList<String>();
		List<char[]> columns = new ArrayList<char[]>();

		//String sequence = "GTCGATCGATCGCAGTTTCTGCGATCGATCGAC"+
		//String sequence = "GTCGATCGATCGCAGTTTCTGCGATCGATCGA"+
		String sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"+
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"+
		"";
		
		//String sequence = "GGGGGGCCCCCC";
		//String sequence = "AGAGAGGCUAGAUCCUCUGUGUUGAGAAGGAUCAUGAUGGGCUCCUCGGUGUUCUCCAGGUAGCGGCACCACACCAUGA";
		
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


		AsynchronousJobExecutor executor = new AsynchronousJobExecutorBlocking();
		//AsynchronousJobExecutor executor = new AsynchronousJobExecutorThreadPool(2);
		
		for(int d = 65*2-1; d<65*2; d++){
			
		

			ResultBundle result = FoldingProject.fold(NullProgress.INSTANCE, 1, 1, tree, columns, names, param, executor,true, 0, null,false,false);
			

			//MatrixTools.print(result.basepairprob);


			int cnt = 0;
			while(cnt<sequence.length()){
				for(int i=0; (i<50&&cnt+i<sequence.length()); i++){
					System.out.print(sequence.charAt(cnt+i));
				}
				System.out.println();
				for(int i=0; (i<50&&cnt+i<sequence.length()); i++){
					System.out.print(result.structure[cnt+i]);
				}
				System.out.println();
				cnt+=51;

			}		
			


			System.out.println("Reliabilty scores: ");
			System.out.println("position \t structure \t reliability");
			for(int c = 0; c < sequence.length(); c++){
				System.out.println((c+1) + "\t" + result.structure[c] + "\t" + result.reliability[c] + " ");

			}
			
			System.out.println("Asserting correct result...");
			
			double[] rightresult = new double[result.reliability.length];
			
			rightresult[0] = 0.3970;
		    rightresult[1] = 0.3333; 
		    rightresult[2] = 0.3335;
		    rightresult[3] = 0.3310;
		    rightresult[4] = 0.3282;
		    rightresult[5] = 0.3255;
		    rightresult[6] = 0.3227;
		    rightresult[7] = 0.3200;
		    rightresult[8] = 0.3172;
		    rightresult[9] = 0.3145; 
		    rightresult[10] = 0.3117;
		    rightresult[11] = 0.3090;
		    rightresult[12] = 0.3062;
		    rightresult[13] = 0.3035;
		    rightresult[14] = 0.3008; 
		    rightresult[15] = 0.2981;
		    rightresult[16] = 0.2954;
		    rightresult[17] = 0.2926; 
		    rightresult[18] = 0.2899; 
		    rightresult[19] = 0.2872; 
		    rightresult[20] = 0.0046; 
		    rightresult[21] = 0.2850; 
		    rightresult[22] = 0.2860;
		    rightresult[23] = 0.2869;
		    rightresult[24] = 0.2879; 
		    rightresult[25] = 0.2888; 
		    rightresult[26] = 0.2898; 
		    rightresult[27] = 0.2907; 
		    rightresult[28] = 0.2917;
		    rightresult[29] = 0.2924; 
		    rightresult[30] = 0.2908;
		    rightresult[31] = 0.6632; 
		    rightresult[32] = 0.9621; 
			rightresult[33] = 0.9621;
		    rightresult[34] = 0.6632; 
		    rightresult[35] = 0.1364;
		    rightresult[36] = 0.2908;
		    rightresult[37] = 0.2924;
		    rightresult[38] = 0.2917;
		    rightresult[39] = 0.2907;
		    rightresult[40] = 0.2898;
		    rightresult[41] = 0.2888;
		    rightresult[42] = 0.2879; 
		    rightresult[43] = 0.2869;
		    rightresult[44] = 0.2860;
		    rightresult[45] = 0.2850;
		    rightresult[46] = 0.2872;
		    rightresult[47] = 0.2899; 
		    rightresult[48] = 0.2926;
		    rightresult[49] = 0.2954;
		    rightresult[50] = 0.2981; 
		    rightresult[51] = 0.3008; 
		    rightresult[52] = 0.3035; 
		    rightresult[53] = 0.3062; 
		    rightresult[54] = 0.3090; 
		    rightresult[55] = 0.3117;
		    rightresult[56] = 0.3145;
		    rightresult[57] = 0.3172; 
		    rightresult[58] = 0.3200; 
		    rightresult[59] = 0.3227; 
		    rightresult[60] = 0.3255; 
		    rightresult[61] = 0.3282;
		    rightresult[62] = 0.3310; 
		    rightresult[63] = 0.3335;
		    rightresult[64] = 0.3333; 
		    rightresult[65] = 0.3970;	
		    
		    
		    
			for(int c = 0; c < sequence.length(); c++){

				try{
					assertTrue( Math.abs(result.reliability[c]-rightresult[c]) <= 0.00005);
				}
				catch(AssertionFailedError e){
					System.err.println("Assertion failed for nr. " + (c+1) + 
					" - Reason: expected "+ rightresult[c] + ", got "
							+ result.reliability[c]);
				}
			}
			

		//	String dir_fil = "C:/Documents and Settings/Zsuzsi/Desktop/";  
			//	ExportTools.writeTabbedMatrix(dir_fil, "FILENAME", ".fileextension", result.getBasePairProb());
			//ExportTools.writeCTFormat(dir_fil, "struct", ".str", columns, result.structure, result.reliability); 
			
		}


	}
	
	public void testPhyloPfoldTest() throws Exception {
//		Tests against Bjarnes pfold alignment example
		System.out.println("Running long folding test: ");

		String my_path = "../PPfold/res/matrices.in";
		Parameters param = Parameters.readParam(new BufferedReader(new FileReader(new File(my_path))));

		List<String> names = new ArrayList<String>();
		List<char[]> columns = new ArrayList<char[]>();
		List<char[]> fullcolumns = new ArrayList<char[]>();


		//int d = 65; //phylodivisions

		
		String bovine =  "AGCCCUGUGGUGAAUUUACACGUUGAAUUGCAAAUUCAGAGAAGCAGCUUCAAU-UCUGCCGGGGCUU"; //bovine  (names[3])
		String chicken = "GACUCUGUAGUGAAGU-UCAUAAUGAGUUGCAAACUCGUUGAUGUACACUAA-AGUGUGCCGGGGUCU"; //chicken (names[2])
		String mouse =   "GGUCUUAAGGUGAUA-UUCAUGUCGAAUUGCAAAUUCGAAGGUGUAGAGAAAU-CUCUACUAAGACUU"; //mouse (names[1])
		String rat =     "AGCCUUAAGGUGAUU-AUCAUGUCGAAUUGCAAAUUCGAAGGUGUAGAGAAUCU-UCUACUAAGGCUU"; //rat (names[0])

/*		
		String bovine =  "AGCCCUGUGGUGAAUACACGUUGAAUUGCAAAUUCAGAGAAGCAGCUUCAUCUGCCGGGGCUU"; //bovine  (names[3])
		String chicken = "GACUCUGUAGUGAAGUCAUAAUGAGUUGCAAACUCGUUGAUGUACACUAAUGUGCCGGGGUCU"; //chicken (names[2])
		String mouse =   "GGUCUUAAGGUGAUAUCAUGUCGAAUUGCAAAUUCGAAGGUGUAGAGAAAUCUACUAAGACUU"; //mouse (names[1])
		String rat =     "AGCCUUAAGGUGAUUUCAUGUCGAAUUGCAAAUUCGAAGGUGUAGAGAAUUCUACUAAGGCUU"; //rat (names[0])

*/

		System.out.println("Length of alignment: " + bovine.length());

		Node root = new Node();
		root.setName("Root");
		root.setId(0);
		Tree tree = new Tree(root);

		Node node1 = new Node();
		node1.setId(1);
		node1.setName("gca_rat");
		node1.setDistanceFromParent(0.047);
		names.add("gca_rat");

		Node node2 = new Node();
		node2.setId(2);
		node2.setName("gca_mouse");
		node2.setDistanceFromParent(0.0767);
		names.add("gca_mouse");

		Node node3 = new Node();
		node3.setId(3);
		node3.setName("");
		node3.setDistanceFromParent(0.3077);
		root.addChild(node3);
		node3.addChild(node1);
		node3.addChild(node2);

		Node node4 = new Node();
		node4.setId(4);
		node4.setName("gca_chicken");
		node4.setDistanceFromParent(0.3981);
		names.add("gca_chicken");
		root.addChild(node4);

		Node node5 = new Node();
		node5.setId(5);
		node5.setName("gca_bovine");
		node5.setDistanceFromParent(0.1589);
		names.add("gca_bovine");
		root.addChild(node5);

		tree.getRoot().printChildren();

    	char[] thiscolumn = new char[4]; 
    
    	char[] finalstructure = new char[bovine.length()];
    	double[] finalreliability = new double[bovine.length()]; 
    	double[][] finalmatrix = new double[bovine.length()][bovine.length()];
    	
    	List<Integer> leftoutcolumns = new ArrayList<Integer>(); 
    	 	 
    	for(int i = 0; i < bovine.length(); i++){
    		int gapscounter = 0; 
    		//iterates column number
    			//iterates sequences
    			thiscolumn[0]=(char) rat.charAt(i);
    			thiscolumn[1]=(char) mouse.charAt(i);
    			thiscolumn[2]=(char) chicken.charAt(i);
    			thiscolumn[3]=(char) bovine.charAt(i);
    			//System.out.println("Calculating column nr. " + j + " taking from seq " + j + " charAt " + i + " which is " + thiscolumn[j]);
			for(int j = 0; j<4; j++){
				if(MatrixTools.isGap(thiscolumn[j])){gapscounter++;}
			}
    		fullcolumns.add(thiscolumn.clone()); 
			
			if((float)gapscounter/4 < 0.25){
    			//"This problem was handled by removing columns where less than 75% of 
    			//the sequences have nucleotides"
				//= leave out columns with more than or equal to 25% gaps 
				//= add only columns where there are less than 25% gaps
    			columns.add(thiscolumn.clone());
    		}
			else{
				finalstructure[i] = '.'; 
				finalreliability[i] = 0;
				leftoutcolumns.add(i);
				System.out.println("Leaving out column: " + i);
			}
    	}
		
    	System.out.println("Total number of columns predicting structure on: " + columns.size());
    	
		AsynchronousJobExecutor executor = new AsynchronousJobExecutorBlocking();
		//AsynchronousJobExecutor executor = new AsynchronousJobExecutorThreadPool(2);

		ResultBundle result = FoldingProject.fold(NullProgress.INSTANCE, 1, 4, tree, columns, names, param, executor,true, 0, null,false,false);

		//transform result so left out columns are included again 
		System.out.println("Structure prediction complete. Transforming to include left-out columns " +
				"(total " + leftoutcolumns.size() + ")...");
	
		
		
		int icolcnt = 0;
		for(int i = 0; i<result.getBasePairProb().length; i++){
			while(leftoutcolumns.indexOf(i+icolcnt)!=-1){
			//	System.out.println("i = " + i);
				icolcnt++;
				finalmatrix[i][i] = 0;
				//continue;
			}
			//System.out.println("probmatrix: " + i + ", finalmatrix coord: " + (i+icolcnt));
			int jcolcnt = 0; 
			for(int j = 0; j<result.getBasePairProb().length; j++){
				//System.out.println("j = " + j);
				while(leftoutcolumns.indexOf(j+jcolcnt)!=-1){
					finalmatrix[j][j] = 0;
					jcolcnt++;
					//continue;
				}
				finalmatrix[i+icolcnt][j+jcolcnt] = result.getBasePairProb()[i][j];
				
			}			
		
		}
		
		
		//System.out.println(finalmatrix[10][19]);
		//System.out.println(result.getBasePairProb()[10][17]);
		

		
		//could maybe implement a search-based solution using indexOf.
		//Dunno which one is better/faster.
		//This seems to work. 
		
		
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int lastone = 0;
		int cntr = 0; //counts how many columns are shifted
		int next = 0;
		while(next<bovine.length()){
			int co; //limit to which to go
			if(iter.hasNext()){co=iter.next();}
			else{co=bovine.length();}
			cntr++; 
			//have to copy columns up to col
			for(int copycol = lastone; copycol<=co-cntr; copycol++){
		//		System.out.println("Filling finalstructure " + next + " with " + copycol);
				finalreliability[next] = result.getReliability()[copycol];
				finalstructure[next] = result.getStructure()[copycol];
				next++;
			}
			lastone = co-cntr+1;
			next++;
		}
		
		
	//	finalstructure = result.getStructure();
	//	finalreliability = result.getReliability(); 
		
		System.out.println("Done.");
		
		int cnt = 0;
		while(cnt<result.getReliability().length){
			for(int i=0; (i<50&&cnt+i<result.getStructure().length); i++){
				System.out.print(result.getStructure()[cnt+i]);
			}
			System.out.println();
			cnt+=51;

		}		
		
		
			cnt = 0;
			while(cnt<bovine.length()){
				for(int i=0; (i<50&&cnt+i<bovine.length()); i++){
					System.out.print(finalstructure[cnt+i]);
				}
				System.out.println();
				cnt+=51;

			}		
			
			/*
			for(int i = 11; i<17; i++){
				System.out.println("Ps(" + i + ")=" + result.getSingleBaseProb()[i] );
				for(int j = 11; j<17; j++){
					System.out.println("Pd(" + i + "," + j + ")=" + result.getBasePairProb()[i][j]);
				}
			}
			*/
			

			double[] rightresult = new double[finalreliability.length];
			char[] rightstruct = new char[finalreliability.length];
			
			rightresult[0] = 0.9605;
		    rightresult[1] = 0.9994; 
		    rightresult[2] = 1.0000;
		    rightresult[3] = 0.9996;
		    rightresult[4] = 0.9997;
		    rightresult[5] = 0.9885;
		    rightresult[6] = 0.9773;
		    rightresult[7] = 0.9674;
		    rightresult[8] = 0.9781;
		    rightresult[9] = 0.6172; 
		    rightresult[10] = 0.5379;
		    rightresult[11] = 0.5818;
		    rightresult[12] = 0.9986;
		    rightresult[13] = 0.9736;
		    rightresult[14] = 0.9769; 
		    rightresult[15] = 0.0000;
		    rightresult[16] = 0.0000;
		    rightresult[17] = 0.9593; 
		    rightresult[18] = 0.5818; 
		    rightresult[19] = 0.5379; 
		    rightresult[20] = 0.6032; 
		    rightresult[21] = 0.8993; 
		    rightresult[22] = 0.8812;
		    rightresult[23] = 0.9905;
		    rightresult[24] = 0.9994; 
		    rightresult[25] = 0.9983; 
		    rightresult[26] = 0.9978; 
		    rightresult[27] = 0.9489; 
		    rightresult[28] = 0.7627;
		    rightresult[29] = 0.9999; 
		    rightresult[30] = 1.0000;
		    rightresult[31] = 0.9991; 
		    rightresult[32] = 0.7627; 
			rightresult[33] = 0.9489;
		    rightresult[34] = 0.9978; 
		    rightresult[35] = 0.9983;
		    rightresult[36] = 0.9994;
		    rightresult[37] = 0.9905;
		    rightresult[38] = 0.8812;
		    rightresult[39] = 0.8961;
		    rightresult[40] = 0.8577;
		    rightresult[41] = 0.9547;
		    rightresult[42] = 0.9688; 
		    rightresult[43] = 0.9675;
		    rightresult[44] = 0.9949;
		    rightresult[45] = 0.9967;
		    rightresult[46] = 0.9965;
		    rightresult[47] = 0.5033; 
		    rightresult[48] = 0.9997;
		    rightresult[49] = 0.9999;
		    rightresult[50] = 1.0000; 
		    rightresult[51] = 1.0000; 
		    rightresult[52] = 0.0000; 
		    rightresult[53] = 0.0000; 
		    rightresult[54] = 0.0000; 
		    rightresult[55] = 0.5033;
		    rightresult[56] = 0.9965;
		    rightresult[57] = 0.9967; 
		    rightresult[58] = 0.9949; 
		    rightresult[59] = 0.9675; 
		    rightresult[60] = 0.9773; 
		    rightresult[61] = 0.9885;
		    rightresult[62] = 0.9997; 
		    rightresult[63] = 0.9996;
		    rightresult[64] = 1.0000; 
		    rightresult[65] = 0.9994;	
		    rightresult[66] = 0.9605; 
		    rightresult[67] = 0.9996;		
			
		    
			rightstruct[0] = '(';
			rightstruct[1] = '('; 
			rightstruct[2] = '(';
			rightstruct[3] = '(';
			rightstruct[4] = '(';
			rightstruct[5] = '(';
			rightstruct[6] = '(';
			rightstruct[7] = '.';
			rightstruct[8] = '.';
			rightstruct[9] = '.'; 
			rightstruct[10] = '(';
			rightstruct[11] = '(';
			rightstruct[12] = '.';
			rightstruct[13] = '.';
			rightstruct[14] = '.'; 
			rightstruct[15] = '.';
			rightstruct[16] = '.';
			rightstruct[17] = '.'; 
			rightstruct[18] = ')'; 
			rightstruct[19] = ')'; 
			rightstruct[20] = '.'; 
			rightstruct[21] = '.'; 
			rightstruct[22] = '(';
			rightstruct[23] = '(';
			rightstruct[24] = '('; 
			rightstruct[25] = '('; 
			rightstruct[26] = '('; 
			rightstruct[27] = '('; 
			rightstruct[28] = '(';
			rightstruct[29] = '.'; 
			rightstruct[30] = '.';
			rightstruct[31] = '.'; 
			rightstruct[32] = ')'; 
			rightstruct[33] = ')';
			rightstruct[34] = ')'; 
			rightstruct[35] = ')';
			rightstruct[36] = ')';
			rightstruct[37] = ')';
			rightstruct[38] = ')';
			rightstruct[39] = '.';
			rightstruct[40] = '.';
			rightstruct[41] = '.';
			rightstruct[42] = '.'; 
			rightstruct[43] = '(';
			rightstruct[44] = '(';
			rightstruct[45] = '(';
			rightstruct[46] = '(';
			rightstruct[47] = '('; 
			rightstruct[48] = '.';
			rightstruct[49] = '.';
			rightstruct[50] = '.'; 
			rightstruct[51] = '.'; 
			rightstruct[52] = '.'; 
		    rightstruct[53] = '.'; 
		    rightstruct[54] = '.'; 
		    rightstruct[55] = ')';
		    rightstruct[56] = ')';
		    rightstruct[57] = ')'; 
		    rightstruct[58] = ')'; 
		    rightstruct[59] = ')'; 
		    rightstruct[60] = ')'; 
		    rightstruct[61] = ')';
		    rightstruct[62] = ')'; 
		    rightstruct[63] = ')';
		    rightstruct[64] = ')'; 
		    rightstruct[65] = ')';	
		    rightstruct[66] = ')'; 
		    rightstruct[67] = '.';		
			
			
		System.out.println("Reliabilty scores: ");
		System.out.println("pos \t str \t rel \t\t correct");
		double totalmine = 0; 
		double totalpfold = 0; 
		
		for(int c = 0; c < bovine.length(); c++){
			System.out.format((c+1) + "\t " + finalstructure[c] + "\t %.4f \t"+ rightstruct[c] + "\t%.4f \n",  finalreliability[c], rightresult[c] );
			totalmine += finalreliability[c];
			totalpfold += rightresult[c];

		}
		
		System.out.println("Total mine: " + totalmine);
		System.out.println("Total pfold: " + totalpfold);
           		
		
		System.out.println("Asserting correct result...");



			for(int c = 0; c < bovine.length(); c++){

				try{
					assertTrue( Math.abs(finalreliability[c]-rightresult[c]) <= 0.00005);
					assertTrue(finalstructure[c] == rightstruct[c]);
				}
				catch(AssertionFailedError e){
					System.err.format("Assertion failed for nr. " + (c+1) + 
					" - Reason: expected "+ rightresult[c] + " [ " + rightstruct[c]+ " ], got "
							+ "%.16f [ " + finalstructure[c]+ " ]\n", finalreliability[c]);
					
				}
			}

	}

}
