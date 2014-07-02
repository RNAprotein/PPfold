package com.ppfold.main;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.ppfold.algo.AsynchronousJobExecutor;
import com.ppfold.algo.AsynchronousJobExecutorThreadPool;
import com.ppfold.algo.ExportTools;
import com.ppfold.algo.FoldingProject;
import com.ppfold.algo.MatrixTools;
import com.ppfold.algo.NeighbourJoining;
import com.ppfold.algo.NullProgress;
import com.ppfold.algo.Parameters;
import com.ppfold.algo.Progress;
import com.ppfold.algo.ResultBundle;
import com.ppfold.algo.Tree;
import com.ppfold.algo.extradata.BinaryStringData;
import com.ppfold.algo.extradata.BinaryStringDataDiffBp;
import com.ppfold.algo.extradata.ExtraData;
import com.ppfold.algo.extradata.ExtraDataBars;
import com.ppfold.algo.extradata.ExtraDataProbMapping;
import com.ppfold.algo.extradata.ForcedConstraints;
import com.ppfold.algo.extradata.SHAPEData;
import com.ppfold.algo.extradata.SHAPEDataFctDiffBp;


public class PPfoldMain implements Runnable {	
	//Global variables are set by GUI or command-line arguments.  
	static String versionnumber = "3.0";
	static String treefilename;
	static String alignmentfilename; 
	static String exportfilehandle;	
	static String seqexportname = null; //if not null, export stuff specific to this sequence.
	static String paramfilename;
	static String paramresname = "matrices.in";
	static String outputdir;
	static String defaultDataDistfile = "dist.dat"; //Default data file

	static boolean auxdata = false;
	static List<DataInfo> datainfo = new ArrayList<DataInfo>();	
	
	static boolean createtree = true;
	static boolean optimizetree = true;
	static boolean verbose = false;	//syang: true;
	static boolean exportson = false; 
	static boolean onlyCT = false;
	static boolean specialname = false;
	static int scfgdivisions = Runtime.getRuntime().availableProcessors()*8;
	static int phylodivisions = Runtime.getRuntime().availableProcessors()*2;
	static int nrprocessors = Runtime.getRuntime().availableProcessors();
	static int iterlimit = 10; 
	static boolean entropycalc = true;
	static String entropyString = null;
	static boolean gui = false;
	
	static boolean userfinished = false; //this will be triggered when the user presses the START button
	static private Progress progress; //Progressbar; either NullActivity (if no GUI), or the PPfoldProgressBar (if GUI) 
	static public boolean shouldstop = false; //this will be triggered when user presses STOP button 
	static public boolean isstopping = false;
	
	//Reset for each project 
	public boolean foldingfinished = false; //thread will stop when folding is finished.
	public boolean success = true;
	public String errormessage  = ""; //if it fails, this will be the output to the user

	public void cleanUp(){
		//Put things here to save memory 
		errormessage = null;
	}
	
	public void run() {
        try{
        	System.out.println("Attempting to read and parse files...");
        	    
        	System.out.println(alignmentfilename);
        	
        	File alignfile = new File (alignmentfilename);
        	if(!specialname){
        		exportfilehandle = alignfile.getName();
        		int lastdotpos = exportfilehandle.lastIndexOf('.');
        		exportfilehandle = exportfilehandle.substring(0,lastdotpos);
        	}
        	Alignment align = AlignmentReader.readAlignment(alignmentfilename);
        	List<ExtraData> extradata = new ArrayList<ExtraData>();

        	BufferedReader paramFileReader = null;
	        Parameters param; 
			try {

				if (paramfilename != null) {
					File file = new File(paramfilename);
					paramFileReader = new BufferedReader(new InputStreamReader(
							new FileInputStream(file)));
				}else{
					paramFileReader = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(paramresname)));
				}
				
	        	param = Parameters.readParam(paramFileReader);
	        } catch (Exception e) {
	            throw new Exception("Error reading parameter file! Check that the file name and format are OK.");
	        }
	        finally {
	            if (paramFileReader != null) {
	                try {
	                    paramFileReader.close();
	                }
	                catch (IOException e) {
	                	throw new Exception("Error closing parameter file.");
	                }
	            }
	        }
			
        	Tree tree;
        	if(!createtree){
        		System.out.println(treefilename);
        	    tree = NewickReader.readNewick(treefilename);
        	}
        	else{
        		tree=null;
        	}
	            		

    		if(auxdata){
				for(int dataID = 0; dataID<datainfo.size(); dataID++){
					ExtraData sequenceData = null;
					DataInfo data = datainfo.get(dataID);
					switch(data.getType()){
					case 0:
						BufferedInputStream shapeDistReader;
						if(data.getDistFileName().equals(defaultDataDistfile)){
							shapeDistReader = new BufferedInputStream(Thread.currentThread().getContextClassLoader().getResourceAsStream(defaultDataDistfile));
						}
						else{
							shapeDistReader = new BufferedInputStream(new FileInputStream(data.getDistFileName()));
						}
		    			sequenceData = ExtraDataBars.readDistTable_toStream(shapeDistReader);
		    			break;
					case 1:
						sequenceData = new ExtraDataProbMapping();
						break;
					case 2:
						sequenceData = new ForcedConstraints();
						break;
					default:
						throw new IllegalArgumentException("Data type could not be interpreted!");
					}

	    			//Find the sequence for which the SHAPE values are known (but only if the name string is not null)
	    			int seqID = -1;
	    			if(data.getSequenceName() != null){
	    				for(int i = 0; i<align.getNames().size(); i++){
	    					if(data.getSequenceName().trim().equals(align.getNames().get(i).trim())){
	    						seqID = i;
	    						data.setSequenceID(seqID);
	    						break;
	    					}
	    				}

	    				if(seqID==-1){
	    					throw new IllegalArgumentException("Sequence with name " + data.getSequenceName() + " could not be found!");
	    				}
	    			}
	    			if(data.getFileName()!=null){ //if null then only contact distance 
	    				BufferedInputStream shapeFileReader = new BufferedInputStream(new FileInputStream(data.getFileName()));
	    				sequenceData.readData_toStream(shapeFileReader, align.calculateLength(seqID));
	    				sequenceData.transformToAlignment(align.getSequences().get(seqID));
	    			}

	    			if(sequenceData instanceof ForcedConstraints && data.getContactDistance()!=-1){
	    					((ForcedConstraints) sequenceData).setContactDistance(data.getContactDistance());
	    			}
	    			
	    			extradata.add(sequenceData);
				
    		}
    		}
    		
        	System.out.println("All files read.");
        	if(outputdir == null){
        		File dir1 = new File (".");
        		try {
        			outputdir = dir1.getCanonicalPath();
        		} catch (IOException e1) {
        			System.out.println("Couldn't open directory ");
        			throw new Exception(e1);
        		}
        	}
        	
        	
			System.out.println("Results will be written to " + outputdir);
        	System.out.println("Starting algorithm...");        	
        	
        	if(seqexportname!=null){
        		System.out.println("Special sequence exported: " + seqexportname + " with prefix " + exportfilehandle + "-" + seqexportname);
        	}
        	
        	//Generate progress if GUI does not take care of it 
    		if(!gui){
    			progress = NullProgress.INSTANCE;
    		}
    		progress.setCurrentActivity("Preparing folding...");
    		fold(progress, align.getSequences(), align.getNames(), tree, param, extradata);		//align.* get the sequence and header name of a fasta alignment
    
            System.out.println();
            System.out.println("FINISHED ");
            success = true;
            foldingfinished = true;        
            errormessage = "";
        }
        catch(InterruptedException e){
        	System.out.println("Folding stopped.");
        	isstopping=true;
        }
        catch(OutOfMemoryError e){
        	System.out.println("Oops! The Java Virtual Machine ran out of memory!");
        	System.out.println("Please increase the heap size.");
        	success = false;
        	String text = new String("Oops! The Java Virtual Machine ran out of memory! (Currently allocated: "+ 
        			(Runtime.getRuntime().totalMemory() / 1048576) + " MB.) "+
        			"You can increase the virtual machine heap size by adding the argument -Xmx[size]. " +
        			"E.g. to run PPfold with 256 MB heap size: \n\n" +
        			"java -Xmx256m -jar PPfold.jar \n\nFor debugging purposes, the stack trace is below. \n");
        	        	
        	text.concat(e.toString()+":\n");
        	for(StackTraceElement s:e.getStackTrace()){
        		text = text.concat(s.toString()+"\n");
        	}
        	errormessage = text;
        	foldingfinished=true;
        	if(!gui){
				e.printStackTrace();
        	}
        }
        catch(Exception e){
        	System.err.println("Fatal error");
        	success = false;
        	String text = new String("");
        	for(StackTraceElement s:e.getStackTrace()){
        		text = text.concat(s.toString()+"\n");
        	}
        	errormessage += e.toString()+":\n"+text;
        	foldingfinished=true;
        	//if(!gui){
			//	e.printStackTrace();
        	//}
        }
    
	}
	
	
	private static void fold(Progress progress, List<String> sequences, List<String> names, Tree tree, 
			Parameters param, List<ExtraData> extradata) throws  InterruptedException, Exception{

		AsynchronousJobExecutor executor = new AsynchronousJobExecutorThreadPool(nrprocessors); 
		progress.setProgress(0.0);

		
		List<char[]> columns = new ArrayList<char[]>();
		//each element represents a *column* of the alignment

		List<char[]> fullcolumns = new ArrayList<char[]>();
		//this will contain the FULL alignment (for outputting later) 

		int nrseq = sequences.size();			//syang: # of species, nr stands for number

		char[] thiscolumn = new char[nrseq]; 
		char[] finalstructure = new char[sequences.get(0).length()];
		float[] finalreliability = new float[sequences.get(0).length()]; 
		float[][] finalmatrix = new float[sequences.get(0).length()][sequences.get(0).length()];
		float[][] finalexp = new float[sequences.get(0).length()][sequences.get(0).length()];
		float[][] finalsingleProb = new float[sequences.get(0).length()][1];	//syang: get single prob
		List<Integer> leftoutcolumns = new ArrayList<Integer>(); 
		
		if(!entropycalc||entropyString==null||(entropyString!=null&&entropyString.trim().equals(""))){
			for(int i = 0; i < sequences.get(0).length(); i++){
				//iterates column number
				int gapscounter=0;
				for(int j=0; j<nrseq; j++){
					//iterates sequences
					thiscolumn[j]=(char) sequences.get(j).charAt(i);
					if(MatrixTools.isGap(thiscolumn[j])){gapscounter++;}
					//System.out.println("Calculating column nr. " + j + " taking from seq " + j + " charAt " + i + " which is " + thiscolumn[j]);
				}
				boolean isoktoremove = true; 
				fullcolumns.add(thiscolumn.clone());

				//May only remove a column if it is not associated with exp data.
				isoktoremove = true;
				if(!extradata.isEmpty()){
					for(int k = 0; k<extradata.size(); k++){
						if(!extradata.get(k).isEmpty(i)){
							isoktoremove = false;
							break;
						}
					}
				}

				if((float)gapscounter/nrseq >= 0.25 && isoktoremove){		//syang: gaps column with no extraData will be assigned as '.' and left out
					finalstructure[i] = '.'; 
					finalreliability[i] = 0;
					leftoutcolumns.add(i);
					//System.out.println("Leaving out column: " + i);
				}
				else{
					columns.add(thiscolumn.clone());
				}

			}
		}
		else{
			//If entropy is calculated after a specific sequence, then remove columns where there are gaps in 
			//other sequences
			System.out.println("Removing gaps after sequence name: " + entropyString);			
			//Find the sequence: 
			int id = -1;
			for(id = 0; id<names.size(); id++){
				if(entropyString.trim().equals(names.get(id).trim())){
					break;
				}
			}
			if(id==-1){
				throw new Exception("Sequence with name " + entropyString + " could not be found!");
			}
			
			//Sequences.get(id) is the one that shouldn't have gaps 
			for(int i = 0; i < sequences.get(0).length(); i++){
				//iterates column number
				for(int j=0; j<nrseq; j++){
					//iterates sequences
					thiscolumn[j]=(char) sequences.get(j).charAt(i);
				}
				fullcolumns.add(thiscolumn.clone());
				//Remove column 
				if(MatrixTools.isGap(sequences.get(id).charAt(i))){
					finalstructure[i] = '.'; 
					finalreliability[i] = 0;
					leftoutcolumns.add(i);
				}
				else{
					//Keep column.
					columns.add(thiscolumn.clone());
				}
			}
		}
		
		//Remove left out columns from Extra data as well 
		for(ExtraData extrad:extradata){
			extrad.removeColumns(leftoutcolumns);
		}
		System.out.println("Total number of columns predicting structure on: " + columns.size());
		
		
		
		
		//Convert character array columns to integer array columns
		List<int[]> columns_int = new ArrayList<int[]>();
		for (int i = 0; i < fullcolumns.size(); i++) {
			columns_int.add(MatrixTools.convertColumn(fullcolumns.get(i)));
		}

		Progress activity = progress.getChildProgress(0.05);
		if(tree==null){				//syang: build NJ tree
			try {				
			//System.out.println("Creating tree by neighbour joining... ");
				tree = NeighbourJoining.generateTreeNJ(activity, sequences,columns_int, names, param);
				//tree.print();
			}catch (InterruptedException e) {
				System.out.println("Process interrupted by user! Stopping...");
				executor.shutDown();
				while(!executor.isTerminated()){
					Thread.sleep(100);
				}
				throw new InterruptedException();
			}
		}
		activity.setProgress(1.0);
		
		System.out.println("Checking data...");

		boolean isInputOK = FoldingProject.checkInput(tree, columns, names);
		if(!isInputOK){
			throw new IllegalArgumentException("There are problems with the input. Check that: " +
					"\n - Alignment rows all have the same length. " +
					"\n - The names of nodes in the tree match the names in the alignment. " +
			"\n - There are no illegal characters in any of the sequences."); 
		}
		System.out.println("Checking finished");
		
		Progress activity2 = progress.getChildProgress(0.05);
		if(optimizetree){		//syang: optimize tree branch lengths
			//System.out.println("Optimizing branch lengths...");
			try {
				tree.optimizeBranchLengths(activity2, columns_int,fullcolumns,names,param, iterlimit);		//syang: iterlimit==10
			} catch (InterruptedException e) {
				System.out.println("Process interrupted by user! Stopping...");
				executor.shutDown();
				while(!executor.isTerminated()){
					Thread.sleep(100);
				}
				throw new InterruptedException();
			}
			//if(!onlyCT){
			
//syang: START of commenting out file outputs		
/*			
				ExportTools.writeTree(outputdir, exportfilehandle, ".newick", tree);
*/
//syang: END of commenting out file outputs				
			
			//}
				tree.print();
		}
		activity2.setProgress(1.0);
		
		//Leave 2% for outputting
		Progress activity3 = progress.getChildProgress(0.88);
		ResultBundle result = null;
		try{
			result = FoldingProject.fold(activity3, phylodivisions,scfgdivisions, tree, 		//syang: the folding job
				columns, names, param, executor, verbose,1,extradata,false,entropycalc);
		}
		catch(InterruptedException e){
			System.out.println("Process interrupted by user! Stopping...");
			executor.shutDown();
			while(!executor.isTerminated()){
				Thread.sleep(100);
			}
			throw new InterruptedException();
		}
		activity3.setProgress(1.0);
		
		//ExportTools.writeCTFormat(outputdir, exportfilehandle+"_RED", ".ct", columns, result.getStructure(), result.getReliability());
		//ExportTools.writeStructureReliability(outputdir, exportfilehandle+"_RED", ".st", result.getStructure(), result.getReliability());

	
		
//syang: START of commenting out file outputs		
/*
		System.out.println("Structure prediction complete. Transforming to include left-out columns " +
				"(total " + leftoutcolumns.size() + ")...");

		progress.setCurrentActivity("Finalizing...");
		int icolcnt = 0;
		for(int i = 0; i<result.getBasePairProb().length; i++){
			while(leftoutcolumns.indexOf(i+icolcnt)!=-1){
				//	System.out.println("i = " + i);
				icolcnt++;
				finalmatrix[i][i] = 0;
				finalexp[i][i] = 0;
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
				finalexp[i+icolcnt][j+jcolcnt] = result.getExpectation()[i][j]; 

			}			

		}
		
		//syang: get single prob
		for(int i = 0; i<result.getSingleBaseProb().length; i++){
			finalsingleProb[i][0] = result.getSingleBaseProb()[i];
		}

		//System.out.println(finalmatrix[10][19]);
		//System.out.println(result.getBasePairProb()[10][17]);

		//could maybe implement a search-based solution using indexOf.
		//Dunno which one is better/faster.
		//This seems to work. 

		Iterator<Integer> iterator = leftoutcolumns.iterator();
		int lastone = 0;
		int cntr = 0; //counts how many columns are shifted
		int next = 0;
		while(next<sequences.get(0).length()){
			int co; //limit to which to go
			if(iterator.hasNext()){co=iterator.next();}
			else{co=sequences.get(0).length();}
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
		
		progress.setCurrentActivity("Exporting results...");
		ExportTools.writeCTFormat(outputdir, exportfilehandle, ".ct", fullcolumns, finalstructure, finalreliability);

		
		if(!onlyCT){
			ExportTools.writeFastaStructure(outputdir, exportfilehandle, ".seq", finalstructure, fullcolumns, names);
			ExportTools.writeStructureReliability(outputdir, exportfilehandle, ".st", finalstructure, finalreliability);
			ExportTools.writeFullFastaStructure(outputdir, exportfilehandle, ".lseq", finalstructure, fullcolumns, names);
		}

		
		if(exportson){
			ExportTools.writeTabbedMatrix(outputdir, exportfilehandle, ".ev", finalexp);	
			ExportTools.writeTabbedMatrix(outputdir, exportfilehandle, ".bp", finalmatrix);
			ExportTools.writeTabbedMatrix(outputdir, exportfilehandle, ".sbp", finalsingleProb);	//syang: get single prob
		}	
		
		if(seqexportname!=null){
			//Export sequence-specific files. 
			String filehandle = exportfilehandle + "-" + seqexportname;
			//Find the special sequence ID
			int seqID = -1;
			if(seqexportname!=null){
				for(int i = 0; i<names.size(); i++){
					if(seqexportname.trim().equals(names.get(i).trim())){
						seqID = i;
						break;
					}
				}
				if(seqID==-1){
					throw new IllegalArgumentException("Sequence with name " + seqexportname + " could not be found!");
				}
			}
			
			char[] reducedStructure = ExportTools.reducestructure(sequences.get(seqID), finalstructure).toCharArray();
			char[] reducedSequence = ExportTools.reducesequence(sequences.get(seqID)).toCharArray();
			float[] reducedProbabilities = ExportTools.reducereliabilities(sequences.get(seqID), finalreliability);
			
			List<char[]> seqcolumns = new ArrayList<char[]>();
			for(char c:reducedSequence){
				seqcolumns.add(new char[]{c});
			}

			ExportTools.writeCTFormat(outputdir, filehandle, ".ct", seqcolumns, reducedStructure, reducedProbabilities);
			ExportTools.writeStructureReliability(outputdir, filehandle, ".st", reducedStructure, reducedProbabilities);
			
		}
*/
//syang: END of commenting out file outputs	
		
		
		
		progress.setCurrentActivity("Folding finished.");
		progress.setProgress(1.0);
	}


	public boolean success() {
		return this.success;
	}
	
	public static void setProgressBar(PPfoldProgressBar progressin){
		progress = progressin;
	}

}
