package com.ppfold.algo;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.ppfold.algo.extradata.ExtraData;
import com.ppfold.algo.extradata.ForcedConstraints;

/**
 * After reading and parsing input, algorithm execution starts here. (Note:
 * Columns with many gaps are already removed)
 * 
 * @author Z.Sukosd
 */
public class FoldingProject {
	static final double LOG_TWO = Math.log(2d);

	/**
	 * Makes sure that all leaves of the tree can be found, and creates a list
	 * of leaves for easy access.
	 */
	public static boolean checkInput(Tree tree, List<char[]> columns,
			List<String> names) {
		boolean treematchesnames = true;
		boolean columnshaverightlengths = true;
		boolean noillegalcharacters = true;

		int firstsize = columns.get(0).length;
		for (int col = 0; col < columns.size(); col++) {
			if (columns.get(col).length != firstsize) {
				columnshaverightlengths = false;
				System.err
						.println("Columns in alignment don't have same length!");
			}
			for (int row = 0; row < columns.get(col).length; row++) {
				char thischar = columns.get(col)[row];
				thischar = Character.toLowerCase(thischar);
				if (thischar != 'a' && thischar != 'u' && thischar != 't'
						&& thischar != 'g' && thischar != 'c'
						&& thischar != 'r' && thischar != 'y'
						&& thischar != 's' && thischar != 'w'
						&& thischar != 'k' && thischar != 'm'
						&& thischar != 'b' && thischar != 'd'
						&& thischar != 'h' && thischar != 'v'
						&& thischar != 'n'
						&& !MatrixTools.isGap(thischar)) {
					System.err.println("Illegal character! " + thischar);
					noillegalcharacters = false;
				}
			}
		}

		for (int row = 0; row < columns.get(0).length; row++) {
			// finds the node corresponding to the rownumber
			// rownumber = corresponds to sequences.
			Node node = tree.findSlowlyNodeWithName(names.get(row));
			if (node == null) {
				System.err.println("Can't find node with name "
						+ names.get(row) + "!");
				treematchesnames = false;
			}
		}
		return treematchesnames && columnshaverightlengths
				&& noillegalcharacters;
	}

//	public static ResultBundle fold(int phylojobsnr, int scfgjobsnr, Tree tree,
//			List<char[]> columns, List<String> names, Parameters param,
//			AsynchronousJobExecutor executor, boolean verbose, int execnr,
//			List <ExtraData> extradata_list, boolean diffbp, boolean entropycalc) {
//		try {
//			return fold(NullProgress.INSTANCE, phylojobsnr, scfgjobsnr, tree,
//					columns, names, param, executor, verbose, execnr, extradata_list, diffbp, entropycalc);
//		} catch (InterruptedException e) {
//			// Never happens - we use the NullActivity
//			e.printStackTrace();
//			return null;
//		}
//	}

	public static ResultBundle fold(Progress act, int phylojobsnr,
			int scfgjobsnr, Tree tree, List<char[]> columns,
			List<String> names, Parameters param,
			AsynchronousJobExecutor executor, boolean verbose, int execnr,
			List <ExtraData> extradata_list, boolean diffbp, boolean entropycalc)
			throws InterruptedException {
		
		if(columns.size()<2){
			return ResultBundle.tinyBundle();
		}
		
		double scfg_to_phylo_ratio = 0.095 * columns.size()
				/ columns.get(0).length;
		double phylopart = 0.95 / (scfg_to_phylo_ratio + 1);
		if (phylopart < 0 | phylopart > 1) {
			System.err.println("Time estimation resulted "
					+ "in illegal value for phylogenetic calculations: "
					+ phylopart);
			phylopart = 0.475; // in case the estimation results in a
			// "weird number", just set it to 50-50%.
		}
		double scfgpart = 0.95 - phylopart; // max is 95% because processing
		// must happen before and after
		act.setCurrentActivity("Applying evolutionary model");

		
		//long starttime = System.currentTimeMillis();
		double[][] probmatrix = FoldingProject.createPhyloProb(act
				.getChildProgress(phylopart), phylojobsnr, tree, columns,
				names, columns.size(), param, executor, verbose, execnr);

		//If there are extra data, multiply probmatrix with the relevant values
		//syang: probmatrix[63][63] for gca-alignment_sy.fasta, so probmatrix.length==probmatrix[0].length==68-5, i.e.
		//remove columns where there are gaps in them
		//for gca-alignment_few_sy.fasta, probmatrix[7][7]
		//System.out.println("\nsyang debug: "+probmatrix.length + " "+probmatrix[0].length);
		
		double[][] probmatrix2 = null;
		if(diffbp&&extradata_list!=null){
			probmatrix2= new double[probmatrix.length][probmatrix[0].length];
			for (int a = 0; a<probmatrix2.length;a++){
				for(int b = 0;b<probmatrix2[0].length; b++){
					probmatrix2[a][b] = probmatrix[a][b]; 
				}
			}
		}

		/*		
		//syang: probmatrix is basepairing prob matrix
		for(int i=0; i<probmatrix.length;i++){
			for(int j=0;j<probmatrix[0].length;j++){
				System.out.print(probmatrix[i][j]+" ");
			}
			System.out.println();
		}
		System.out.println();
		//probmatrix2 is null since diffbp is set to be false in PPfoldMain.java
		for(int i=0; i<probmatrix2.length;i++){
			for(int j=0;j<probmatrix2[0].length;j++){
				System.out.print(probmatrix2[i][j]+" ");
			}
			System.out.println();
		}*/
		
		
		
		//syang: save another copy of probmatrix[][] for constraint folding
		double[][] probmatrix0 = null;
		//if(extradata_list!=null){
			probmatrix0= new double[probmatrix.length][probmatrix[0].length];
			for (int a = 0; a<probmatrix0.length;a++){
				for(int b = 0;b<probmatrix0[0].length; b++){
					probmatrix0[a][b] = probmatrix[a][b]; 
				}
			}
		//}		
		
		
		//syang:
		// do not allow basepairs less than 4 nucleotides apart.
		for (int i = 0; i < probmatrix0.length; i++) {
			for (int j = 1; j < 4; j++) {
				if (i + j < probmatrix0.length) {
					probmatrix0[i][i + j] = 0;
				}
				if (i - j > 0) {
					probmatrix0[i][i - j] = 0;
				}
			}
		}
		
		
		// syang: the result of noConstraints topInsideS
		System.out.println("Folding (noConstraints)...");
		act.setCurrentActivity("Applying grammar");
		PointRes topInsideS = FoldingProject.calcSCFG(act
				.getChildProgress(scfgpart), scfgjobsnr, param.getProb(),
				probmatrix0, probmatrix2, executor, verbose, diffbp, entropycalc);
		
		
		if(extradata_list!=null){
			System.out.println("Number of auxiliary data items: " + extradata_list.size());
			ArrayList<ForcedConstraints> fcList = new ArrayList<ForcedConstraints>();
			for(ExtraData data:extradata_list){
				if(data instanceof ForcedConstraints){			//syang: pick ForcedConstraints out from ExtraData?
					//If forcing constraints, collect them into one and remove from the list
					fcList.add((ForcedConstraints) data);
				}
			}
			for(ExtraData data:fcList){
				extradata_list.remove(data);
			}
			if(fcList.size()>0){
				System.out.println("Processing " + fcList.size () + " hard constraint datasets");
				ForcedConstraints fcMerged = new ForcedConstraints();
				fcMerged = fcMerged.combinedForcedConstraints(probmatrix.length, fcList);
				extradata_list.add(fcMerged);
			}

		}
		
		System.out.println("Processing all auxiliary data");
		if(extradata_list!=null){
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix[i][j] *= extradata.getProbabilityGivenInnerPaired(i,j);
						probmatrix[j][i] = probmatrix[i][j];
					}
					probmatrix[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
			/*
			//syang: probmatrix is basepairing prob matrix
			for(int i=0; i<probmatrix.length;i++){
				for(int j=0;j<probmatrix[0].length;j++){
					System.out.print(probmatrix[i][j]+" ");
				}
				System.out.println();
			}
			System.out.println();				
			//probmatrix2 is null since diffbp is set to be false in PPfoldMain.java
			for(int i=0; i<probmatrix2.length;i++){
				for(int j=0;j<probmatrix2[0].length;j++){
					System.out.print(probmatrix2[i][j]+" ");
				}
				System.out.println();
			}*/
		}
		
		if(diffbp&&extradata_list!=null){
			//If inner and outer basepairs are distinguished
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix2.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix2[i][j] *= extradata.getProbabilityGivenOuterPaired(i,j);
						probmatrix2[j][i] = probmatrix2[i][j];
					}
					probmatrix2[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}			
		}		
		
		
		//System.out.println("Time in phylogenetic part: " + (System.currentTimeMillis()-starttime));
		// do not allow basepairs less than 4 nucleotides apart.
		for (int i = 0; i < probmatrix.length; i++) {
			for (int j = 1; j < 4; j++) {
				if (i + j < probmatrix.length) {
					probmatrix[i][i + j] = 0;
				}
				if (i - j > 0) {
					probmatrix[i][i - j] = 0;
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			for (int i = 0; i < probmatrix2.length; i++) {
				for (int j = 1; j < 4; j++) {
					if (i + j < probmatrix2.length) {
						probmatrix2[i][i + j] = 0;
					}
					if (i - j > 0) {
						probmatrix2[i][i - j] = 0;
					}
				}
			}
		}
		

		System.out.println("Folding (Constraints)...");
		act.setCurrentActivity("Applying grammar");
		
		// syang: the result of Constraint topInsideS
		//syang: param is the initial scfg, rate matrix and so
		//starttime = System.currentTimeMillis();
		PointRes topInsideS_constraint = FoldingProject.calcSCFG(act
				.getChildProgress(scfgpart), scfgjobsnr, param.getProb(),
				probmatrix, probmatrix2, executor, verbose, diffbp, entropycalc);
		//System.out.println("Time in SCFG part: " + (System.currentTimeMillis()-starttime));
		
		
		
		//syang: compute the final motif-wise probability
		System.out.println("Top inside S variable:"+
		"\nbefore constraint: "+topInsideS+"\nafter constraint: " +topInsideS_constraint);
		topInsideS_constraint.divide(topInsideS);
		System.out.println("The quotient (i.e. motif-wise prob):"+topInsideS_constraint+
				" = "+topInsideS_constraint.toFloat());
		
		
		
		//Shut down the executor so we aren't hanging at the end
		executor.shutDown();
		
		ResultBundle result=new ResultBundle();
		return result;

	}

	private static PointRes calcSCFG(Progress act, int userjobsnr,
			double[][] prob, double[][] probmatrix, double[][] probmatrix2,
			AsynchronousJobExecutor executor, final boolean verbose, final boolean diffbp, final boolean entropycalc)
			throws InterruptedException {
		PointRes tmp = new PointRes(0, 0);		//syang: first 0 is fraction, second 0 is exponent
		final long starttime = System.nanoTime();

		if (verbose) {
			System.out.println("Timer (SCFG) started. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
			System.out.println("Processing input...");
		}

		act.setCurrentActivity("Applying grammar: processing input");
		

		int length = probmatrix.length;
		
		//syang: userjobsnr==32 for gca-alignment_sy.fasta and gca-alignment_few_sy.fasta
		//System.out.println("userjobsnr:"+userjobsnr);
		
		// First corrections of user input
		//syang: userjobsnr corrected to 6 for gca-alignment_few_sy.fasta
		if (userjobsnr > length) {
			userjobsnr = length - 1;
		} else if (userjobsnr < 1) {
			userjobsnr = 1;
		} else {
		}

		// Now calculate number of jobs
		int distance = length / userjobsnr + 1; // the +1 is to prevent
		// generation of superfluous
		// jobs
		
		//syang: distance==2 for gca-alignment_sy.fasta and gca-alignment_few_sy.fasta
		//System.out.println("distance:"+distance);
		
		int firstrowspaces = length + 1;

		int isthereextra = (firstrowspaces % distance != 0) ? 1 : 0; 
		// triangle is compelete iff isthereextra = 0.
		// syang: nrdivisions==32 for gca-alignment_sy.fasta; 4 for gca-alignment_few_sy.fasta
		int nrdivisions = ((firstrowspaces - firstrowspaces % distance) / distance)
				+ isthereextra;
		// nrdivisions = how many jobs should there be in 1st row?
		// syang: nrsectors==528 for gca-alignment_sy.fasta; 10 for gca-alignment_few_sy.fasta
		int nrsectors = nrdivisions * (nrdivisions + 1) / 2; 
		// total number of jobs in full triangle

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("User wish for number of divisions: "
					+ userjobsnr);		//syang: 6 for gca-alignment_few_sy.fasta
			System.out
					.println("Number of sectors in inside-outside calculations: "
							+ nrsectors);		
			System.out
					.println("Distance for inside-outside has been calculated as: "
							+ distance);
			System.out
					.println("Corrected number of jobs in first row of inside-outside: "
							+ nrdivisions);		
			System.out
					.println("Total number of inside-outside jobs to generate: "
							+ nrsectors);
		}

		// syang: interestingpoints==1984 for gca-alignment_sy.fasta; 24 for gca-alignment_few_sy.fasta
		int interestingpoints = length * length / 2;
		// syang: totalpoints==2112 for gca-alignment_sy.fasta; 40 for gca-alignment_few_sy.fasta
		int totalpoints = nrsectors * distance * distance;
		// syang: extrapoints==128 for gca-alignment_sy.fasta; 16 for gca-alignment_few_sy.fasta
		int extrapoints = totalpoints - interestingpoints;
		// syang: fraction==6.451613% for gca-alignment_sy.fasta; 66.666664% for gca-alignment_few_sy.fasta
		float fraction = (float) extrapoints * 100 / interestingpoints;
		if(verbose){
			System.out.println("Divisions = " + nrdivisions
				+ ", Interesting points: " + interestingpoints
				+ ", Extra points: " + extrapoints + ", Total points: "
				+ totalpoints + ", " + "Fractional extra: " + fraction + "%");
		}
		if (verbose) {
			System.out.println("Generating sectors: ");
		}
		act.setCurrentActivity("Applying grammar: generating sectors");

		Sector top = SectorGenerator.GenerateSectors(nrsectors, distance,
				nrdivisions, length,diffbp);

		final Master master = new Master(top, prob);

		// sector: set basepairs (phylogenetic probabilities) for inside sectors
		PointRes number = new PointRes(0, 0);
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length - i; j++) {
				Sector sectoro = findSector(i, j, master.bottom);	//syang: find which sector point(i,j) belongs to; Note the shape of the sector(Fig3 in paper)
				int so = findPointST(i, j, sectoro, distance)[0];
				int to = findPointST(i, j, sectoro, distance)[1];
				number.setToDouble(probmatrix[i][i + j]);
				sectoro.setBasePairs(so, to, number);
			}
		}
		
		if(diffbp&&probmatrix2!=null){
			//if differentiating inner and outerbasepairs, create a new bp matrix per sector
			for (int i = 0; i < length; i++) {
				for (int j = 0; j < length - i; j++) {
					Sector sectoro = findSector(i, j, master.bottom);
					int so = findPointST(i, j, sectoro, distance)[0];
					int to = findPointST(i, j, sectoro, distance)[1];
					number.setToDouble(probmatrix2[i][i + j]);
					sectoro.setBasePairs2(so, to, number);
				}
			}
		}
		
		
		
		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}


		act.setCurrentActivity("Applying grammar: inside algorithm");
		
		if (verbose) {
			System.out.println("Calculating inside values... ");
		}

		final AtomicInteger finishedinsidejobscount = new AtomicInteger(0); // counts
		// how many inside jobs are done
		master.CreateInsideJobChannel();
		Progress insideAct = act.getChildProgress(0.32);
		//final long gridstarttime = System.nanoTime();
		while (master.unProcessedInsideSectors()) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			CYKJob cYKJob = master.takeNextInsideJob(); // This call will block
			// until a job is ready
			final Progress jobAct = insideAct.getChildProgress(1.0 / nrsectors);
			// System.out.println(cYKJob.sectorid + " " + " 0 " +
			// (System.nanoTime()-starttime));
			final int sectorNumber = cYKJob.getSectorid();
			executor.startExecution(cYKJob, new JobListener() {
				public void jobFinished(JobResults result) {
					master.setInsideResult(sectorNumber, result);
					jobAct.setProgress(1.0);
					// System.out.println(sectorNumber + " " + " 1 " +
					// (System.nanoTime()-starttime));
					finishedinsidejobscount.incrementAndGet();
				}

				public void jobFinished(double[][] result) {
				}// doesn't happen here

				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here
			});
		}

		// wait for last job to finish
		while (finishedinsidejobscount.get() < nrsectors) {
			Thread.sleep(100);
			if(act.shouldStop()){	
				executor.shutDown();
			}
			act.checkStop();
		}

		insideAct.setProgress(1.0);

		if(verbose){
			System.out.println("Top inside: "
				+ master.top.getInsideMatrixS().getProb(distance - 1,
						length - 1 - master.top.pos[1]));
		}
		
/*		//syang: test the change of prob[][]
		System.out.println("for syang, test the change of prob");
		for(int i=0;i<prob.length;i++){
			for(int j=0;j<prob[0].length;j++){
				System.out.print("i,j: "+prob[i][j]+", ");
			}
		}
		System.out.println();
*/		
		//syang: directly output the top inside value
//		System.out.println("For syang, Top inside: "
//				+ master.top.getInsideMatrixS().getProb(distance - 1,
//						length - 1 - master.top.pos[1]));
		
/*		//syang: test
		System.out.println("For syang, distance: "
				+ distance + "; length:" + length + "; master.top.pos[1]:" 
				+ master.top.pos[1] + "; n:" 
				+ master.top.getInsideMatrixS().get_n() + "; dim:" 
				+ master.top.dim);
		for(int i=0; i< master.top.dim; i++){
			for (int j=0; j< master.top.dim; j++){
				System.out.println("For syang, " + i + "," +j + ": "
						+ master.top.getInsideMatrixS().getProb(i,j));
			}
		}
*/

/*		//syang: print the entire inside matrix
 * 		System.out.println("========\nFor syang: inside");
		Sector bottomRight = master.bottom;
		PointRes matSy[][]=new PointRes[7][7];
		while(bottomRight!=null){
			System.out.println("#"+bottomRight.sectorid);
			for(int i=0;i<bottomRight.dim;i++){
				for(int j=0;j<bottomRight.dim;j++){
					System.out.print("("+i+","+j+")="+bottomRight.getInsideMatrixF().getProb(i,j)+"\t");
					
					int x=distance+bottomRight.pos[0]-1-i;
					int y=j-x+bottomRight.pos[0]+bottomRight.pos[1];
					if(x>=0 && y>=0){
						matSy[y][x]=bottomRight.getInsideMatrixF().getProb(i,j);
					}
				}
				System.out.println();
			}
			
			bottomRight=bottomRight.next;
		}
		
		for(int i=matSy.length-1;i>=0;i--){
			for(int j=0;j<matSy.length-i;j++){
				System.out.print(matSy[i][j]+"\t");
			}
			System.out.println();
		}
*/
/*		//syang: print the entire basePair matrix
		System.out.println("========\nFor syang: inside");
		Sector bottomRight = master.bottom;
		PointRes matSy[][]=new PointRes[7][7];
		while(bottomRight!=null){
			System.out.println("#"+bottomRight.sectorid);
			for(int i=0;i<bottomRight.dim;i++){
				for(int j=0;j<bottomRight.dim;j++){
					System.out.print("("+i+","+j+")="+bottomRight.getBasePairs().getProb(i,j)+"\t");
					
					int x=distance+bottomRight.pos[0]-1-i;
					int y=j-x+bottomRight.pos[0]+bottomRight.pos[1];
					if(x>=0 && y>=0){
						matSy[y][x]=bottomRight.getBasePairs().getProb(i,j);
					}
				}
				System.out.println();
			}
			
			bottomRight=bottomRight.next;
		}
		
		for(int i=matSy.length-1;i>=0;i--){
			for(int j=0;j<matSy.length-i;j++){
				System.out.print(matSy[i][j]+"\t");
			}
			System.out.println();
		}
*/		
		
		//long insidetime = System.nanoTime();
		//System.out.print(nrdivisions + " - ");
		//System.out.println("TOTAL TIME ELAPSED IN INSIDE PART (ALL): "
		//		+ (insidetime - starttime) * 1e-9 + " seconds ");
		//double gridtime1 = (insidetime - gridstarttime) * 1e-9;
		//System.out.println("TOTAL TIME ELAPSED IN INSIDE PART (DISTRIBUTED): "
		//		+ gridtime1 + " seconds ");

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		
		//syang: START of commenting out outside algorithm
		/*
		act.setCurrentActivity("Applying grammar: outside algorithm");
		
		if (verbose) {
			System.out.println("Calculating outside values...");
		}
		master.CreateOutsideJobChannel();
		final AtomicInteger finishedoutsidejobscount = new AtomicInteger(0); // counts
		// how many outside jobs are done

		// outside algorithm
		// set basepairs for outside algorithm (they are shifted relative to
		// inside algo)
		number.setToFloat(0);
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length - i; j++) {
				Sector sectoro = findSector(i, j, master.bottom);
				int so = findPointST(i, j, sectoro, distance)[0];
				int to = findPointST(i, j, sectoro, distance)[1];
				if (i > 0 && (j + i + 1) < length) {
					number.setToDouble(probmatrix[i - 1][i + j + 1]);
				} else if (i == 0) {
					number.setToFloat(0);
				} else if (j + i + 1 == length) {
					number.setToFloat(0);
				}
				sectoro.setBasePairs(so, to, number);
			}
		}
		
		if(diffbp){
			for (int i = 0; i < length; i++) {
				for (int j = 0; j < length - i; j++) {
					Sector sectoro = findSector(i, j, master.bottom);
					int so = findPointST(i, j, sectoro, distance)[0];
					int to = findPointST(i, j, sectoro, distance)[1];
					if (i > 0 && (j + i + 1) < length) {
						number.setToDouble(probmatrix2[i - 1][i + j + 1]);
					} else if (i == 0) {
						number.setToFloat(0);
					} else if (j + i + 1 == length) {
						number.setToFloat(0);
					}
					sectoro.setBasePairs2(so, to, number);
				}
			}
		}
		
		//final long outsidegridstarttime = System.nanoTime();
		Progress outsideAct = act.getChildProgress(0.50);
		while (master.unProcessedOutsideSectors()) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			final Progress jobAct = outsideAct
					.getChildProgress(1.0 / nrsectors);
			CYKJob cYKJob = master.takeNextOutsideJob(); // This call will block
			// until a job is
			// ready
			final int sectorNumber = cYKJob.getSectorid();
			executor.startExecution(cYKJob, new JobListener() {
				// System.out.println(sectorNumber + " " + " 2 " +
				// (System.nanoTime()-starttime));
				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here

				public void jobFinished(JobResults result) {
					master.setOutsideResult(sectorNumber, result);
					jobAct.setProgress(1.0);
					// System.out.println(sectorNumber + " " + " 3 " +
					// (System.nanoTime()-starttime));
					finishedoutsidejobscount.incrementAndGet();
				}

				public void jobFinished(double[][] result) {
				}// doesn't happen here
			});
		}

		// wait for last job to finish
		while (finishedoutsidejobscount.get() < nrsectors) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		outsideAct.setProgress(1.0);
		
				
		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
*/	//syang: END of commenting out outside algorithm
		
		
		
		
		
		
		
/*		//syang: directly output the top outside value
		System.out.println("For syang, Bottom outside: "
				+ master.bottom.getOutsideMatrixL().getProb(
						1,1));
				System.out.println("For syang, Bottom outside: "
						+ master.bottom.getOutsideMatrixL().getProb(
								length - 1 - master.bottom.pos[1],distance - 1));
*/				
		//syang: test
//		System.out.println("For syang, distance: "
//				+ distance + "; length:" + length + "; master.top.pos[1]:" 
//				+ master.top.pos[1] + "; n:" 
//				+ master.top.getInsideMatrixS().get_n() + "; dim:" 
//				+ master.top.dim);
/*	for(int i=0; i< master.top.dim; i++){
				for (int j=0; j< master.top.dim; j++){
					System.out.println("For syang, " + i + "," +j + ": "
							+ master.top.getOutsideMatrixL().getProb(i,j));
				}
			}
		System.out.println("========");*/
/*		for(int i=0; i< master.bottom.dim; i++){
			for (int j=0; j< master.bottom.dim; j++){
				System.out.println("For syang, bottomInsideS " + i + "," +j + ": "
						+ master.bottom.getInsideMatrixS().getProb(i,j));
			}
		}
		System.out.println("========");
		for(int i=0; i< master.bottom.dim; i++){
			for (int j=0; j< master.bottom.dim; j++){
				System.out.println("For syang, bottomOutsideL " + i + "," +j + ": "
						+ master.bottom.getOutsideMatrixL().getProb(i,j));
			}
		}	
		
		PointRes bottomOutsideL = master.bottom.getOutsideMatrixL().getProb(1,1);
		bottomOutsideL.multiply(master.bottom.getInsideMatrixL().getProb(1,1));
		PointRes bottomOutsideS = master.bottom.getOutsideMatrixS().getProb(1,1);
		bottomOutsideS.multiply(master.bottom.getInsideMatrixS().getProb(1,1));
		PointRes bottomOutsideF = master.bottom.getOutsideMatrixF().getProb(1,1);
		bottomOutsideF.multiply(master.bottom.getInsideMatrixF().getProb(1,1));
		
		System.out.println("L: "+bottomOutsideL+"\nS: "+bottomOutsideS+"\nF: "+bottomOutsideF);
		bottomOutsideL.add(bottomOutsideS);
		bottomOutsideL.add(bottomOutsideF);
		System.out.println("\nSum: "+bottomOutsideL);
*/		
/*		System.out.println("========\nFor syang:");
		Sector bottomRight = master.top.below.below.below.below;//master.bottom.above.above.above.above.above.above.above.above.above.above;
		//while(bottomRight.next.below.below.below.below.below.below.below.below.below.below.below==null){
		while(bottomRight!=null){
			PointRes bottomOutsideL = bottomRight.getOutsideMatrixL().getProb(1,1);
			bottomOutsideL.multiply(bottomRight.getInsideMatrixL().getProb(1,1));
			PointRes bottomOutsideS = bottomRight.getOutsideMatrixS().getProb(1,1);
			bottomOutsideS.multiply(bottomRight.getInsideMatrixS().getProb(1,1));
			PointRes bottomOutsideF = bottomRight.getOutsideMatrixF().getProb(1,1);
			bottomOutsideF.multiply(bottomRight.getInsideMatrixF().getProb(1,1));
			
			System.out.println("\n#"+bottomRight.sectorid+"\nL: "+bottomOutsideL+"\nS: "+bottomOutsideS+"\nF: "+bottomOutsideF);
			bottomOutsideL.add(bottomOutsideS);
			bottomOutsideL.add(bottomOutsideF);
			System.out.println("Sum: "+bottomOutsideL);
			bottomRight = bottomRight.next;
		}
*/		
/*		for(int i=0; i< bottomRight.dim; i++){
			for (int j=0; j< bottomRight.dim; j++){
				System.out.println("For syang, " + i + "," +j + ": "
						+ bottomRight.getOutsideMatrixL().getProb(i,j));
			}
		}
*/
/*		//syang: print the entire outside matrix
  		System.out.println("========\nFor syang: outside");
		Sector bottomRight = master.bottom;
		PointRes matSy[][]=new PointRes[7][7];
		while(bottomRight!=null){
			System.out.println("#"+bottomRight.sectorid);
			for(int i=0;i<bottomRight.dim;i++){
				for(int j=0;j<bottomRight.dim;j++){
					System.out.print("("+i+","+j+")="+bottomRight.getOutsideMatrixS().getProb(i,j)+"\t");
					
					int x=distance+bottomRight.pos[0]-1-i;
					int y=j-x+bottomRight.pos[0]+bottomRight.pos[1];
					if(x>=0 && y>=0){
						matSy[y][x]=bottomRight.getOutsideMatrixS().getProb(i,j);
					}
				}
				System.out.println();
			}
			
			bottomRight=bottomRight.next;
		}
		
		for(int i=matSy.length-1;i>=0;i--){
			for(int j=0;j<matSy.length-i;j++){
				System.out.print(matSy[i][j]+"\t");
			}
			System.out.println();
		}
		
*/
/*		//syang: print the entire basePair matrix
		System.out.println("========\nFor syang: inside");
		Sector bottomRight4 = master.bottom;
		PointRes matSy4[][]=new PointRes[7][7];
		while(bottomRight4!=null){
			System.out.println("#"+bottomRight4.sectorid);
			for(int i=0;i<bottomRight4.dim;i++){
				for(int j=0;j<bottomRight4.dim;j++){
					System.out.print("("+i+","+j+")="+bottomRight4.getBasePairs().getProb(i,j)+"\t");
					
					int x=distance+bottomRight4.pos[0]-1-i;
					int y=j-x+bottomRight4.pos[0]+bottomRight4.pos[1];
					if(x>=0 && y>=0){
						matSy4[y][x]=bottomRight4.getBasePairs().getProb(i,j);
					}
				}
				System.out.println();
			}
			
			bottomRight4=bottomRight4.next;
		}
		
		for(int i=matSy4.length-1;i>=0;i--){
			for(int j=0;j<matSy4.length-i;j++){
				System.out.print(matSy4[i][j]+"\t");
			}
			System.out.println();
		}
*/				
		
				
		

		PointRes topInsideS = new PointRes(0,0);
		topInsideS=master.top.getInsideMatrixS().getProb(
				distance - 1, length - 1 - master.top.pos[1]);
		
		act.setProgress(1.0);
		return topInsideS;
//syang: END of commenting out finalizing			
		
		
	}

	private static double[][] createPhyloProb(Progress act, int userjobsnr,
			Tree tree, List<char[]> columns_char, List<String> names,
			final int length, Parameters param,
			AsynchronousJobExecutor executor, final boolean verbose, int execnr)
			throws InterruptedException {
		final long starttime = System.nanoTime();
		if(verbose){
		System.out.println("Timer (phylogeny) started. (time: "
				+ (System.nanoTime() - starttime) * 1e-9 + " s)");
		}
		if (verbose) {
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out.println("Processing input...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: processing input");
		}
		List<int[]> columns = new ArrayList<int[]>();
		for (int i = 0; i < columns_char.size(); i++) {
			columns.add(MatrixTools.convertColumn(columns_char.get(i)));
		}
		if (verbose) {
			System.out.println("User wish for number of divisions: "
					+ userjobsnr);
		}

		// correct user input
		if (userjobsnr > length) {
			userjobsnr = length;
		} else if (userjobsnr < 1) {
			userjobsnr = 1;
		}
		int nrjobs = userjobsnr;

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out
					.println("Actual nr. of divisions in phylogenetic calculations: "
							+ nrjobs);
		}
		if (verbose) {
			System.out.println("Generating matrices for each node...");
		}

		tree.generateLeafList(names); // to speed up calculations later

		final double[][] probmatrix = new double[length][length];
		List<PhyloJob> jobs = new ArrayList<PhyloJob>();
		// First generate "exp(Rt)" matrix for each node.
		// This is the same for all columns of that node.
		tree.getRoot().calculateChildrenMatrix(param.getSD(), param.getSV(),
				param.getSV1());

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("Generating jobs...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: dividing tasks");
		}

		// Create jobs for SINGLE columns
		// now all single columns are in one job
		// Add the only single-column job
		int colcnt = 0;
		PhyloJob lastjob = new PhyloJob();
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob (for
		// debugging only)
		lastjob.startcol = 0;
		lastjob.jobid = 0;
		lastjob.endcol = length - 1;
		lastjob.type = false;
		lastjob.param = param;
		for (int col = 0; col < length; col++) {
			// send columns
			lastjob.columns.add(columns.get(col));
			colcnt++;
		}
		jobs.add(lastjob);

		// Create jobs for column PAIRS
		// count total number of column pairs
		int paircnt = 0;
		for (int i = 0; i < length; i++) {
			for (int j = i + 1; j < length; j++) {
				paircnt++;
			}
		}
		if(verbose){
			System.out.println("Total number of pairs: " + paircnt);
			System.out.println("Pairs in a job: " + (paircnt) / nrjobs);
		}
		// have to calculate matrices, but only once
		tree.getRoot().calculateChildrenMatrix(param.getDD(), param.getDV(),
				param.getDV1());

		int lastendcol = -1;
		int paircolcnt = 0;
		// add all jobs except last
		for (int jobnr = 0; jobnr < nrjobs - 1; jobnr++) {
			PhyloJob job = new PhyloJob();
			job.tree = Tree.copyTree(tree); // must copy the entire tree into
			// each phylojob
			job.names = names; // must copy the names into each phylojob
			job.type = true;
			job.param = param;
			job.jobid = jobnr + 1;

			int col1 = lastendcol + 1;
			job.startcol = col1;
			paircolcnt += length - col1;
			col1--;
			while (paircolcnt < ((long)(jobnr + 1) * (paircnt)) / nrjobs) {
				col1++;
				if (col1 == length) { // no more columns left, finish
					col1--;
					break;
				}
				job.columns.add(columns.get(col1));
				paircolcnt += length - col1;
			}

			job.endcol = col1;
			lastendcol = col1;
			for (int col2 = job.startcol + 1; col2 < length; col2++) {
				// send pairing columns
				job.columns2.add(columns.get(col2));
			}
			jobs.add(job);
			paircolcnt -= length - col1;
		}
		// Add the last job
		lastjob = new PhyloJob();
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob
		lastjob.startcol = lastendcol + 1;
		lastjob.endcol = length;
		lastjob.type = true;
		lastjob.param = param;
		lastjob.jobid = nrjobs;
		for (int col1 = lastjob.startcol; col1 < length; col1++) {
			// send columns
			lastjob.columns.add(columns.get(col1));
		}
		for (int col2 = lastjob.startcol + 1; col2 < length; col2++) {
			// send pairing columns
			lastjob.columns2.add(columns.get(col2));
		}
		jobs.add(lastjob);

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out
					.println("Total number of jobs in phylogenetic calculations: "
							+ jobs.size());
		}
		if (verbose) {
			System.out.println("Executing jobs...");
		}

		// start executing the jobs
		// counts how many phylojobs are done		
		final AtomicInteger finishedphylojobscount = new AtomicInteger(0); 

		act.setCurrentActivity("Evolutionary model: " +
					"calculating column probabilities");
		
		long gridstarttime = System.nanoTime();
		// execute single column jobs
		Progress singleColAct = act.getChildProgress(0.1);
		for (int jobnr = 0; jobnr < 1; jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJob job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;

			final Progress jobAct = singleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) colcnt);
			executor.startExecution(job, new JobListener() {
				public void jobFinished(JobResults result) {
				} // doesn't happen here

				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here

				public void jobFinished(double[][] result) {
					for (int col = 0; col < result.length; col++) {
						probmatrix[col + startcol][col + startcol] = result[col][0];
					}
					finishedphylojobscount.incrementAndGet();
					jobAct.setProgress(1.0);
				}
			});
		}

		Progress doubleColAct = act.getChildProgress(0.9);

		// execute double column jobs
		for (int jobnr = 1; jobnr < jobs.size(); jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJob job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;
			final int endcol = job.endcol;
			final Progress jobAct = doubleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) length);
			if (job.columns2.size() != 0) {
				executor.startExecution(job, new JobListener() {
					public void jobFinished(JobResults result) {
					} // doesn't happen here

					public void jobFinished(List<ResultBundle> result) {
					} // doesn't happen here

					public void jobFinished(double[][] result) {
						for (int col1 = startcol; col1 < endcol + 1; col1++) {
							for (int col2 = col1 + 1; col2 < length; col2++) {
								probmatrix[col1][col2] = result[col1 - startcol][col2
										- startcol - 1];
								probmatrix[col2][col1] = probmatrix[col1][col2]; // make
								// it
								// symmetric
							}
						}
						finishedphylojobscount.incrementAndGet();
						jobAct.setProgress(1.0);
					}
				});
			} else {
				finishedphylojobscount.incrementAndGet();
				jobAct.setProgress(1.0);
			}
		}

		// wait for last thread to finish
		// wait for jobs to finish
		while (finishedphylojobscount.get() < jobs.size()) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		singleColAct.setProgress(1.0);
		doubleColAct.setProgress(1.0);
		act.setProgress(1.0);
		//act.setCurrentActivity("Evolutionary model: done");

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		System.out.println("TOTAL TIME ELAPSED IN PHYLOGENETIC PART: "
				+ (int)((System.nanoTime() - starttime) * 1e-9) + " seconds ");
		
		if(verbose){
			System.out.println("                    ...of which distributed: "
				+ (System.nanoTime() - gridstarttime) * 1e-9 + " seconds");
		}

		//System.out.println();

		// Result contains the a priori probability distribution matrix.
		return probmatrix;
	}

	static Sector findSector(int i, int j, Sector bottom) {
		// finds the sector in which the point i,j is located
		Sector current = bottom;
		// first find vertical position
		while (current.next != null) {
			if (i >= current.pos[0]
					&& (i < current.next.pos[0] || current.next.pos[0] == 0)) {
				break;
			}
			current = current.next;
		}
		// current now has the right column
		while (current.above != null) {
			if (j >= current.pos[1] - (i - current.pos[0])
					&& j < current.above.pos[1] - (i - current.pos[0])) {
				break;
			}
			current = current.above;
		}
		return current;

	}

	static int[] findPointST(int i, int j, Sector sector, int distance) {
		int result[] = new int[2];
		// lower half of the triangle, just return whatever sector is the one
		result[0] = sector.pos[0] + distance - i - 1;
		result[1] = j + i - (sector.pos[1] + sector.pos[0]);
		return result;
	}

	static int pairs(int a, int b) {
		int answer = 0;
		for (int c1 = 0; c1 < a; c1++) {
			for (int c2 = c1 + 1; c2 < b; c2++) {
				answer++;
			}
		}
		return answer;
	}
	
	static double log2(double val){
		return Math.log(val)/LOG_TWO;
	}

}