package com.ppfold.main;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

public class AlignmentReaderTest extends TestCase {

	public void setUp() {
	}

	public void testAlignment() throws Exception{
		
		List<String> lines = new ArrayList<String>();
		
		lines.add(">sequence1");
		lines.add("GCTGTGTGTGAAGGGTGCCCA");
		lines.add("AACGGTGCGCTAGCGAAC");
		lines.add("");
		lines.add("> sequence2");
		lines.add("GCTGTGTAAGGGGCTCAAGGGTGCCCA");
		lines.add("GCTGTGTAAAAAGCTCAAGGGTGCCCA");
		lines.add("GCTGTGTAATTTGCTCAAGGGTGCCCA");
		lines.add(">sequence3");
		lines.add("CCGTAAGGGGCTCAAGGGTGCCCA");
		lines.add("CCTGTGTAAAAAGCTCAAGGGTGCCCA");
		lines.add("CCTGTGTAATTTGCTCAAGGGTGCCCA");
		
		Alignment align = AlignmentReader.parse(lines); 
		
		
	}

}
