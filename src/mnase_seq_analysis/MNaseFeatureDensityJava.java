package mnase_seq_analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord; 
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader;

import genomics_functions.DensityEst;
import genomics_functions.TF;

public class MNaseFeatureDensityJava {

	public static void main(String[] args) throws IOException {

		// Read in the parameters
		String bam_file_name = args[0];
		String feature_file_name = args[1];
		String output_file_name = args[2];
		int win = Integer.parseInt(args[3]);
		int fragL = Integer.parseInt(args[4]);
		int fragH = Integer.parseInt(args[5]);
		double bw = Double.parseDouble(args[6]);	
		
		// Read in the feature file name
		ArrayList<TF> feature_list = TF.read_in_tf_list(feature_file_name);
			
		// Open the BAM File
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), 
												   new File(bam_file_name + ".bai")
												  );
		
		// Open the output file
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,pos,strand,");
		String sep = ",";
		for(int w = -win; w <= win; w++){
			if(w == win)
				sep = "\n";
			output.write(w + sep);
		}
		
		// Iterate through each TF
		for(TF feature : feature_list){
		
			// Get the genomic coordinates for the feature
			String chr = feature.getChr();
			int start_pos = feature.getPos() - win;
			int end_pos = feature.getPos() + win;
			
			// Set up the array around the ACS position
			double[] storage = new double[end_pos - start_pos + 1];
						
			// Get the SAM Record Iterator
			SAMRecordIterator bam_itr = bam_file.queryOverlapping(chr, start_pos - 250, end_pos);
						
			while(bam_itr.hasNext()){
	
				// Get each read (the SAMRecord object has functions to get at the SAM attributes)
				SAMRecord r = bam_itr.next();
	
				// Get the read midpoint
				int read_start = r.getAlignmentStart();
				int read_width = r.getInferredInsertSize();
				int read_mid = read_start + (read_width/2);
				
				// Check if read falls within the window coordinates
				if((read_width >= fragL) && (read_width <= fragH)){
				
					// Add the read to the array
					DensityEst.add_to_chr_array(storage, 1, 
												read_mid - start_pos, bw
											   );
				
				}
	
			}
	
			// Close the bam iterator
			bam_itr.close();

			// Write the output
			feature.write_reads_output(output, storage, 1);
			
		}
		
		// Close the bam file
		bam_file.close();
		
		// Close the output buffer
		output.close();
				
	}
		
}
