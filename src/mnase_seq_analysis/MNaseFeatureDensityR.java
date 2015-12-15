package get_mnase_density;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.SAMRecord; 
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader;

import jbfunctions.BAMInput;
import jbfunctions.DensityEst;

public class MNaseFeatureDensityR {

	public static void main(String[] args) throws IOException {

		// Read in the parameters
		String bam_file_name = args[0];
		String chr = args[1];
		int start_pos = Integer.parseInt(args[2]);
		int end_pos = Integer.parseInt(args[3]);
		int fragL = Integer.parseInt(args[4]);
		int fragH = Integer.parseInt(args[5]);
		double bw = Double.parseDouble(args[6]);	
		
		// Get the read scale factor
		double read_scale_factor = (10.0E6 / BAMInput.get_number_aligned_reads(bam_file_name, new String[]{chr})) *
								   (BAMInput.get_chr_length(bam_file_name, chr) / 12156677.0);
			
		// Open the BAM File
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), 
												   new File(bam_file_name + ".bai")
												  );
		
		// Get the SAM Record Iterator
		SAMRecordIterator bam_itr = bam_file.queryOverlapping(chr, start_pos - 250, end_pos);
		
		// Set up the array around the ACS position
		double[] storage = new double[end_pos - start_pos + 1];

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

		// Close the bam file
		bam_file.close();
		
		// Print out the storage array
		for(Double d : storage)
			System.out.println(d * read_scale_factor);
		
	}
		
}
