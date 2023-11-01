# Welcome to Seq2Dx!
The following tools can be used independently or as a pipeline for analyzing sequence data for diagnostic targeting. 

### Sequence to Diagnostic Pipeline
Step 1: 
Generate sequence data for target and off-target species.  
***Sequence Search*** can be used to search for publicly available sequences.  
Step 1B: 
Optional - Clean up sequence data and characterize sequence variability.  
Use a sequence visualization and alignment software to screen the sequences for errors or low quality.  
***Sequence Haplotypes*** can be used to collapse sequences into groups based on patterns of variability.  
Step 2: 
Compare target and off-target sequences to identify top candidates for designing primers.  
***Find Assay Target Area*** can identify sections of the sequence data where target sequences are most similar to each other, but most different from off-target sequences.  
Step 2B: 
Optional - Characterize inter and intra-species or strain variation to inform primer site locations.  
***Sequence Variability*** can be used to generate a table that lists all positions where variation exists across sequences.  
Step 3: 
Design primers using the target area candidates.  

### Sequence Input
This tool integrates with [GenBank](https://www.ncbi.nlm.nih.gov/nucleotide/) and can search for sequences from the nucleotide database (more databases coming soon). Search terms can be general (e.g. cougar mitochondrial) or can be more specific using the NCBI [Search field descriptions](https://www.ncbi.nlm.nih.gov/books/NBK49540/).  
This tool goes one step further by allowing the search results to be filtered in multiple layers, unlocking the ability to refine and clean up the sequence list before downloading. There is the option of downloading the search result table or the sequence fasta file.  

### Sequence Haplotypes
This tool analyzes aligned multiple sequence fasta files and collapses them into groups of variability. In cases where sequences have varying length, a base or reference sequence can be used to fill in sequence gaps and reduce the number of outliers that have to be collapsed manually.  

### Sequence Variability
This tool analyzes aligned multiple sequence fasta files to generate a variability table of base pair composition along all positions. A threshold can be set to only show positions where variability is below a certain percentage - for example a threshold of 99 would exclude all positions where there is a 100% consensus across all sequences. Gaps are not scored as differences and are ignored by the algorithm.   

### Find Assay Target Area
This tool analyzes a set of sequences labeled as Target (sequences that the diagnostic assay should detect) against a set of sequences labeled as Off-Target (sequences that the diagnostic assay should NOT detect) to identify areas where targets are conserved but dissimilar to the off-targets. This comparison is done along the entire length of the reference target sequence by splitting up the pairwise comparison into smaller windows that slide along the sequences. Further details on this algorithm can be found [here](https://github.com/kmceres/Thylacine_Design/tree/main/general_primer_design).  


Powered by Conservation X Labs and the University of Washington Scientific Software Engineering Center at the eSciences Institute.
