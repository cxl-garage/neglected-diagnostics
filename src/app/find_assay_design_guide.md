### QuickGuide

1. Drag and drop, or Click on the [Browse files] button to choose and upload the individual sequence files or non-aligned multiple sequence files (FASTA or FAS) for the Assay Target(s). *These are the sequences that the assay should be able to detect.*    
2. Drag and drop, or Click on the [Browse files] button to choose and upload the individual sequence files or non-aligned multiple sequence files (FASTA or FAS) for the Assay Off-Target(s). *These are the sequences that the assay should not be able to detect.*  
3. Input Assay Design Area Parameters  
    a. From the set of sequences included in the Target sequences, choose one to be the reference sequence. Input the sequence header (not the name of the file) under [Select Reference Sequence].  
    b. Enter a number into the [assay design region size] box, or use the +/- symbols to adjust the default number. *This should be equal to or larger than the intended fragment size for the assay. It is the total area where targets and off targets will be compared.*  
    c. Enter a number into the [assay design region slide size] box, or use the +/- symbols to adjust the default number. This designates the shift in base pair distance between each design area along the total sequence length which is defined by the reference sequence.  
    d. Enter a number into the [maximum differences allowed between assay area and target sequences] box, or use the +/- symbols to adjust the default number. *This controls the calculated average differences between the target sequences and the assay design area candidates.* 
    e. Enter a number into the [maximum differences allowed between assay area and off-target sequences] box, or use the +/- symbols to adjust the default number. *This controls the calculated average differences between the off-target sequences and the assay design area candidates.*
4. Click on the [Find Assay Design Area] button to run the process. *Note that large files, small target region sizes or a combination of both will require more processing time.
5. Select the download option from the dropdown.  
    a. Select [Download as CSV] to retrieve the result table displayed above  
    b. Select [Download sequences as a fasta file] to retrieve the assay design area sequences only.  
6. Click on the [Download] button. 