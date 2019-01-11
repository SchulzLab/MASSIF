# MASSIF - motif association with domain information
TODO: Overview and short explanation of MASSIF

# Install MASSIF

Necessary installed software and packages:

- [The Meme Suite](http://meme-suite.org/doc/download.html) (version 5.0.2) to use CentriMo
- compiler that is able to use openMP (and omp.h file which is part of the GNU OpenMP Library)

Download the repository. Notice, that the code for PASTAA is enclosed. 
To compile the C++ code of MASSIF and PASTAA perform the following commands: 
```
cd MASSIF/src
make
```
# Run test cases
To check if MASSIF is working correctly, we provide two test cases, a small example with only two sequence sets  and 11 motifs as well as a bigger one with 19 sequence sets and the corresponding 19 motifs.
To start the test cases run
``` 
cd ../
bash testSmall.sh path_to_meme_suite
bash testBig.sg path_to_meme_suite
```
where *path_to_meme_suite* is the path to the meme suite (something like /Home/.../meme-2.0.5/).

# Required input of MASSIF

** Using the domain information as prediction **
 To run the script where MASSIF apply the domain information as prediction the following input is required:
 
 - *motif_file* file that contains all consider motifs as transfer format (see tests/transfac_testSmall.txt for an example).
 - *path_to_meme_suite* path to the meme directory  
 - *path_to_seq_dir* path to a directory that contain for each considered TF a fasta file (see tests/seq_testSmall/)
 - *name* name for the output files
 - *path_to_biological_signal* path to a directory that contains for each TF a file with the biological signal (see tests/biologicalSignal_testSmall/).  TODO: explain how this files must look like
 
** Using the domain information as a filter **
The script that uses the domain information as a filer needs additionally:
- *thresholds_domainInfo* a file that gives for a specific pvalue threshold the corresponding domain information pro DNA-binding domain. In the directory RandomMotifs there are several possible files for different pvalue thresholds provided. For the results shown in our paper we used as a pvalue threshold 0.001 (meaning, we used the file pvalue_0.001_ThresholdDomainInfo.txt)

# Output of MASSIF







