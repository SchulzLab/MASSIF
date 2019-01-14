# MASSIF - motif association with domain information
TODO: Overview and short explanation of MASSIF

# Installation

Necessarily installed software and packages:

- [The Meme Suite](http://meme-suite.org/doc/download.html) (version 5.0.2 or 5.0.3) to use CentriMo
- C++ compiler that is able to use openMP (and omp.h file which is part of the GNU OpenMP Library)
(Notice, for Mac OS the clang++ compiler and the packages libopm is necessary)

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

# Required input

**Using the domain information as prediction**

 To run the script where MASSIF apply the domain information as prediction the following input is required:
 
 - **motif_file** file that contains all consider motifs as TRANSFAC format (see [tests/transfac_testSmall.txt](tests/transfac_testSmall.txt for an example).
 - **path_to_meme_suite** path to the meme directory  
 - **path_to_seq_dir** path to a directory that contains for each considered TF a fasta file (see [tests/seq_testSmall/](tests/seq_testSmall/)
 - **name** name for the output files
 - **path_to_biological_signal** path to a directory that contains for each TF a file with the biological signal (see [tests/biologicalSignal_testSmall/](tests/biologicalSignal_testSmall/).  TODO: explain how this files must look like
 
 For instance, to run the small example the following command is required:
 ```
 python mainScript_DomainInfoPrediction.py  tests/transfac_testSmall.txt path_to_meme_suite/ tests/seq_testSmall/ testSmall tests/biologicalSignal_testSmall/
 ```
 where *path_to_meme_suite* is the path to the meme suite (something like /Home/.../meme-2.0.5/).
 
**Using the domain information as a filter**

The script that uses the domain information as a filer needs additionally:
- **thresholds_domainInfo** a file that gives for a specific pvalue threshold the corresponding domain information pro DNA-binding domain. The directory [RandomMotifs](RandomMotifs/) provides  several possible files for different pvalue thresholds. For the results shown in our paper, we used as a pvalue threshold 0.001 (meaning, we used the file [pvalue_0.001_ThresholdDomainInfo.txt](RandomMotifs/pvalue_0.001_ThresholdDomainInfo.txt))

To run the small example we need the following command:
 ```
  python mainScript_DomainInfoFilter.py  tests/transfac_testSmall.txt path_to_meme_suite/ tests/seq_testSmall/ testSmallFisher tests/biologicalSignal_testSmall/ RandomMotifs/pvalue_0.001_ThresholdDomainInfo.txt
 ```
 with *path_to_meme_suite* is the path to the meme suite (something like /Home/.../meme-2.0.5/).

# Output 
Both variations of MASSIF produce the following output: 

- **CentriMo_name** original result from CentriMo.
- **result_CentriMo_name.txt** parsed result from CentriMo.
- **PASTAA_name** original result from PASTAA.
- **result_PASTAA_name.txt** parsed result from PASTAA.
- **result_DomainInfo_name.txt** ranking of the domain information (the bigger the value the better).
- **result_fisherMethod_name.txt** final result. For each TF a ranking of motifs is given. The higher the ranking of the motif the more likely this motif corresponds to the TF. 

Using the domain information as prediction outputs additionally:

- **result_DomainInfoPvalues_name.txt** ranking  of the domain information interpreted as pvalue.

Applying the domain information as filter also provides as output:

- **significantMotifs_name.txt** file that contains for each TF the reduced motif set after the domain information is applied as filter.
- **setOfPWMs_name** directory that contains for each TF the motif_file for the reduced motif set.
