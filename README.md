# MASSIF - motif association with domain information
TODO: Overview and short explanation of MASSIF

# Install MASSIF

Necessary installed software and packages:

- [The Meme Suite](http://meme-suite.org/doc/download.html) (version 5.0.2) to use CentriMo
- compiler that is able to use openMP (and omp.h file which is part of the GNU OpenMP Library)

Download the repository.
To compile the C++ code of MASSIF  perform the following commands: 
```
cd MASSIF/src
make
```
# Run test cases
To check if MASSIF is working correctly, we provide two test cases, a small example with only two sequence sets  and X motifs as well as a bigger one with 19 sequence sets and the corresponding 19 motifs.
To start the test cases run
``` 
cd ../
bash testSmall.sh path_to_meme_suite
bash testBig.sg path_to_meme_suite
```
where path_to_meme_suite is the path to the meme suite (something like /Home/.../meme-2.0.5/).


# Required input of MASSIF

# Output of MASSIF
