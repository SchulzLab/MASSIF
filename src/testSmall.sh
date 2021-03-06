#!/bin/bash

#TODO PWMs muss leer sein sonst geht das mit dem kleinen bsp nicht

	echo "test small example with ARNT and REST"

	#ensures that PFMs is really empty
	mkdir PFMs
	rm PFMs/*
	#call main script (domain information is used as prediction
	echo "using domain info as prediction"
 	time python mainScript.py  tests/transfac_testSmall.txt ../finalData/meme-5.0.2 tests/seq_testSmall/ testSmall tests/biologicalSignal_testSmall/

	echo "check if result is equal to provided results - If it print any line before ..... something went wrong"
	diff --suppress-common-lines result_fisherMethod_testSmall.txt  tests/result_fisherMethod_testSmall.txt
	echo "....."

	echo "using domain info as filter"
 	time python mainScript_FilteringFisherMethod.py  tests/transfac_testSmall.txt ../finalData/meme-5.0.2 tests/seq_testSmall/ testSmallFisher tests/biologicalSignal_testSmall/ RandomMotifs/pvalue_0.001_ThresholdDomainInfo.txt 

	echo "check if result is equal to provided results - If it print any line before ..... something went wrong"
	diff --suppress-common-lines result_fisherMethod_testSmallFisher.txt  tests/result_fisherMethod_testSmallFisher.txt
	echo "....."
