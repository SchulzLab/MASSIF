#!/bin/bash


path_to_meme=$1


	#ensures that PFMs is really empty
	mkdir PFMs
	rm PFMs/*

	echo "test bigger example with 21 seq sets and 21 motifs"

	echo "using domain info as prediction"
 	time python mainScript.py tests/transfac_testBig.txt  ${path_to_meme} tests/seq_testBig/ testBig tests/biologicalSignal_testBig/

	
	echo "using domain info as filter"
	rm PFMs/*
 	time python mainScript_FilteringFisherMethod.py tests/transfac_testBig.txt  ${path_to_meme} tests/seq_testBig/ testBigFilter tests/biologicalSignal_testBig/ RandomMotifs/pvalue_0.001_ThresholdDomainInfo.txt
