import sys
import os
import operator
import subprocess
import math

#Global variables
NUMBER_CLUSTERS = 35
NUMBER_THREADS = 5
CLUSTER_JASPAR = "clusterJASPAR/"
NUMBER_MOTIFS = 100000.0 #random motifs
#PATH_RANDOM_MOTIFS = "RandomMotifs/RandomMotifs_100000_AverageEntropy_0.6_2/"  

def writeOutput(output, result, TF):

	#write output
	output.write("TF:\t" + TF + '\n')
	output.write(".\tpredicted motif\tscore\n")
	output.write('\n')
	counter = 0
	for t in result:
		counter = counter +1
		output.write(str(counter) + '\t' + str(t[0]) + '\t' + str(t[1]) + '\n')
	output.write("\n")

def evaluationPASTAA(PASTAA, name):

	output = open("result_PASTAA_" + name + ".txt", "w")
	store_result = {}
	#iterates over the results of the different seq sets	
	for PASTAA_dir in os.listdir(PASTAA):
		pos = PASTAA_dir.find("_")
	    	if PASTAA_dir[:pos] == "enrichment": 
			actual_result = open(PASTAA + "/"  + PASTAA_dir , 'r')
			result = []
			for line in actual_result:
				line = line.split('\t')
				name = line[0].strip()
				value = min(float(line[1]), 1.0)
				result.append((name,value))

			name_TF = PASTAA_dir[pos+1:-4]
			#write output
			writeOutput(output, result, name_TF)
			store_result[name_TF] = result
			actual_result.close()
	output.close()
	return(store_result)

def evaluationCentriMo(centrimo_output_dir, name):
	
	output = open("result_CentriMo_" + name + ".txt", "w")
	store_result = {}
	#iterates over the results of the different seq sets	
	for centrimo_dir in os.listdir(centrimo_output_dir):
		if centrimo_dir == "PWMsInMemeFormat.txt":
			continue
		else:
			for filename in os.listdir(centrimo_output_dir + "/" +centrimo_dir):
	    			if filename == "centrimo.tsv": 
					first_line = True
					actual_result = open(centrimo_output_dir  + "/" + centrimo_dir +'/'+ filename, 'r')
					result = []
					for line in actual_result:
						if line[0] == '#' or first_line == True or line.strip() == "":
							first_line = False
							continue
						else:	
							line = line.strip()
							line = line.split('\t')
							name = line[2].split(' ')
							number = float(line[5]) #adjusted pvalue
							#number = float(line[4]) #E value
							result.append((name[-1], number))

					result = sorted(result, key= lambda element:(element[1]), reverse = False)
			
					pos = centrimo_dir.find("_")
					name_TF = centrimo_dir[pos+1:]
					#write output
					writeOutput(output, result, name_TF)
					store_result[name_TF] = result
					actual_result.close()
	output.close()
	return (store_result)

def call(command):

	try:
		output = subprocess.check_output(command, shell = True)
		print(output)
	except subprocess.CalledProcessError as e:
		print(e.output)
	#	 e.output.startswith('error: {'):
    	#		error = json.loads(e.output[7:]) # Skip "error: "
    	#		print error['code']
    	#		print error['message']

def lookupPvalue(actual_DBD, pvalue, values):

	counter = 1
	pos =-1.0
	for i in values:
		i = float(i)
		if i > pvalue:
			counter = counter + 1 
		if i == pvalue:
			pos = counter
			break
		if i < pvalue:			
			pos = counter -1
			break
	if pos == 0:
		return 1/NUMBER_MOTIFS
	elif pos == -1.0:
		return 1
	else:
		return pos/NUMBER_MOTIFS

def determinePvaluesOfDomainInfo(inputFile, DBD_info_dir,DBD_file, outputFile):
	output = open(outputFile, "w")
	input_ = open(inputFile, "r")

	dis = {}

	# stores all DomainInfovalues from the randomMotifs sorted in inreascing order
	for i in range(30):
		DBD_number = i + 1
		f = open ( DBD_info_dir + "/DomainInfo_" + str(DBD_number) + "_" + str(DBD_number) + ".txt","r")
		DBD = []
		first_line =  True
		for line in f:
			if not first_line: 
				line = line.strip()
				line = line.split("\t")
				DBD.append(line[3])
			else: 
				first_line = False
	
		DBD.sort(reverse = True)
		dis["DBD_" + str(DBD_number)] = DBD
	#store for each TF the DBD	
	TFtoDBD = open(DBD_file, "r") 
	lookupDBD = {}
	for line in TFtoDBD:
		line = line.strip()
		line = line.split("\t")
		lookupDBD[line[0]] = line[1]	


	results = {}
	each_set = []
	TF = ""
	for line in input_:
		line = line.strip()
		if not line:
			output.write("\n")
		else:
			line = line.split("\t")
			if line[0] == "TF:":
				if each_set != []:
					results[TF] = each_set
					each_set = []
				TF = line[1]
				actual_DBD = lookupDBD[TF]
				helper = dis[actual_DBD]
				output.write(line[0] + "\t"  + line[1] + "\n") 
			elif line[0] == ".":
				output.write(line[0] + "\t"  + line[1] + "\t" + line[2] +  "\n") 
			else:
				pvalue = lookupPvalue(actual_DBD, float(line[2]),helper)
				each_set.append((line[1], pvalue))
				output.write(line[0] + "\t"  + line[1] + "\t" + str(pvalue) +  "\n") 

	results[TF] = each_set
	each_set = []
	return results
	

def determineSignificantMotifs(result_file, pvalue_cutoff_file, TF_to_DBD_file, PWM_dir, transfac_file, name ):

	result = open(result_file, "r")
	pvalue_cutoff = open(pvalue_cutoff_file, "r")
	TF_to_DBD = open(TF_to_DBD_file, "r")
	transfac = open(transfac_file, "r")
	output_ = open("significantMotifs_" + name + ".txt", "w")

	#read transfac_file
	pwms = {}
	pwm = ""
	end = False
	for line in transfac:
		if line[0] == "/" and end == True:
			end = False
			pwms[motif] = pwm
			pwm = ""
			motif = ""	
		if line[0] =="I":
			helper = line.strip()
			helper = helper.split("\t")
			motif = helper[1]
			end = True
		pwm = pwm + line
	pwms[motif] = pwm

	#read TF_to_DBD
	TFtoDBD = {}
	for line in TF_to_DBD: 
		line = line.strip()
		line = line.split("\t")
		TFtoDBD[line[0]] = line[1]
	
	#read Pvalue cutoff
	pvalueCutoff = {}
	for line in pvalue_cutoff:
		line = line.strip()
		line = line.split("\t")
		pvalueCutoff[line[0]] = line[1] 

	#read result
	first = True
	for line in result:
		line = line.strip()
		if line == "" or line[0] == ".":
			continue
		else:
			line = line.split("\t")
			if line[0] == "TF:":
				TF = line[1]
				output_.write("............")
				output_.write("TF:\t" + TF + "\n")
				if not first:
					output.close()
				else:
					first = False
				output = open(PWM_dir + "transfac_" + TF + ".txt", "w")
			else:
				motif = line[1]
				DBD = TFtoDBD[motif]
				cutoff = pvalueCutoff[DBD]
				value = line[2]
				if value >= cutoff:
					output.write(pwms[motif])
					output_.write(motif + "\n")	

#---------------
# specification of the input files:
#	transfac_file
#	path_to_centrimo_dir meme_?
#	path_to_dir containing the fasta_seq
#	name of the run
#	path_to_biological_signal (for PASTAA)
#---------------
def main(transfac_file, CentriMo, fasta_dir, name , biological_signal, pvalue_cutoff_file):

	fasta_files = [f for f in os.listdir(fasta_dir)] #if isfile(join(fasta_dir, f))]

	#---------------
	#TFtoDoamin 
	#---------------
	print("----------\nTFtoMOtifAnalysis\n----------")

	#---------------
	#Step 1: create PFMs 
	#---------------
	command = "mkdir PFMs"
	call(command)
	
	command = "python src/convertCountMatrixToFreqMatrix.py "  + transfac_file +  " PFMs/" 
	call(command)
	
	#---------------
	#Step2: call TFtoMOtifDomainInfo
	#---------------
	command = "./src/TFtoMotifDomainInfo PFMs/ " + fasta_dir + " result_DomainInfo_" + name + ".txt " + CLUSTER_JASPAR + "sum_TFs.txt " + CLUSTER_JASPAR + "DBDs.txt " + CLUSTER_JASPAR + "TF_to_DBD.txt"
	call(command)
	#---------------
	#Step3: determine pvalues for resulting DomainInfo values
	#---------------
	print("----------\ndeterminePvaluesOfDomainInfo\n----------")
	
	command = "mkdir setOfPWMs_" + name  + "/"
	call(command)
	determineSignificantMotifs("result_DomainInfo_" + name + ".txt", pvalue_cutoff_file, CLUSTER_JASPAR + "TF_to_DBD.txt", "setOfPWMs_" + name + "/", transfac_file, name)

#	result_DomainInfo = determinePvaluesOfDomainInfo("result_DomainInfo_" + name + ".txt", PATH_RANDOM_MOTIFS, CLUSTER_JASPAR + "TF_to_DBD.txt", "result_DomainInfoPvalues_" + name + ".txt")


	#---------------
	#CentriMo
	#---------------
	print("----------\nCentriMo\n----------")

	command = "mkdir CentriMo_" + name
	call(command)
	
	#TODO: hier koennte man denke auch parallelisieren
	for f in fasta_files:
		TF = f[:-3]
		#step1: parse transfac PWMs in meme format
		command = CentriMo + "/scripts/transfac2meme setOfPWMs_" + name + "/transfac_" + TF + ".txt" +  " >CentriMo_" + name + "/PWMsInMemeFormat.txt"
		call(command)
		#step2: call CentriMo
		command =  CentriMo + "/src/centrimo --oc " +"CentriMo_"+ name + "/centrimo_"+ TF + " --ethresh 1000000000 " + fasta_dir + f + " CentriMo_" + name + "/PWMsInMemeFormat.txt"
		call(command)
			
	result_centrimo = evaluationCentriMo("CentriMo_" + name, name)
	
	#---------------
	#PASTAA 
	#---------------
	print("----------\nPASTAA\n----------")


	command = "mkdir  PASTAA_" + name
	call(command)

	for f in fasta_files:
		TF = f[:-3]
		command = "./src/PSCM_to_PSEM setOfPWMs_" + name + "/transfac_" + TF + ".txt" +  " >PASTAA_" + name  + "/energy.txt"
		call(command)
		command =  "./src/TRAP PASTAA_"+ name + "/energy.txt " +  fasta_dir + f + " > PASTAA_"+ name + "/affinity_" + TF + ".txt"	
		call(command)

		command =  "./src/PASTAA PASTAA_" + name  + "/affinity_" + TF + ".txt " + biological_signal + "gene_list_" + TF +  ".txt |  sort -k2,2 -g  > PASTAA_"+ name + "/enrichment_" + TF + ".txt"
		print (command)
		call(command)

	result_pastaa = evaluationPASTAA("PASTAA_"+ name, name)

	#---------------
	#fisher test
	#---------------
	print("----------\nfishers test\n----------")

	#---------------
	#Step1: normalize PASTAA using Bonferroni
	#---------------
	output_pastaa_normalized = open("result_PASTAA_normalized_"+ name + ".txt", "w")
	keys = result_pastaa.keys()
	for k in keys:
		helper = result_pastaa[k]
		for e in range(len(helper)):
			helper[e]= (helper[e][0], min(float(helper[e][1]) * 2413, 1.0)) # number of test that pastaa performs
		result_pastaa[k] = helper
		writeOutput(output_pastaa_normalized, result_pastaa[k], k)	
	output_pastaa_normalized.close()	

	#---------------
	#Step2: fisher analysis
	#---------------

	result_fisher = {}
	for k in keys:
		result_p = sorted(result_pastaa[k] , key=lambda tup: tup[0])
		result_c = sorted(result_centrimo[k], key = lambda tup :tup[0])
		result_f = []
		counter_centrimo = -1
		for i in range(len(result_p)):
			counter_centrimo = counter_centrimo + 1
			if (result_p[i][0] != result_c[counter_centrimo][0]): # or result_c[counter_centrimo][0] != result_d[i][0]):
				motif = result_p[i][0]
				current_value = -2 * (math.log(max(float(result_p[i][1]),2.2250738585072014e-308 )))	
				result_f.append((motif, current_value))
				counter_centrimo = counter_centrimo -1
			else:
				motif = result_p[i][0]
				#fishers method
				current_value = -2 * (math.log(max(float(result_p[i][1]),2.2250738585072014e-308 )) + math.log(max(float(result_c[counter_centrimo][1]), 2.2250738585072014e-308)))	
				result_f.append((motif, current_value))
		
		result_fisher[k] = sorted(result_f, key=lambda tup : tup[1], reverse = True)

	final_output = open("result_fisherMethod_" + name + ".txt", "w")
	for k in keys:
		writeOutput(final_output, result_fisher[k], k)	

# call main
if(len(sys.argv))< 7:
	print("Usage python ./mainScript_domainInfoFilter.py transfac_file, CentriMo_source_code_dir, fasta_dir,name,  biological_signal, pvalue_cutoff")
else:
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
