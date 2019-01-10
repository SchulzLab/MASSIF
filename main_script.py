import sys
import os
import operator
import subprocess
import math

#Global variables
CLUSTER_INFO_FILE = "clustering/info.txt"
NUMBER_CLUSTERS = 35
NUMBER_THREADS = 5
CLUSTER_JASPAR = "clusterJASPAR/"
NUMBER_MOTIFS = 100000.0 #random motifs
PATH_RANDOM_MOTIFS = "RandomMotifs/RandomMotifs_100000_AverageEntropy_0.6_2/"  
#reads from file which TF belongs to which cluster
def readClusterInfo():

	cluster_info = open(CLUSTER_INFO_FILE, 'r')
	cluster = {}
	for line in cluster_info:
		line = line.strip()
		line = line.split('\t')
		for elem in line:
			if elem[0] == 'c':
				c = elem
			else:
	#			pos = elem.find('-')
	#			elem = elem[pos+1:]
				cluster[elem] = c	
	return cluster

def evaluation(result, TF,counter,  cluster, wrong_predicted):
	cluster_list = []
	index = -1
	number_c = [0] * NUMBER_CLUSTERS
	correct_predicted = False
	counter_ = 1
	actual_value = 0
	end = False
	for t in result:
		c = cluster[t[0].strip()] #looks up cluster of the current motif
		if not c in cluster_list:
			index = index + 1 #counter_
			cluster_list.append(c)
		if cluster[TF.strip()] == c and index < 40:
			#if not c in cluster_list:
			for i in range(index, NUMBER_CLUSTERS):
				number_c[i] = number_c[i] + 1
			if index < 1:
				correct_predicted = True
			end = True  
		if end ==True:
			break
	if correct_predicted == False:
		wrong_predicted.append(TF)
	for i in range(len(counter)):
		counter[i] = counter[i] + number_c[i]
	return counter

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
		
def writeOutputEval(output, number_files, counter_eval, wrong_predicted):
	output.write("with cluster-info\tTop 1:\t" + "Top 2:\t" +"Top 3:\t" + "Top 4:\t" + "Top 5:\t" + "Top 7:\tTop 10" + "\twrong predicted Top 1\n")
	output.write("result total number:" +"\t" + str(counter_eval[0]) + " , " + str(counter_eval[1]) + " , " + str(counter_eval[2]) + " , " + str(counter_eval[3]) + " , " +  str(counter_eval[4]) + " , " + str(counter_eval[6]) + " , " + str(counter_eval[9])+ "\t" +str(wrong_predicted) +"\n")	
	output.write("result in %:\t" + str(float(counter_eval[0]/ number_files)) + " , " + str(float(counter_eval[1]/number_files)) + " , " + str(float(counter_eval[2]/number_files)) + " , " + str(float(counter_eval[3]/number_files)) + " , " +  str(float(counter_eval[4]/number_files)) + " , " + str(float(counter_eval[6]/number_files)) + " , " + str(float(counter_eval[9]/number_files))+ "\n")	

def evaluationPASTAA(PASTAA, name, cluster):

	output = open("result_PASTAA_" + name + ".txt", "w")
	output2 = open("evaluaton_PASTAA_" + name + ".txt", "w")	
	counter_eval = [0] *NUMBER_CLUSTERS
	wrong_predicted = []
	number_files = 0.0
	store_result = {}
	#iterates over the results of the different seq sets	
	for PASTAA_dir in os.listdir(PASTAA):
		pos = PASTAA_dir.find("_")
	    	if PASTAA_dir[:pos] == "enrichment": 
			number_files = number_files + 1
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
			#evaluate result
			counter_eval = evaluation(result, name_TF, counter_eval, cluster, wrong_predicted)
			actual_result.close()
	output.close()
	writeOutputEval(output2, number_files, counter_eval, wrong_predicted)
	return(store_result)

def evaluationCentriMo(centrimo_output_dir, name, cluster):
	
	output = open("result_CentriMo_" + name + ".txt", "w")
	output2 = open("evaluaton_CentriMo_" + name + ".txt", "w")	
	counter_eval = [0] *NUMBER_CLUSTERS
	wrong_predicted = []
	number_files = 0.0
	store_result = {}
	#iterates over the results of the different seq sets	
	for centrimo_dir in os.listdir(centrimo_output_dir):
		if centrimo_dir == "PWMsInMemeFormat.txt":
			continue
		else:
			number_files = number_files + 1
			#NOTICE: hier muss iwo das / dabei sein
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
	
					#evaluate result
					counter_eval = evaluation(result, name_TF, counter_eval, cluster, wrong_predicted)
					actual_result.close()
	output.close()
	writeOutputEval(output2, number_files, counter_eval, wrong_predicted)
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

#	score = 100000
	counter = 1
	pos =-1.0
#	print(pvalue)
	for i in values:
		i = float(i)
		if i > pvalue:
#			print("smaller")
			counter = counter + 1 
			#score = i
		if i == pvalue:
#			print("equal")
			pos = counter
			break
		if i < pvalue:			
#			print("i smaller")
			pos = counter -1
			break
#	print(pos)
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
	

#---------------
# specification of the input files:
#	transfac_file
#	path_to_centrimo_dir meme_?
#	path_to_dir containing the fasta_seq
#	name of the run
#	path_to_biological_signal (for PASTAA)
#---------------
def main(transfac_file, CentriMo, fasta_dir, name , biological_signal):


	fasta_files = [f for f in os.listdir(fasta_dir)] #if isfile(join(fasta_dir, f))]
	#---------------
	#CentriMo
	#---------------

	#step1: call CentriMo
	command = "mkdir CentriMo_" + name
	call(command)
	#step2: parse transfac PWMs in meme format
	command = CentriMo + "/scripts/transfac2meme " + transfac_file + " >CentriMo_" + name+ "/PWMsInMemeFormat.txt"
	print(command)
	call(command)
	
	#TODO: hier koennte man denke auch parallelisieren
	for f in fasta_files:
		TF = f[:-3]
		print(name)
		command =  CentriMo + "/src/centrimo --oc " +"CentriMo_"+ name + "/centrimo_"+ TF + " --ethresh 1000000000 " + fasta_dir + f + " CentriMo_"+  name +"/PWMsInMemeFormat.txt"
		call(command)
			
	cluster = readClusterInfo()
	result_centrimo = evaluationCentriMo("CentriMo_" + name, name, cluster)
	#print ("result CentriMo : " + str(result_centrimo)) 
	
	#---------------
	#PASTAA 
	#---------------
	command = "mkdir  PASTAA_" + name
	call(command)
	command = "./src/PSCM_to_PSEM " + transfac_file +  " >PASTAA_" + name  + "/energy.txt"
	print(command)
	call(command)

	for f in fasta_files:
		TF = f[:-3]
		#print("TRAP")
		command =  "./src/TRAP PASTAA_"+ name + "/energy.txt " +  fasta_dir + f + " > PASTAA_"+ name + "/affinity_" + TF + ".txt"	
		print(command)
		call(command)

		#print("PASTAA")
		#TODO: biological_signal files muessen iwie anders benannt werden -> so ist es schwierig fuer user es zu benutzen
	#	command =  "./src/PASTAA PASTAA_" + name  + "/affinity_" + TF + ".txt " + biological_signal + "gene_id_with_signal_" + TF +  ".txt |  sort -k2,2 -g  > PASTAA_"+ name + "/enrichment_" + TF + ".txt"
		command =  "./src/PASTAA PASTAA_" + name  + "/affinity_" + TF + ".txt " + biological_signal + "gene_list_" + TF +  ".txt |  sort -k2,2 -g  > PASTAA_"+ name + "/enrichment_" + TF + ".txt"
		print (command)
		call(command)

	result_pastaa = evaluationPASTAA("PASTAA_"+ name, name, cluster)
	#print ("result PASTAA : " + str(result_pastaa)) 

	#---------------
	#TFtoDoamin 
	#---------------
	#TODO: hier dann alles hin um TFtoDomainInfo aufzurufen, nochmal genau nachsehen was man alles braucht um das file aufzurufen	
	# TODO: noch rausfinden wie das mit RandomMotifs war. Wo das gespeichert wurde und vorallem auf was wir uns da geeinigt hatten wie der pvalue sein soll
	# muss das ein anwender ueberhaupt benutzen oder geben wir da das file vor?!
	
	print("TFtoMotifAnalysis")
	#---------------
	#Step 1: create PFMs 
	#---------------
	command = "mkdir PFMs"
	call(command)
	
	command = "python convertCountMatrixToFreqMatrix.py "  + transfac_file +  " PFMs/" 
	print(command)
	call(command)
	
	#---------------
	#Step2: call TFtoMOtifDomainInfo
	#---------------
	command = "./src/TFtoMotifDomainInfo PFMs/ " + fasta_dir + " result_DomainInfo_" + name + ".txt " + CLUSTER_JASPAR + "sum_TFs.txt " + CLUSTER_JASPAR + "DBDs.txt " + CLUSTER_JASPAR + "TF_to_DBD.txt"
	print(command)
	call(command)
	#---------------
	#Step3: determine pvalues for resulting DomainInfo values
	#---------------
	print("determinePvaluesOfDomainInfo")
	result_DomainInfo = determinePvaluesOfDomainInfo("result_DomainInfo_" + name + ".txt", PATH_RANDOM_MOTIFS, CLUSTER_JASPAR + "TF_to_DBD.txt", "result_DomainInfoPvalues_" + name + ".txt")

	#---------------
	#fisher test
	#---------------
	print("fishers test")

	#---------------
	#Step1: normalize PASTAA using Bonferroni
	#---------------
	output_pastaa_normalized = open("result_PASTAA_normalized_"+ name + ".txt", "w")
	keys = result_pastaa.keys()
	for k in keys:
		helper = result_pastaa[k]
		for e in range(len(helper)):
			helper[e]= (helper[e][0], min(float(helper[e][1]) * 2413, 1)) # number of test that pastaa performs
		result_pastaa[k] = helper
		writeOutput(output_pastaa_normalized, result_pastaa[k], k)	
	output_pastaa_normalized.close()	
	
	#print("check normalized pastaa result: " + str(result_pastaa))

	#---------------
	#Step2: fisher analysis
	#---------------

	result_fisher = {}
	for k in keys:
		result_p = sorted(result_pastaa[k] , key=lambda tup: tup[0])
		result_c = sorted(result_centrimo[k], key = lambda tup :tup[0])
		result_d = sorted(result_DomainInfo[k], key = lambda tup : tup[0])
		result_f = []
		counter_centrimo = -1
		for i in range(len(result_p)):
			counter_centrimo = counter_centrimo + 1
			if (result_p[i][0] != result_c[counter_centrimo][0] or result_c[counter_centrimo][0] != result_d[i][0]):
				if result_p[i][0] != result_c[counter_centrimo][0]:
					motif = result_p[i][0]
					current_value = -2 * (math.log(max(float(result_p[i][1]),2.2250738585072014e-308 )) + math.log(max(float(result_d[i][1]), 2.2250738585072014e-308 )))	
					#current_value = -2 * (math.log(max(float(result_p[i][1]),2.2250738585072014e-308 )))	
					result_f.append((motif, current_value))
					counter_centrimo = counter_centrimo -1
				else:
					print("NOPE\t" + result_p[i][0] + " " + result_c[counter_centrimo][0] + " " + result_d[i][0] + "\n" )
					break
			else:
				motif = result_p[i][0]
				#fishers method
				current_value = -2 * (math.log(max(float(result_p[i][1]),2.2250738585072014e-308 )) + math.log(max(float(result_c[counter_centrimo][1]), 2.2250738585072014e-308)) + math.log(max(float(result_d[i][1]), 2.2250738585072014e-308 )))	
				#current_value = -2 * (math.log(max(float(result_p[i][1]),2.2250738585072014e-308 )) + math.log(max(float(result_c[counter_centrimo][1]), 2.2250738585072014e-308)))	
				result_f.append((motif, current_value))
		
		result_fisher[k] = sorted(result_f, key=lambda tup : tup[1], reverse = True)

	#print("result_fisher: " + str(result_fisher))
	final_output = open("result_fisherMethod_" + name + ".txt", "w")
	output2 = open("evaluation_fisherMethod_" + name + ".txt", "w")	
	counter_fisher = [0] *NUMBER_CLUSTERS
	wrong_predicted = []
	number_files = 0.0
	for k in keys:
		number_files = number_files + 1
		writeOutput(final_output, result_fisher[k], k)	
		#evaluation
		counter_fisher = evaluation(result_fisher[k], k, counter_fisher, cluster, wrong_predicted)
	print(counter_fisher)	
	writeOutputEval(output2, number_files, counter_fisher, wrong_predicted)
	

# call main
if(len(sys.argv))< 6:
	print("Usage python ./main_script.py transfac_file, CentriMo_source_code_dir, fasta_dir,name,  biological_signal")
else:
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

#---------------
#transfer transfac matrices to meme requried meme format (rewriten because CentriMo uses a .pl script for that)
#---------------
#def transfac2meme(transfac_file):

#	transfac_PWMs = open(transfac_file, "r")
#	output = open("PWMsInMemeFormat.txt", "w")
	#write heading
#	output.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n")

#	name = ""
#	motif_A =  []
#	motif_C =  []
#	motif_G =  []
#	motif_T =  []

#	first_time = True
#	for line in transfac_PWMs:
#		line = line.strip()
#		line = line.split("\t")
#		if line[0] == "//": #new motif starts
#			if first_time == True:
#				first_time = False
#			else:			
#				output.write("MOTIF\t" + name + "\n\n")
#				output.write("letter-probability matrix: alength= 4 w= " + str(len(motif_A))+ " nsites= " + str(sum(motif_A))+ " E= 0\n")
#				for i in range(len(motif_A)):
#					sum_ = motif_A[i] + motif_C[i] + motif_G[i] + motif_T[i]
#					output.write(str(round(motif_A[i]/sum_, 6)) + "\t" + str( round(motif_C[i]/sum_, 6)) + "\t" + str(round(motif_G[i]/sum_, 6)) + "\t" + str(round(motif_T[i]/sum_, 6)) + "\n")
#
#			output.write("\n")
#			name = ""
#			motif_A = []
#			motif_C = []
#			motif_G = []
#			motif_T = []
#
#		elif line[0] == "XX" or line[0] == "P0":
#			continue
#
#		elif line[0][0] == "I":
#			name = line[1]
#
#		else:
#			motif_A.append(float(line[1]))
#			motif_C.append(float(line[2]))
#			motif_G.append(float(line[3]))
#			motif_T.append(float(line[4]))


