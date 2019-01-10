import sys
import os
import operator


def readPvalues(info_file):

	result = {}
	TF = ""
	helper = {}

	for line in info_file:
		line = line.strip()
		if not line:
			continue
		else:
			line = line.split("\t")
			if line[0] == "TF:":
				if TF != "":
					result[TF] = helper
					helper = {}
				TF = line[1]
				print(TF)
			elif line[0] == ".":
				continue
			else:
				helper[line[1]] = line[2]							

	

	result[TF] = helper
	#print(result)
	return result
			

def main(result_PASTAA, result_CentriMo, result_domainInfo, outputDir):


	PASTAA_file = open(result_PASTAA, "r")
	centrimo_file = open(result_CentriMo, "r")
	domainInfo_file = open(result_domainInfo, "r")

	PASTAA = readPvalues(PASTAA_file)
	centrimo = readPvalues(centrimo_file)	
	domainInfo = readPvalues(domainInfo_file)	

#	TFs = domainInfo.keys()
	TFs = PASTAA.keys()
	for k in TFs:
		output = open(outputDir + "/" + k + ".txt", "w")
		output.write("motif\tPASTAA\tCentriMo\tDomainInfo\n")
		p = PASTAA[k]
		c = centrimo[k]	
		d = domainInfo[k]
		for motifs in p.keys():
		#for motifs in d.keys():
			if motifs in c and  motifs in d:
			#if motifs in d:
			#	output.write(motifs + "\t" + p[motifs] + "\t" + c[motifs] + "\t" + d[motifs] + "\n")
				#output.write(motifs + "\t" + p[motifs] + "\t" + d[motifs] + "\n")
				output.write(motifs + "\t" + p[motifs] + "\t" + c[motifs] +  "\n")
				#output.write(motifs + "\t" + c[motifs] + "\t" + d[motifs] + "\n")
			else:
				if motifs not in c: 
					print("Centrimo\tTF: " + k + " motif: " + motifs)
			#		output.write(motifs + "\t" + p[motifs] + "\t1.0\t" + d[motifs] + "\n")
					output.write(motifs + "\t" + p[motifs] + "\t1.0" + "\n")
			#		#output.write(motifs + "\t1.0\t" + d[motifs] + "\n")
				else:
					print("DomainInfo\tTF: " + k + " motif: " + motifs)
		output.close()


if(len(sys.argv))< 5:

	print("Usage python ./combinePvalues.py  result_PASTAA, result_CentriMo, result_DomainInfo ,output")
else:
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
