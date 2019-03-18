import sys 
import os
import operator

def main(info_file , all_TFs_file):

	info = open(info_file, 'r')
	all_TFs = open(all_TFs_file, 'r')
	output = open( "TF_to_DBD.txt", 'w')

	#read all TFs 
	sort_TFs = {}
	first_line = True
	for line in all_TFs:
		if first_line == True:
			first_line = False
		else:
			line = line.strip()
			line = line.split('\t')
		#	line = line[:-4]
			sort_TFs[line[1].upper()] = line[0]
	
#	sort_TFs.sort() #sort them lexicograical
#	print(sort_TFs)
	#read DBDs	
	DBDs = {}
	DBD = ""
	TFs = []
	first = True
	for line in info:
		line = line.strip()
		line = line.split('\t')
		if line[0][0] == 'D':
			if not first: 
				DBDs[DBD] = TFs
				TFs = []
			else: 
				first = False
			DBD = line[0]
		else:
			for item in line:
				if item[0] != 'c':
					TFs.append(item)

	DBDs[DBD] = TFs
	print(DBDs)
	
	#for output JASPAR TFs	
	for item in sort_TFs.keys():
		value = sort_TFs[item]
		#print(item)
		#output.write(item + '\t')
		
		for item_2 in DBDs:
		#	print (item_2)
			helper = DBDs[item_2]
			if item in helper:
			#for e in helper:
				#if item == e:
				output.write(value + '\t' +item_2 + '\t-1\n')


	motifs = []
	for filename in os.listdir("PFMs/"):
    		if filename.endswith(".mat"): 
			filename = filename[:-4]
			motifs.append(filename)	


	motifs.sort() #sort them lexicograical

	#for output motifs	
	counter = 0
	for item in motifs:
		output.write(item + '\tX\t' + str(counter) +  '\n')
		counter = counter + 1


if (len(sys.argv))< 4:
    print("Usage python ./TF_DBDs.py info_all_DBDs.txt , ENSG_HNCG.tx")
else:
        main(sys.argv[1], sys.argv[2])
