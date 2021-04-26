import sys 
import os
import operator

#parste similarity file, sodass es ein file gibt das pro DBD  fuer jeden TF den maximalen 
# Wert ueber alle cluster der DBD enthaelt 
#Und ein zweites file das pro DBD die Summe ueber alle TFs enthaelt ( also pro TF den maximal
# Wert 





def main(info_all_DBD_file,similarity_file):
	
	result = {} #enthaelt als schlussel die DBD und als Wert eine Liste von Paaren, Paar besteht aus TF und max Wert ueber alle Cluster der DBD

	sim = open(similarity_file, 'r')
	info_all_DBD = open(info_all_DBD_file, 'r')
	

	#------------------------------------------------
	DBDs = {}# enthaelt DBDs als schluessel und als value liste der cluster
	key = ""
	value = ""
	first = False
	for line in info_all_DBD:
		line = line.strip()
		if line[0] == 'D':
			key = line	
#			print (line)
		else:
			line = line.split('\t')	
			for e in line:
				e = e.strip()
				if e[0] == 'c':
					value =  e + '-'
					first = True
				else:
					if(first):
						value = value  +  e
						first = False			
					else:	
						value = value + ':' +  e
			if DBDs.has_key(key) == True:
				helper = DBDs.get(key)
				DBDs[key] = helper + [value]
			else:
				#print(key)	
				DBDs[key] = [value]
#	print(DBDs)
#	print(len(DBDs))

	print("DBDs gelesen")
	#-----------------------------------------------

	#fill result
	cluster = ""
	TF = ""
	DBD = ""
	for line in sim:
		line = line.strip()
		line = line.split('\t')
		if((line[0][0] == 'c' and line[1][0] !='c') or (line[1][0] == 'c' and line[0][0] !='c')): #checks if line is a line of interest 
#			print("kandidaten: " + line[0] + " " + line[1])
				
			if line[0][0] == 'c':
				cluster = line[0]		
				TF = line[1]
			else: 
				cluster = line[1]
				TF = line[0]

			
			#zu welcher DBD gehoert cluster
			for item in DBDs: 
#				print(item)
				helper = DBDs[item]
				if (helper.count(cluster) > 0):
					DBD = item
#					print("Jippppiiiie " + cluster + " "+ DBD)
					if DBD in result: #wird geschaut ob DBD schon in result vorhanden ist
						helper = result[DBD]
						insert = True
						for a, b in helper:#schaut in allen tupeln ob TF schon vorhanden ist
							if a == TF:	
	#							print("check max: " + b + " " + line[2])
								if float(line[2]) > float(b): #wert des actuellen clusters ist fur TF hoeher 
									helper.remove((a,b)) #loeschen des alten tupels (sind statisch kann man also nich veraendern)
									#insert = True
	#								print(helper)	
								else:
									insert = False
						if insert: #tuple wird mit neuem s max wert wieder eingefuegt /bzw tuple wird eingefuegt wenn TF noch nicht betrachtet wurde
							result[DBD] = helper + [(TF, line[2])]
	#						print("insert new tuple: " + TF + " " + line[2])
					else:
						result[DBD] = [(TF, line[2])] # wenn DBD noch nicht in result vorhanden ist
	#					print("insert new DBD: " + TF + " " + line[2])
	
#	print(result)
	print("Similarity gelesen")
	#------------------------------
	# write output
	# for each DBD and each TF write max value over all cluster
	sorted_TFs = {}
	hihi = []
 	output = open("DBDs.txt", 'w')
 	#output = open("/MMCI/MS/EpiregDeep/work/TFtoMotifs/mosta_src/RandomMotifs/DBDs.txt", 'w')
 	#output = open( "/DBDs.txt", 'w')
	first = True
	for item in result:
#		print(item)
#		if item == "DBD_22":
#			print("ZZZZZZZZZZZZZ")
		helper = result[item]
		helper.sort(key = lambda x : x[0]) #sortiert tuple in alphabetischer reihnfolge der TFs (-> spaeter auch wichtig fuer c++ programm)
		if first:#schreib header
			for a,b in helper:
				if first:
					output.write('.')
				hihi.append(a)
				output.write('\t' + a)
				first = False
			output.write('\n')
		output.write(item )

		for a,b in helper:# fuer alle tuple in der liste 

			#--------
			#fuer zweites output file
			if a in sorted_TFs:
				sorted_TFs[a] = sorted_TFs[a] + float(b)
			else:
				sorted_TFs[a] = float(b)
			#-------

			#a muss man nicht schreiben, weiss ja das es alphabetisch sortiert ist und dementsprechend die header beschriftet hat
			output.write( '\t' + b)
		output.write('\n')	
	output.close()
	#--------------------------------
#	print(sorted_TFs)
#	print(hihi)
	output = open("sum_TFs.txt", 'w')
	#output = open( "/sum_TFs.txt", 'w')
#	output = open("/MMCI/MS/EpiregDeep/work/TFtoMotifs/mosta_src/RandomMotifs/sum_TFs.txt", 'w')
	for item in hihi:
		output.write(item + '\t'+ str(sorted_TFs[item]) + '\n')

	output.close()






if (len(sys.argv))< 3:
    print("Usage python ./parse_result_sstat.py info_all_DBD.txt , result_similarity.txt")
else:
        main(sys.argv[1], sys.argv[2])
