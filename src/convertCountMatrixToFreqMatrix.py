import sys 
import os
import operator

EPSILON = 0.001 

def compansatAddEpsilon(zero_A, A, no_zeros, i ,e):

	if zero_A == False and A[i] > e and no_zeros > 0: 
		A[i] = A[i] - e
		no_zeros = no_zeros - 1
	else:
		zero_A = True


def convertCountToFreq(A,C,G,T):
	
	sum_ = 0.0
	number_of_zeros = 0
	zero_A = False
	zero_C = False
	zero_G = False
	zero_T = False

	for i in range(len(A)):
		sum_ = A[i] + C[i] + G[i] + T[i]
#		print("sum: " + str(sum_))
		if A[i] != 0.0 and (A[i]/sum_) > 2*EPSILON:	
			A[i] = A[i]/ sum_ 
		else:
			A[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_A = True

		if C[i] != 0.0 and (C[i]/sum_) > 2*EPSILON:
			C[i] = C[i]/ sum_ 
		else:
			C[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_C = True

		if G[i] != 0.0 and (G[i]/sum_) > 2*EPSILON:
			G[i] = G[i]/ sum_
		else:
			G[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_G = True

		if T[i] != 0.0 and (T[i]/sum_) > 2*EPSILON:
			T[i] = T[i]/ sum_
		else:
			T[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_T = True
	
#		print(number_of_zeros)	 
		e = (number_of_zeros * EPSILON)/ (4-number_of_zeros)
		no_zeros = 4 - number_of_zeros
		while  no_zeros > 0:
			#if number_of_zeros == 3:
			#	e = (number_of_zeros * EPSILON)
			#else:
			if zero_A == False and A[i] > e and no_zeros > 0: 
				A[i] = A[i] - e
				no_zeros = no_zeros - 1
			else:
				zero_A = True

			if zero_C == False and C[i] > e and no_zeros> 0: 
				C[i] = C[i] - e	
				no_zeros = no_zeros - 1
			else: 
				zero_C = True
	
			if zero_G == False and G[i] > e and no_zeros > 0: 
				G[i] = G[i] - e
				no_zeros = no_zeros - 1
			else:
				zero_G = True

			if zero_T == False and T[i] > e and no_zeros > 0: 
				T[i] = T[i] - e		
				no_zeros = no_zeros - 1
			else:
				zero_T = True

		number_of_zeros = 0	
		zero_A = False
		zero_C = False
		zero_G = False
		zero_T = False
#		print("A: " + str(A))
#		print("C: " + str(C))
#		print("G: " + str(G))
#		print("T: " + str(T))
	


def check(A,C,G,T):
	check = True
	for i in range(len(A)):
		if  A[i] +  C[i] +  G[i] +  T[i] != 1.0:
			check = False
#			print(  "i: " + str(i)+ " " + str(A[i] + C[i] + G[i] +  T[i]))

	return check




def writeOutput(A,output ):

	for i in range(len(A)):
		if i != len(A) -1:
			output.write(str(A[i]) + "\t")
		else: 	
			output.write(str(A[i]) + "\n")
def main(input_file, output_dir):

	countMatrix = open(input_file, 'r')
	
	A = []
	C = []
	G = []
	T = []

	name = ""
	for line in countMatrix:
		line = line.split()
		if line[0] == "ID":
			if name != "":
				convertCountToFreq(A,C,G,T)
#				print("A: " + str(A))
#				print("C: " + str(C))
#				print("G: " + str(G))
#				print("T: " + str(T))
#				if check(A,C,G,T) == False:
#					print ("AAAAAHHHH: " + name)
				output = open(output_dir + "/" + name + ".txt", 'w')
				writeOutput(A, output)
				writeOutput(C, output)
				writeOutput(G, output)
				writeOutput(T, output)
				output.close()
			
			A = []
			C = []
			G = []
			T = []
			name = line[-1]
			#name = line[2]
#			print(name)

		elif line[0] == "XX" or line[0] == "//" or line[0] == "P0" or line[0] == "AC" or line[0] == "DE":
			continue

		else:
			A.append(float(line[1]))
			C.append(float(line[2]))
			G.append(float(line[3]))
			T.append(float(line[4]))
			
		

	convertCountToFreq(A,C,G,T)
#	if check(A,C,G,T) == False:
#		print ("AAAAAHHHH: " + name)
#		print("A: " + str(A))
#		print("C: " + str(C))
#		print("G: " + str(G))
#		print("T: " + str(T))
		
	output = open(output_dir + "/" + name + ".txt", 'w')
	writeOutput(A, output)
	writeOutput(C, output)
	writeOutput(G, output)
	writeOutput(T, output)
	output.close()







if (len(sys.argv))< 3:
    print("Usage python ./convertCounMatrixToFReqMatrix.py humanJaspar_ver7.txt, output_dir")

else:
        main(sys.argv[1], sys.argv[2])
