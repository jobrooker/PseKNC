#!/usr/bin/python

# -----------------------------------------
#  PseKNC - 06/11/2014
#  Jordan Brooker 
#  Jrdnbrkr@gmail.com
# -----------------------------------------

from __future__ import division

import math,string,re,sys,getopt

#  _mean: Calculates mean value of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the mean value of listy

def _mean(listy):
	return float(sum(listy))/len(listy)

#  _std(listy): Calculates standard deviation of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the standard deviation of listy
 
def _std(listy,ddof=1):
	mean = _mean(listy)
	# for i in listy, (i - mean)^2 => list
	temp=[math.pow(i-mean,2) for i in listy]
	# temp is added together, divided by length of listy,
	# and the square root is taken
	res = math.sqrt(float(sum(temp))/len(listy))
	return res

#  sepSequence
# -----------------------------------------
# inputs: seq, string and k, int
# output: list of k-tuples

def sepSequence(seq, k):
	i = k-1
	seqq = []
	while i < len(seq):
		j = 0
		nuc = ''
		while j < k:
			nuc = seq[i-j] + nuc
			j = j + 1
		seqq.append(nuc)
		i += 1
	return seqq

#  getValues: Returns a line of values for one property
#             from physicochemical property files
# -----------------------------------------
# input: prop = string of one property and supInfo = string
# output: values = a string representing a list of property values

def getValues(prop,supInfo):
	values = ""
	name = re.search(prop,supInfo)
	if name:
		strr = prop + '\s*\,(.+)'
		b = re.search(strr, supInfo)
		if b:
			values = b.group(1)
	return values

#  getSpecificValue: Returns a property value for a specific di- or tri-
#                    nucleotide
# -----------------------------------------
# input: olinuc = string, prop = string, supInfo = string
# output: value = an int that is property value of an olinuc

def getSpecificValue(olinuc,olinucs,prop,values,supInfo):
	values = values.split(",")
	#valueS = [float(x) for x in values.split(",")]
	count = olinucs.index(olinuc)
	value = values[count]
	return float(value)

#  hn: Hn function 
# -----------------------------------------
# inputs: olinuc = string, prop = string, supInfo = string
# output: temp = int

def hn(olinuc,olinucs,prop,supInfo,values):
	#values = getValues(prop,supInfo).rstrip()
	h0 = float(getSpecificValue(olinuc,olinucs,prop,values,supInfo)) 
	valueS = [float(x) for x in values.split(",")]
	temp = float((h0 - _mean(valueS)) / _std(valueS))
	return temp

#  theta2: Theta(i,i+j) function, for Type I
# -----------------------------------------
# input: seq = string, props = string, i = int,
#	 Type = int, k = int, j = int, supInfo = string
# output: summ = int

def theta2(seq,olinucs,props,i,k,j,supInfo):
	summ = 0
	values = ''
	for prop in props:
		values = getValues(prop,supInfo).rstrip()
		hn1 = hn(seq[i-1],olinucs,prop,supInfo,values)
		hn2 = hn(seq[i+j-1],olinucs,prop,supInfo,values)
		subsqr = math.pow(float(hn1-hn2),2)
		summ = summ + subsqr
	return float(summ)/len(props)

#  J: J(i,i+j) function, for Type II 
# --------------------------------------
# inputs: seq = string, prop = string, i = int, 
#         k = int, j = int, supInfo = string
# output: product = int

def J(seq,olinucs,prop,i,k,j,supInfo):
	values = getValues(prop,supInfo)
	hn1 = hn(seq[i-1],olinucs,prop,supInfo,values)
	hn2 = hn(seq[i+j-1],olinucs,prop,supInfo,values)
	return float(hn1*hn2)

#  theta1: Theta(j) and Tau(LamGam) function
# -----------------------------------------
# input: seq = string, props = string, Type = int
#        k = int, j = int, supInfo = string
# output: final = int

def theta1(seq,olinucs,props,Type,k,j,supInfo):
	k = int(k)
	gamma = len(props)
	seqq = sepSequence(seq,k)
	i = 1
	a = 0
	if Type == 1:
		var = len(seq) - int(k) - int(j) + 1
		while i <= var:
			b = 0
			b = theta2(seqq,olinucs,props,i,k,j,supInfo)
			a = a + b
			i = i + 1
	else:
		ii = 0
		var = len(seq) - int(k) - int(j/gamma)
		while i <= var:
			b = 0
			if ii == gamma:
				ii = 0
				b = J(seqq,olinucs,props[ii],i,k,int(j/gamma),supInfo) 
				a = a + b
			else:
				b = J(seqq,olinucs,props[ii],i,k,int(j/gamma),supInfo)
				a = a + b
			ii = ii + 1
			i = i + 1
	final = float(a)/var
	return final

#  pseKNCHelper: Creates list of adjusted frequency values for 4^k 
#  oligonucleotides and the lambda terms
# -----------------------------------------
# input: seq = string, props = string, Type = int, k = int, j = int, w = int
# output: freqs = list of ints

def pseKNCHelper(seq,olinucs,props,Type,k,j,w,geneticMaterial,supInfo):
	gamma = len(props)
	freqs = []
	seqq = sepSequence(seq,k)
	olinucs = olinucs.split(",")
	for olinuc in olinucs:
		freq = seqq.count(olinuc)
		freq = float(freq/len(seqq))
		total = 0
		i = 1
		if Type == 2:
			j = j/gamma
		while i <= j:
			total = total + theta1(seq,olinucs,props,Type,k,i,supInfo)
			i = i + 1
		total = float(freq/(1 + (float(w)*total)))
		total = int(total*1000) / 1000.0
		freqs.append(total)
	#Computing Lambda terms...
	fourK = math.pow(4,k)
	mu = fourK + 1
	while (fourK+1) <= mu <= (fourK + j):
		top = float(w) * theta1(seq,olinucs,props,Type,k,int(mu-fourK),supInfo)
		bottomTheta = 0
		bottom = 0
		i = 1
		while 1 <= i <= j:
			bottomTheta += theta1(seq,olinucs,props,Type,k,i,supInfo)
			i += 1
		bottom = 1 + (float(w) * bottomTheta)
		term = float(top / bottom)
		term = int(term * 1000) / 1000.0
		freqs.append(term)
		mu += 1
	return freqs

#  pseKNC: Opens input and output files, calls functions to calculate
#  values and writes outputs to the output file in appropriate format
# -----------------------------------------
# input: inputFile = string, outputFile = string, propNames = string,
#        Type = int, k = int, j = int, w = int, formatt = string
# output: nothing, write to output file

def pseKNC(inputFile,outputFile,props,Type,k,j,w,formatt,geneticMaterial):
	j = float(j)
	listy = ''
	output = ''
	outputs = ''
	gamma = len(props)
	#Getting supporting info from files
	if k == 2:
		SupFileName = 'Supporting Information S1 ' + geneticMaterial + '.txt'
	else:
		SupFileName = 'Supporting Information S3 DNA.txt'
	SupFile = open(SupFileName,'r')
	supInfo = SupFile.read()
	o = re.search('Physicochemical properties\,(.+)\n',supInfo)
	olinucs = ''
	if o:		
		olinucs = o.group(1).rstrip()
	SupFile.close()
	# Calculating frequencies
	tooples = getInputSeqs(inputFile)
	for toople in tooples:	
		if Type == 1:
			listy = pseKNCHelper(toople[1],olinucs,props,Type,k,j,w,geneticMaterial,supInfo)
		else:
			listy = pseKNCHelper(toople[1],olinucs,props,Type,k,j*gamma,w,geneticMaterial,supInfo)
		i = 0
		z = 1
		listz = ''
		# Generating output file format
		while i < len(listy):
			if formatt == 'csv':
				if i < 1:
					listz = listz + str(listy[i])
					i = i + 1
				else:
					listz = listz + "," + str(listy[i])
					i = i + 1
			elif formatt == 'tab':
				if i < 1:
					listz = listz + str(listy[i])
					i = i + 1
				else:
					listz = listz + "\t" + str(listy[i])
					i = i + 1
			elif formatt == 'svm':
				if i < 1:
					listz = listz + str(z) + ":" + str(listy[i])
					z = z + 1
					i = i + 1
				else:
					listz = listz + " " + str(z) + ":" + str(listy[i])
					z = z + 1
					i = i + 1
		output = toople[0] + listz + '\n'
		outputs = outputs + output
	# Writing to output file
	OutFileName = outputFile
	OutFile = open(OutFileName, 'w')
	OutFile.write(outputs)
	OutFile.close()

#  getInputSeqs
# ----------------------------------------

#Think about how to return labels and sequences
def getInputSeqs(inputFile):
	label = ''
	labels = []
	sequence = ''
	sequences = []
	tooples = []
	InFileName = inputFile
	InFile = open(InFileName, 'r')
	for line in InFile:
		f = re.search('\>',line)
		g = re.search('(^[A*a*C*c*T*t*G*g*]*$)',line)
		h = re.search('(^[A*a*C*c*G*g*U*u*]*$)',line)
		if f and not(sequence == ''):
			label = line
			labels.append(label)
			sequences.append(sequence)
			sequence = ''
		elif f:
			label = line
			labels.append(label)
		elif g and not(g.group(1)==''):
			sequence = sequence + g.group(1).upper()
		elif h and not(h.group(1) == ''):
			sequence = sequence + h.group(1).upper()
	sequences.append(sequence)
	for label in labels:
		tooples.append((label,sequences[labels.index(label)]))
	return tooples
	InFile.close()

#  generate_permutations
# ----------------------------------------
# inputs: chars = int
# outputs: all possible oligonucleotides of length k

def generate_permutations(chars):
	allowed_chars = ['A','T','G','C']
    	status = []
    	for tmp in range(chars) :
        	status.append(0)
    	last_char = len(allowed_chars)
    	rows = []
    	for x in xrange(last_char ** chars) :
       		rows.append("")
        	for y in range(chars - 1 , -1, -1) :
            		key = status[y]
            		rows[x] = allowed_chars[key] + rows[x]
        	for pos in range(chars - 1, -1, -1) :
           		if(status[pos] == last_char - 1) :
                		status[pos] = 0
          		else :
                		status[pos] += 1
                		break;
        
    	return rows

#  simplePseKNC: Calculates the frequencies of all possible
#  oligonucleotides in the input sequence
# ----------------------------------------
# inputs: inputFile = a string, outputFile = a string, k = int, 
#         formatt = string
# output: nothing, writing to a file

def simplePseKNC(inputFile,outputFile,k,formatt,geneticMaterial):
	OutFileName = outputFile
	OutFile = open(OutFileName, 'w')
	tooples = getInputSeqs(inputFile)
	listy = ''
	# Need to generate all possible oligonucleotides
	olinucz = generate_permutations(k)
	for toople in tooples:
		OutFile.write(toople[0] + toople[1] + '\n')
		listy = sepSequence(toople[1],k)
		z = 1
		for olinuc in olinucz:
			listz = ''
			freq = listy.count(olinuc)
			freq = int(freq/len(listy)*1000)/1000.0
			if formatt == 'csv':
					if olinucz.index(olinuc) < 1:
						print olinucz.index(olinuc)
						print "olinuc: ",olinuc
						listz = listz + str(freq)
					else:
						listz = listz + "," + str(freq)
			elif formatt == 'tab':
					if olinucz.index(olinuc) < 1:
						listz = listz + str(freq)
					else:
						listz = listz + "\t" + str(freq)
			elif formatt == 'svm':
					if olinucz.index(olinuc)< 1:
						listz = listz + str(z) + ":" + str(freq)
						z = z + 1
					else:
						listz = listz + " " + str(z) + ":" + str(freq)
						z = z + 1
			OutFile.write(listz)
		OutFile.write('\n')
	OutFile.close()

#  main: Gets arguments from the command line, calls either
#  simpePseKNC or pseKNC
# -----------------------------------------

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h?'help'i:x:o:t:k:w:j:f:ps",[])
	except getopt.GetoptError:
		print 'pseknc.py -i <inputfile> -o <outputfile> -x <propertiesfile>'
		sys.exit(2)
	Type = 1
	kay = 2
	weight = 1
	lam = 1
	formatt = 'csv'
	s = 0
	u = 0
	inputFile = 'test.txt'
	outputFile = 'out.txt'
	geneticMaterial = 'DNA'
	props = ''
	for opt,arg in opts:
		if opt == '-h' or opt == '-?' or opt == '-help':
			print ""
			print "pseknc.py [options]"
			print "options: "
			print "-i   The input file, in valid FASTA format"
			print "-o   The results will be printed to this output file."
			print "-x   File containing list of properties to be used in calculations,"
			print "      with each property on a separate line."
			print "      Example:"
			print "        Tilt"
			print "        Shift"
			print "        ..."
			print "-f  The output format (default = csv):"
			print "      tab -- Simple format, delimited by TAB. This is designated "
			print "              as the default format."
			print "      svm -- The libSVM training data format."
			print "      csv -- The format that can be loaded into a spreadsheet program."
			print "-t  PseKNC_Type (default = 1):"
			print "      1 - Type 1 PseKNC"
			print "      2 - Type 2 PseKNC"
			print "-k  Kind of oligonucleotide (default = 2):"
			print "      2 - Dinucleotide"
			print "      3 - Trinucleotide"
			print "-j  Set the value of lambda parameter in the PseKNC algorithm."
			print "      Must be smaller than the length of any query sequence (default = 1)."
			print "-w  Set the value of weight parameter in the PseKNC algorithm."
			print "      It can be a value between (0,1] (default = 1.0)."
			print "-s  Calculate only the frequency of each oligonucleotide in the"
			print "      input sequence. Unless otherwise specified, the default value"
			print "      of k is 2."
			print "-p  Lists the file names of physicochemical data information."
			print "More detailed parameter descriptions can be found in the README file."
			print ""
			sys.exit()
		elif opt == '-i':
			InFileName = arg
			InFile = open(InFileName,'r')
			y = -1
			for line in InFile:
				g = re.search('\>',line)
				h = re.search('^[A*a*C*c*G*g*T*t*]+$', line)
				i = re.search('^[A*a*C*c*G*g*U*u*]+$',line)
				if g:
					z = 1
				elif h:
					if not(y == 0):
						if z == 1:
							z = 0
					else:
						print "  -i: The input file sequences contain DNA and RNA. Please only enter either DNA or RNA at a time."
						sys.exit()
					y = 1
					geneticMaterial = 'DNA'
				elif i:
					if not(y == 1):
						if z == 1:
							z = 0
					else:
						print "  -i: The input file sequences contain DNA and RNA. Please only enter either DNA or RNA at a time."
						sys.exit()
					y = 0
					geneticMaterial = 'RNA'
				else:
					print "  -i: The input file sequences can only include valid nucleotide"
					print "      characters, [AaTtCcGg] OR [AaUuCcGg], and must be in FASTA format."
					print "      Make sure that there are no empty lines at the end of the file."
					sys.exit()
			if z == 0:
				inputFile = arg
			else:
				print "  -i: One of the ID lines in the input file has no sequence following it."
				sys.exit()
		elif opt == '-x':
			# Getting list of properties (ex. Tilt, shift)
			properties = ''
			props = []
			propFileName = arg
			propFile = open(propFileName, 'r')
			for line in propFile:
				e = re.search('(.+)', line)
				if e:
					properties = properties + e.group(1) + ","
			propFile.close()
			props = properties.split(",")
			props = props[0:len(props)-1]
		elif opt == '-o':
			outputFile = arg
		# Checking if type argument is 1 or 2
		elif opt == '-t':
			if arg == '1' or arg == '2':
				Type = arg
			else:
				print "  -t: The Type argument must be 1 or 2."
				sys.exit()
		# Checking if k argument is 2 or 3
		elif opt == '-k':
			kay = int(arg)
		# Checking if weight argument is between 0 and 1.0
		elif opt == '-w':
			if 0 < float(arg) <= 1:
				weight = arg
			else:
				print "  -w: The weight factor argument must be "
				print "       between 0.1 and 1.0."
				sys.exit()
		# Checking if the lambda argument is a whole number and smaller than L-k
		elif opt == '-j':
			InFileName = inputFile
			InFile = open(InFileName,'r')
			listy = ''
			for line in InFile:
				g = re.search('\>',line)
				if g:
					label = line
				else:
					listy = listy + line.rstrip()
			InFile.close()
			ell = len(listy)
			if (float(arg) % 1) == 0 and (1 <= float(arg) < (ell - kay)):
				lam = arg
			else:
				print "  -j: Lambda must be a whole number and smaller"
				print "       than the length of the query sequence minus"
				print "       the k-tuple number (lambda < L-k)."
				print "    Length of query sequence: ",ell
				print "    k: ",kay
				sys.exit()
		# Checking if valid format option is entered
		elif opt == '-f':
			if arg in ('csv','tab','svm'):
				formatt = arg
			else:
				print "  -f: The output format can be"
				print "        csv, comma delimited"
				print "        tab, tab delimited"
				print "        svm, libSVM training data format"
				sys.exit()
		elif opt == '-p':
			print "Physicochemical Property Text Files:"
			print "Dinucleotides: Supporting Information S1 DNA.txt OR"
			print "		      Supporting Information S1 RNA.txt"
			print "Trinucleotides: Supporting Information S3.txt"
			sys.exit()
		elif opt == '-s':
			s = 1
	# Checking if -k and -x arguments match
	if kay == 2:
		SupFileName = "Supporting Information S1 " + geneticMaterial + ".txt"
	else:
		SupFileName = "Supporting Information S3 DNA.txt"
	SupFile = open(SupFileName,'r')
	tt = 0
	for prop in props:
		for line in SupFile:
			t = re.search(prop,line)
			if t:
				tt=1
		if not(tt==1):
			print "  '" + prop + "' was not found in the supporting information file."
			print "  Please check that the k value and the query properties correspond."
			sys.exit()
	SupFile.close()
	if s == 1:
		simplePseKNC(inputFile,outputFile,int(kay),formatt,geneticMaterial)
	else:
		if kay == 2 or kay == 3:
				kay = int(kay)
		else:
			print "  -k: The k-tuple argument must be 2 or 3"
			print "       (dinucleotides or trinucleotides)."
			sys.exit()
		pseKNC(inputFile,outputFile,props,Type,kay,lam,weight,formatt,geneticMaterial)
	print "DONE."
	

if __name__ == "__main__":
	main(sys.argv[1:])

