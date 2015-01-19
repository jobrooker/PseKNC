#!/usr/bin/python

# -----------------------------
#  Project 2 - Three Autocorrelations
#  Jordan Brooker
#  06/11/14 
#  jrdnbrkr@gmail.com
# -----------------------------

from __future__ import division

import math,string,re,sys,getopt

#  _mean: Calculates mean value of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the mean value of listy

def _mean(listy):
	return sum(listy)/len(listy)

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
	res = math.sqrt(sum(temp)/len(listy))
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
# input: PropFileName = string, and prop = string of one property
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
# input: olinuc = string, prop = string, propFileName = string
# output: value = an int that is property value of an olinuc

def getSpecificValue(olinuc,olinucs,prop,supInfo):
	#olinucs = getNucSeq(SupFileName).split(",")
	olinucs = olinucs.split(",")
	values = getValues(prop,supInfo).rstrip()
	values = values.split(",")
	#valueS = [float(x) for x in values.split(",")]
	count = olinucs.index(olinuc)
	value = values[count]
	return float(value)

#  avgP
# --------------------------------------
# inputs: seq = string,length = int,k = int, prop = string, SupFileName = string
# outputs: sum = int

def avgP(seq,olinucs,length,k,prop,supInfo):
	limit = length - k + 1
	i = 1
	sum = 0
	while i < limit or i == limit:
		#value = hn(seq[i - 1],prop,SupFileName)
		value = getSpecificValue(seq[i - 1],olinucs,prop,supInfo)
		sum = sum + value
		i = i + 1
	sum = sum/limit
	return sum

#  Moreau
# -------------------------------------
# inputs: seq = string, length = int, k = int, l = int, prop = string, 
#         supFileName = string
# output: final = int

def moreau(seq,olinucs,length,k,l,prop,supInfo):
	limit = length - k - l + 1
	d = 1
	prod = 0
	while d < limit or d == limit:
		current = getSpecificValue(seq[d-1],olinucs,prop,supInfo)
		#hn(seq[d-1],prop,SupFileName)
		next = getSpecificValue(seq[d+l-1],olinucs,prop,supInfo)
		#hn(seq[d+l-1],prop,SupFileName)
		prod = prod + (current * next)
		d = d + 1
	final = prod / limit
	return final

#  geary
# --------------------------------------
# inputs: seq = string, length = int, k = int, l = int, prop = string,
#         SupFileName = string
# output: final = int

def geary(seq,olinucs,length,k,l,prop,supInfo):
	lim = length - k + 1
	limit = length - k - l + 1
	b = 1
	sqr = 0
	while b < limit or b == limit:
		current = getSpecificValue(seq[b-1],olinucs,prop,supInfo)
		#hn(seq[b-1],prop,SupFileName)
		next = getSpecificValue(seq[b+l-1],olinucs,prop,supInfo)
		#hn(seq[b+l-1],prop,SupFileName)
		sqr = sqr + ((current - next) * (current - next))
		b = b + 1
	top = sqr * lim
	limit2 = (length - k - l + 1)
	c = 1
	sqr2 = 0
	while c < limit2 or c == limit2:
		#current = hn(seq[c-1],prop,SupFileName)
		current = getSpecificValue(seq[c-1],olinucs,prop,supInfo)
		avg = avgP(seq,olinucs,length,k,prop,supInfo)
		sqr2 = sqr2 + (current - avg) * (current - avg)
		c = c + 1
	bottom = sqr2 * limit * 2
	final = float((top/bottom)*1000)/1000.0
	return final
		

#  moran
# --------------------------------------
# inputs: seq = string, length = int, k = int, l = int, prop = string,
#         SupFileName = string
# output: final = int

def moran(seq,olinucs,length,k,l,prop,supInfo):
	limit = length - k - l + 1
	j = 1
	top = 0
	avg = avgP(seq,olinucs,length,k,prop,supInfo)
	while j < limit or j == limit:
		current = getSpecificValue(seq[j-1],olinucs,prop,supInfo)
		#hn(seq[j-1],prop,SupFileName) 
		partOne = current - avg
		next = getSpecificValue(seq[j+l-1],olinucs,prop,supInfo)
		#hn(seq[j+l-1],prop,SupFileName)
		partTwo = next - avg
		top = top + (partOne * partTwo)
		j = j + 1
	top = top/limit 
	limit2 = length - k + 1
	bottom = 0
	b = 1
	while b < limit2 or b == limit2:
		current = getSpecificValue(seq[b-1],olinucs,prop,supInfo)
		#hn(seq[b-1],prop,SupFileName)
		bottom = bottom + ((current - avg) * (current - avg))
		b = b + 1
	bottom = bottom/limit2
	final = top/bottom
	return final 

#  autocorrelation
# -------------------------------------------
# inputs: autoc = string, sequence = string, k = int, l = int,
#         propFileName = string
# output: values = list of ints

def autocorrelation(autoc,inputfile,outputFile,props,k,l,geneticMaterial):
	l = int(l)
	InFileName = inputfile
	OutFileName = outputFile
	InFile = open(InFileName,'r')	
	OutFile = open(OutFileName,'w')
	sequence = ''
	sequences = []
	label = ''
	labels = []
	#Getting sequences
	for line in InFile:
		f = re.search('(\>.+)',line)
		g = re.search('(^[A*a*T*t*C*c*G*g*]*$)',line)
		if f and not(sequence == ''):
			label = f.group(1)
			labels.append(label)
			sequences.append(sequence)
			sequence = ''
		elif f:
			label = f.group(1)
			labels.append(label)
		elif g and not(g.group(1)==''):
			sequence = sequence + g.group(1).upper()
	sequences.append(sequence)
	# Getting supporting info from files
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
	#Writing to output file
	m = 0
	for sequence in sequences:
		label = labels[m]
		OutFile.write(label)
		length = len(sequence)
		OutFile.write('\n' + sequence + '\n')
		OutFile.write("----------------------------------\n")
		seq = sepSequence(sequence,k)
		if k == 2:
			SupFileName = 'Supporting Information S1' + geneticMaterial + '.txt'
		else:
			SupFileName = 'Supporting Information S3 DNA.txt'
		OutFile.write("Properties: " + str(props) + '\n')
		autocs = autoc.split(",")
        for auto in autocs:
            OutFile.write(auto.upper() + ": ")
            values = []
            for prop in props:
                if auto in ('moran','Moran'):
                    #value = float("%.3f" % moran(seq,olinucs,length,k,l,prop,supInfo))
                    value1 = str(moran(seq,olinucs,length,k,l,prop,supInfo))
                    value2 = ""
                    b = re.search('(.+)\.(\d\d\d)', value1)
                    if b:
                        value2 = b.group(1) + "." + b.group(2)
                    values.append(value2)
                elif auto in ('geary','Geary'):
                    #value = float("%.3f" % geary(seq,olinucs,length,k,l,prop,supInfo))
                    value1 = str(geary(seq,olinucs,length,k,l,prop,supInfo))
                    value2 = ""
                    c = re.search('(.+)\.(\d\d\d)',value1)
                    if c:
                        value2 = c.group(1) + "." + c.group(2)
                        
                    values.append(value2)
                elif auto in ('moreau','Moreau'):
                    #value = float("%.3f" % moreau(seq,olinucs,length,k,l,prop,supInfo))
                    value1 = str(moreau(seq,olinucs,length,k,l,prop,supInfo))
                    value2 = ""
                    d = re.search('(.+)\.(\d\d\d)', value1)
                    if d:
                        value2 = d.group(1) + "." + d.group(2)
                    values.append(value2)
            OutFile.write(str(values) + '\n')
		OutFile.write('======================================\n')
		m += 1
	InFile.close()
	OutFile.close()
	

#  main: Gets arguments from the command line
# -----------------------------------------

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h?'help'a:i:o:k:j:x:p",[])
	except getopt.GetoptError:
		print './autocorrelations.py -a <Geary,Moran,and/or Moreau> -i <inputfile> -o <outputfile> -x <propertiesfile> '
		sys.exit(2)
	autoc = 'geary,moran,moreau'
	kay = 2
	jay = 1
	geneticMaterial = 'DNA'
	for opt,arg in opts:
		if opt == '-h' or opt == '-?' or opt == '-help':
			print ""
			print "autocorrelations.py [options]"
			print "options: "
			print "-a  One or more autocorrelations to be calculated"
			print "      Geary"
			print "      Moran"
			print "      Moreau"
			print "-i  The input file, in valid FASTA format"
			print "-o  The results will be printed to this output file."
			print "-x  File containing list of properties to be used in calculations,"
			print "      with each property on a separate line."
			print "      Example:"
			print "        Tilt"
			print "        Shift"
			print "        ..."
			print "-k  Kind of oligonucleotide (default = 2):"
			print "      2 - Dinucleotide"
			print "      3 - Trinucleotide"
			print "-j  lambda (default = 1)"
			print "-p  Lists the file names of physicochemical property data files."
			print "    Takes no input."
			print "More detailed parameter descriptions can be found in the README file."
			print ""
			sys.exit()
		elif opt == '-a':
			for a in arg.split(","):
				if a in ('Geary','geary','Moran','moran','Moreau','moreau'):
					autoc = arg
				else:
					print "  -a: Please only enter one or more of these autocorrelations:"
					print "      Geary"
					print "      Moran"
					print "      Moreau"
					sys.exit()
		elif opt == '-i':
			InFileName = arg
			InFile = open(InFileName,'r')
			y = -1
			for line in InFile:
				g = re.search('\>',line)
				h = re.search('^[A*a*C*c*G*g*T*t*]*$', line)
				i = re.search('^[A*a*C*c*G*g*U*u*]*$',line)
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
					print "      characters, [AaTtCcGg] OR [AaUuCcGg], and must be in FASTA format"
					sys.exit()
			if z == 0:
				inputFile = arg
			else:
				print "  -i: One of the ID lines in the input file has no sequence following it."
				sys.exit()
		#Check if in SupFiles?
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
		elif opt == '-k':
			if arg == '2' or arg == '3':
				kay = int(arg)
			else:
				print "  -k: The k-tuple argument must be 2 or 3"
				print "       (dinucleotides or trinucleotides)."
				sys.exit()
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
			if (float(arg) % 1) == 0 and float(arg) >= 1 and float(arg) < (ell - kay):
				jay = arg
			else:
				print "  -j: Lambda must be a whole number and smaller"
				print "       than the length of the query sequence minus"
				print "       the k-tuple number (lambda < L-k)."
				print "    Length of query sequence: ",ell
				print "    k: ",kay
				sys.exit()
		elif opt == '-p':
			print "Physicochemical Property Text Files:"
			print "Dinucleotides: Supporting Information S1 DNA.txt OR"
			print "               Supporting Information S2 RNA.txt"
			print "Trinucleotides: Supporting Information S3.txt"
			sys.exit()
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
	autocorrelation(autoc,inputFile,outputFile,props,kay,jay,geneticMaterial)
	print "DONE."
	
if __name__ == "__main__":
	main(sys.argv[1:])
""



































