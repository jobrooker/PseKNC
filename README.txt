PseKNC Python Program Version 1.0 06/02/2014

General Usage Notes
-----------------------------

- This program can be used to calculate pseudo K-tuple nucleotide composition
  for one, or many, query DNA sequences. This compression is useful because a 
  DNA sequence can be represented by a list of numbers, a vector made by 
  looking at physicochemical properties of the K-tuple nucleotides (1).

- Currently, only dinucleotides and trinucleotides can be processed. Different 
  uses of the program can be obtained by entering optional parameters (See 
  parameter descriptions below). 

- Program also includes a "simple" option: frequencies of the oligonucleotides
  (either di- or tri-) can be calculated from a query sequence, without looking
  at any physicochemical properties.

-------------------------------------------------------------------------------

Command Line Arguments
-----------------------------

Example: ./pseknc.py -i test.txt -o out.txt -x propNames.txt -f csv -t 1 -k 2 -j 1 -w 1

Example: ./pseknc.py -i test.txt -o out.txt -x propNames.txt -k 2 -s 

-h,-?,--help  help, will print brief descriptions of the different arguments, including
              possible values

**Required**
-i  input file name, must be in FASTA format. An input file should contain only DNA or only RNA.

     >Example1
     AGTCAGTTATGACATGAC
     >Example2
     AGTCAGTTA
     TGACATGAC
     >Example3
     AGTCAGTTatgacatgac

-o  output file name

-x  property file, containing a list of properties to be used in calculations

     *** NOTE: Program is case-sensitive for query properties. 'tilt' instead of 'Tilt' 
     will not work. To avoid spelling errors, copying and pasting property names into 
     the query propery file (ex. propNames.txt) will have the best results. ***

      Example:
	Tilt
	Shift
	...

**Optional**
-f change the format of the output file
      tab - output delimited by tab 
      svm - the libSVM training data format
      csv - format that can be loaded into a spreadsheet program (DEFAULT)

-t Type of PseKNC
      1 - Type I PseKNC, used to represent a DNA sequence with a vector 
          containing (4^k + λ) items (DEFAULT)
      2 - Type II PseKNC, used to represent a DNA sequence with a vector
          containing (4^k + λΛ) items, where Λ represents the number of 
          physical properties being analyzed

-k K-tuple, the length of the oligonucleotides
      2 - dinucleotides (ex. AA) (DEFAULT)
      3 - trinucleotides (ex. AAA)

-p Will list the text files including with the program that list the 
    physicochemical properties, for which data is available for use in this 
    program. Nothing needs to be entered after '-p'.
      
-j Sequence order correlation, represented by λ or j. λ must be a whole number
    (ex. 1,2,3,...) and cannot be larger than L - k,
    (length of the query sequence) - (the length of the oligonucleotides).
      1 - (DEFAULT)

-w Weight factor can be any value between (0,1].
      1 - (DEFAULT)

-s "Simple" PseKNC only calculates the frequency of all the different oligo-
    nucleotides in the query sequence (k must also be specified, or left as the
    default value of 2). If this argument is entered, the regular pseKNC 
    program will not run and the results of the "simple" calculation will be 
    written to the specified output file. As with '-p', nothing else needs 
    to be entered after '-s'.


-----------------------------
EXAMPLES
-----------------------------

Demo Run #1
-----------------------------

-i test.txt
       
      >Example1
      AGTCAGTTATGACATGAC
      >Example2
      AGTCAGTTA
      TGACATGAC
      >Example3
      AGTCAGTTatgacatgac

-o out.txt (will be created by the program during its run)

-x propNames.txt

      MW-Daltons
      Nucleosome-Rigid

-k 3

      k = 3, so the query sequence will be split into trinucleotides

-w 0.5

      w = 0.5, the weight factor in the program's algorithm will be 0.5

Make sure you are in the directory (in this case, 'dir') where the program is.
At the command line:

 user:~/dir$ ./pseknc.py -i test.txt -o out.txt -x propNames.txt -k 3 -w 0.5

The program will print 'DONE.' and the results will be printed to the output
file named in the function call (in this case, out.txt).

 user:~/dir$ ./pseknc.py -i test.txt -o out.txt -x propNames.txt -k 3 -w 0.5
 DONE.
 user:~/dir$

The output file generated from this command should look like this (csv is the 
default for the output format):

       >Example1
       0.0,0.0,0.0,0.0,0.0,0.0,0.067,...
       >Example2
       0.0,0.0,0.0,0.0,0.0,0.0,0.067,...
       >Example3
       0.0,0.0,0.0,0.0,0.0,0.0,0.067,...

Installation
-----------------------------

Download the program file to a directory. Command line arguments must be called 
from this directory to run the program.

Questions
-----------------------------
Questions about use or content can be sent to:

Jordan Brooker
Department of Computer Science
Vassar College, Poughkeepsie, NY 12604, USA
jrdnbrkr@gmail.com


References:

1. Chen, W., Lei, T., Jin, D., Lin, H., and Chou, K., 2014, PseKNC: A flexible
	web server for generating pseudo K-tuple nucleotide composition: 
	Analytical Biochemistry, v. 456, p. 53-60.








