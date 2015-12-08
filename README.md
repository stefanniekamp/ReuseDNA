##The MIT License

Copyright (c) 2015 University of California, San Francisco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[http://www.opensource.org/licenses/mit-license.php](http://www.opensource.org/licenses/mit-license.php)

##Overview

####Extract from publication:

Using [caDNAno](https://github.com/sdouglas/cadnano2), a computer-aided design tool for DNA origami, we first routed the scaffold to approximate a 3D shape and exported the design from caDNAno as a text file in JSON format. Next, we used a custom Python script (provided here) and input the JSON file and desired number of unique staples (e.g. 10). We included the option to specify specific staple sequences as input, or to input a range of desired staple lengths and automatically generate sequences. The script determines each custom scaffold sequence by generating a random staple layout, repeatedly assigning a set of staple sequences to that layout, and then assigning the appropriate complementary bases to the scaffold. Because highly repetitive scaffold sequences can be difficult to synthesize, we used 10,000 iterations to create a library of custom scaffold sequences. We ranked the sequences according to total fraction of nucleotides that appear in a “repetitive” motif, defined as a 12 base window that appears morethan once in the scaffold. Then we selected the least repetitive scaffold sequence for synthesis.

####Userguide:

This custom python script consists basically of four parts (a more detailed description about each module can be found in the four scripts themselves):

1. **random_staple_sequence_generator.py**: This script calculates random staple (oligonucleotide) sequences if desired. In addition to generating these sequences this part also has a 'cleanup' function where sequences can be specified that should be avoided. As of now these sequences are:
	A.  aaaaa = 'AAAAA'
	B.  ccccc = 'CCCCC'
	C.  ggggg = 'GGGGG'
	D.  ttttt = 'TTTTT'
	E.  SacII = 'CCGCGG'
	F.  SalI = 'GTCGAC'
	G.  SacI = 'GAGCTC'
    H.  BtsCI = 'GGATG'
    I.  BtsCIComplement = 'CCTAC'
2. 	**Staple_breaking_and_sequence_analysis.py**: This script supervises the staple breaking (it calls staple_combinations.py) and finally analyses the repetitiveness of the final scaffold sequences. Therefore it uses an external script (see external folder (sequence was obtained from [https://code.google.com/p/py-rstr-max/](https://code.google.com/p/py-rstr-max/))). 
3. 	**staple_combinations.py**: This script generates the various staple combinations for Staple_breaking_and_sequence_analysis.py
4. 	**workflow10.py**: This is the main script that calls all other scripts. This is the only one that needs to be modified for application to custom design. A description of how that is done is provided below.

#####1. How to use the script?

...

#####2. Which parameters can be chosen? 

...

#####3. How to modify these input parameters? 

...

#####4. What is the input (requirement for cadnano/json file)?

...

#####5. What is the output? 

...


------------------

Below (Figure 1), we are showing an example result where we used 10 different staple sequences and:

1. cadnanoFile = "Example_json-files/24helix.json"
2. numTests = 1000 (= iterations)
3. Random staple sequences = True
4. minRepeatLength = 12


| [![](Example_json-files/10Lcd-24helix-random-1000-iter-56-77.png?raw=true =580x)](Example_json-files/10Lcd-24helix-random-1000-iter-56-77.png) |
|:-:|
|  Figure 1: Example output - Degree of repetitiveness for scaffold sequences for each of the 1000 designs.|

----------------

####Contact:

For questions please contact: [stefan.niekamp@ucsf.edu](stefan.niekamp@ucsf.edu)
	

	
    



