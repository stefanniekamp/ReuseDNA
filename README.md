##The MIT License

Copyright (c) 2015, Regents of the University of California

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[http://www.opensource.org/licenses/mit-license.php](http://www.opensource.org/licenses/mit-license.php)

##Overview

###Extract from publication:

Using [caDNAno](https://github.com/sdouglas/cadnano2), a computer-aided design tool for DNA origami, we first routed the scaffold to approximate a 3D shape and exported the design from caDNAno as a text file in JSON format. Next, we used a custom Python script (provided here) and input the JSON file and desired number of unique staples (e.g. 10). We included the option to specify specific staple sequences as input, or to input a range of desired staple lengths and automatically generate sequences. The script determines each custom scaffold sequence by generating a random staple layout, repeatedly assigning a set of staple sequences to that layout, and then assigning the appropriate complementary bases to the scaffold. Because highly repetitive scaffold sequences can be difficult to synthesize, we used 10,000 iterations to create a library of custom scaffold sequences. We ranked the sequences according to total fraction of nucleotides that appear in a “repetitive” motif, defined as a 12 base window that appears morethan once in the scaffold. Then we selected the least repetitive scaffold sequence for synthesis.

###Userguide:

This custom Python script consists basically of four parts (a more detailed description about each module can be found in the four scripts themselves):

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

####How to use the script?

The first thing you need to do is to download the folder from github. You will also need Python installed on your computer in order to run the scripts. 

Then you can modify the script so that it runs with your specific settings. Here, you only need to make shanges in the workflow10.py,  workflow15.py, or workflow20.py script. 

#####1. Which parameters can be chosen? 

These three scripts represent the different numbers of staples to use for a design. The parameters that can be changed besides the number of staples are:

1. The input cadnano / json-file
2. The number of iterations that should be used (= number of different versions that will be generated and you can chose from)
3. The number of best design(s) in terms of degree of repetitivness in scaffold sequence ranked from lowest to highest that should be shown (say you chose 1000 iterations then you might want to know which 10 of these 1000 designs/sequences are the best and what their degree of repetitivness is)
4. If you want to use predefined staple sequences or random sequences
5. The minimum repeat length for repetitive motifs (default is 12 and we do not recommend to change it (at least not to go lower than 12))
6. Staple lengths and color (Note: The staple length has to be a multiple of 7. More details below in "requirement for cadnano / json-file")

#####2. How to modify these input parameters?  

Open the workflow document in your favorite python editor and modify parameters. A detailed desciption of how to change these parameters in the script can be found in the file.

It will look like this 

	#########################################################################
	###                                                                   ###
	###        """This is the area where you can make changes."""         ###
	###                                                                   ###
	#########################################################################

	"""Here the input location for the cadnano / json-file can be specified. Default is 24 helix bundle.""" 
	cadnanoFile = "Example_json-files/24helix.json"


	"""Here the number of iterations / number of different versions that will be generated can be specified. Default is 1000."""
	numTests = 1000


	"""Here the number of best design(s) in terms of degree of repetitivness in scaffold sequence ranked from lowest to highest can be specified (say you chose 1000 iterations then you might want to know which 10 of these 1000 designs/sequences are the best 	and what their degree of repetitivness is). NOTE: This number always has to be smaller or equal to the number of iterations!"""
	numberofbestones = 10


	"""Select if you want to use random staple sequences (=True) or if you prefer predefined sequences (=False). These predefined sequences can be specified below."""
	use_random_staples = True


	"""The minimum repeat length for repetitive motifs of the scaffold sequence can be chosen here (default is 12 and we do not recommend to change it)."""
	minRepeatLength = 12
	

	"""Here the staple length can be specified below. NOTE: Every length will be used twice (will have two different sequences) for the design with 10 unique staple sequences."""
	StapleLength01 = 56
	StapleLength02 = 63
	StapleLength03 = 70
	StapleLength04 = 77
	PolyTStapleLength = 38


	"""Here the sequences can be specified below (only if you picked 'use_random_staples = False' above. NOTE: Every staple length will be used twice (will have two different sequences). Thus, Seq01 and Seq02 are the sequences for StapleLength01, Seq03 and Seq04 are the sequences for StapleLength02, ..."""
	Seq01 = "GAGTTTTACGTCTAGTCTCCGCTACAAATGGAGTCACGAAATAGGGCACCATCGTC"
	Seq02 = "TAGCCTCTCAGTGTAGTTAAGATTTAGGAGTTCGCAACTGTGAGGACTTCGTGCGA"
	Seq03 = "TACGTGCATTCGCTTTGGAGGCATTCTCGCTTCCAAACCATCGATGTTTACGTAGGCGCTGTT"
	Seq04 = "GACAAATATCTTCTGCACAAATCCCGTCAGAGAGCCGCGTGTACTGGATTTATCGGCCGACAT"
	Seq05 = "CATAATACTCATATGTGATGCTCGAAACTGCTGAACGGTGTTAACTGCTATGAAGACCATAAGTCATGAC"
	Seq06 = "CGCGCACCCATCCGCCCTATTGAAACGGGTTGTTGCGAAGCGTAAGGAGCACAGCGAGGGGCGGGAGCGC"
	Seq07 = "ATACAAGCAATCCACGCCGACCGGCCGATCGAAAGGACGGTCATATACCCGTATTGTCCTGTTAGTCAAACTGGGAC"
	Seq08 = "ACAATCCACGGCAAATACTCCTGATGATCATATGCACGGTCTCCTTCGCTCGCAGGCCTCAACAACCGGCCATACTG"
	Seq09 = "AGCATACGTACCCTGATCCCAGTGTAGATATACAGAAT"
	Seq10 = "AACAGCTGGCCATTGCAGGGTATGCCCATAGACGCGAA"

#####3. What is the input (requirement for cadnano/json file)?

As of now, the script only works for the honeycomb lattices. For the cadnano design, there are a few requirements: The first thing is that all staples have to form loops and should not be broken. In addition, we found that the highest folding yield is achieved if you use a high staple crossover density and a lower scaffold crossover density. Thus, use all available staple crossovers when forming loops. For the loop size, larger loops are usually easier to solve and allow for a higher flexibility for staple breaking by the Python script. If loops are too small, it might happen, that there is only a small or a too small solution space, which will either lead to high degree of repetitiveness in scaffold sequences (small space) or cause the script to crash (too small space) since it will not be able to solve the problem. An example for a not successful calculation would be a loop size of 42 bp with a shortest staple length of 56 bp (ignoring the “end piece” staples here”; more about them below). Below you can see what we mean by staple loops.

| [![](Example_json-files/24helixbundle_small-loop-size.png?raw=true =580x)](Example_json-files/24helixbundle_small-loop-size.png) |
|:-:|
|  **Figure 1:** Example design for a 24 helix bundle with small loop size (red, 504 bp).|

| [![](Example_json-files/24helixbundle_medium-loop-size.png?raw=true =580x)](Example_json-files/24helixbundle_medium-loop-size.png) |
|:-:|
|  **Figure 2:** Example design for a 24 helix bundle with medium loop size (red, 2016 bp).| 

| [![](Example_json-files/24helixbundle_large-loop-size.png?raw=true =580x)](Example_json-files/24helixbundle_large-loop-size.png) |
|:-:|
|  **Figure 3:** Example design for a 24 helix bundle with large loop size (red, 5040 bp).|    

Second, for each design you will have "end pieces", meaning the position where a helix ends. For these end pieces you can have different lengths. For the cadnano honeycomb lattice design you will usually have 16 bp, 38 bp, and 58 bp long staples if you use 3 bp on each side of the staple for polyT passivation. We highly recommend to use 38 bp long "end pieces”, which would look like this:

| [![](Example_json-files/end-pieces_polyTpassivation_38bp.png?raw=true =280x)](Example_json-files/end-pieces_polyTpassivation_38bp.png) |
|:-:|
|  **Figure 4:** Example "end pieces" for cadnano honeycomb lattice design with polyT passivation (length = 38 bp).|  

If you decide to use "end piece" staples of a different lengths make sure to design them in a similar fashion as we did. This means, use one staple to connect two helices via crossover but never connect more than two helices with "end piece" staples (would cause script to crash). But using "end piece" staples, which are not connected (no crossover; only bound to one helix), works just fine. If you change the length in the cadnano (json) file make sure to adjust the length for polyTstaples in the Python script accordingly.  

#####4. What is the output? 

Below in Figure 5 we are showing an example result where we used 10 unique staple sequences and:

1. cadnanoFile = "Example_json-files/24helix.json"
2. numTests = 1000 (= iterations)
3. numberofbestones = 10
4. Random staple sequences = True
5. minRepeatLength = 12
6. Staple lengths (each twice): 56, 63, 70, 77, 38

| [![](Example_json-files/10Lcd-24helix-random-1000-iter-56-77.png?raw=true =580x)](Example_json-files/10Lcd-24helix-random-1000-iter-56-77.png) |
|:-:|
|  **Figure 5:** Example output for design with 10 unique staple sequences - Degree of repetitiveness for scaffold sequences for each of the 1000 designs.|

In addition, the output will be the following text in your terminal: "*Best design(s) in terms of degree of repetitivness in scaffold sequence: 302, 134, 439, 789, 211, 54, 589, 690, 559, 388 - ranked from lowest to highest with the following degree(s) of repetitivness: 0.44455922865013775, 0.46143250688705234, 0.4712465564738292, 0.4758953168044077, 0.48364325068870523, 0.4839876033057851, 0.4852827364518031, 0.48674242424242425, 0.48794765840220383, 0.48829201101928377*"  and the json as well as sequence files for all 1000 versiosn are saved in the "cadnano_and_sequence-files" folder. In Figure 6 you can see an example cadnano file (here version 302). The color scheme is as follows:

A. Shades of red (56 bp long staples)
B. Shades of yellow (63 bp long staples)
C. Shades of green (70 bp long staples)
D. Shades of blue (77 bp long staples)
E. Shades of orange (38 bp long staples)

| [![](Example_json-files/Example-output-for-24helixbundle_cadnano-design.png?raw=true =580x)](Example_json-files/Example-output-for-24helixbundle_cadnano-design.png) |
|:-:|
|  **Figure 6:** Example (version 302) output for design with 10 unique staple sequences - Cadnano plot.|

We recommend to look at the 10 best versions in cadnano and then pick one.

Just for comparison, in Figure 7 you can see the result where we used 15 unique staple sequences and:

1. cadnanoFile = "Example_json-files/24helix.json"
2. numTests = 1000 (= iterations)
3. numberofbestones = 10
4. Random staple sequences = True
5. minRepeatLength = 12
6. Staple lengths (each twice): 56, 63, 70, 77, 38

| [![](Example_json-files/15Lcd-24helix-random-1000-iter-56-77.png?raw=true =580x)](Example_json-files/15Lcd-24helix-random-1000-iter-56-77.png) |
|:-:|
|  **Figure 7:** Example output for design with 15 unique staple sequences - Degree of repetitiveness for scaffold sequences for each of the 1000 designs.|

And in Figure 8 you can see the result where we used 20 unique staple sequences and:

1. cadnanoFile = "Example_json-files/24helix.json"
2. numTests = 1000 (= iterations)
3. numberofbestones = 10
4. Random staple sequences = True
5. minRepeatLength = 12
6. Staple lengths (each twice): 56, 63, 70, 77, 38

| [![](Example_json-files/20Lcd-24helix-random-1000-iter-56-77.png?raw=true =580x)](Example_json-files/20Lcd-24helix-random-1000-iter-56-77.png) |
|:-:|
|  **Figure 8:** Example output for design with 20 unique staple sequences - Degree of repetitiveness for scaffold sequences for each of the 1000 designs.|

#####5. Run the script

Just open a terminal and go to the directory where you saved the folder with all the scripts. Then type:

	python workflow10.py
	
Click enter. Depending on the cadnano file and the number of iterations as well as your computational power the script might run for a while. On my Mac with a 2.5 GHz processor and 8GB 1600 MHz DDR3 Memory 1000 iterations for the 24 helix bundle takes about an hour.  

If you just want to see if the script is running properly use the example cadnano / json-files which are provided in the "Example_json-files" folder.


----------------

####Contact:

For questions please [contact the autors](stefan.niekamp@ucsf.edu).
	

	
    



