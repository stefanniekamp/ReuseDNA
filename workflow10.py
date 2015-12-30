from __future__ import division
import pickle
import json
import os
import matplotlib.pyplot as plt
import Staple_breaking_and_sequence_analysis as main
import random_staple_sequence_generator as random_staple_sequences

outputFolder = "cadnano_and_sequence-files/"
if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
    
tempicklefolder = "temppickle/"
if not os.path.exists(tempicklefolder):
    os.makedirs(tempicklefolder)

picklefile = "temppickle/temppickle.p"


"""------------------------------------------------------------------------------------"""
#########################################################################
###                                                                   ###
###        """This is the area where you can make changes."""         ###
###                                                                   ###
#########################################################################

"""Here the input location for the cadnano / json-file can be specified. Default is 24 helix bundle.""" 
cadnanoFile = "Example_json-files/24helix.json"


"""Here the number of iterations / number of different versions that will be generated can be specified. Default is 1000."""
numTests = 1000


"""Here the number of best design(s) in terms of degree of repetitivness in scaffold sequence 
ranked from lowest to highest can be specified (say you chose 1000 iterations then 
you might want to know which 10 of these 1000 designs/sequences are the best and what 
their degree of repetitivness is). NOTE: This number always has to be smaller or equal to
the number of iterations!"""
numberofbestones = 10


"""Select if you want to use random staple sequences (=True) or if you prefer predefined sequences (=False).
These predefined sequences can be specified below."""
use_random_staples = True


"""The minimum repeat length for repetitive motifs of the scaffold sequence can be chosen here
(default is 12 and we do not recommend to change it (at least not to go lower than 12))."""
minRepeatLength = 12  


"""Here the staple length can be specified below. NOTE: Every length will be used twice 
(will have two different sequences) for the design with 10 unique staple sequences."""
StapleLength01 = 56
StapleLength02 = 63
StapleLength03 = 70
StapleLength04 = 77
PolyTStapleLength = 38


"""Here the sequences can be specified below (only if you picked 
'use_random_staples = False' above. NOTE: Every staple length will be used twice 
(will have two different sequences). Thus, Seq01 and Seq02 are the sequences for StapleLength01,
Seq03 and Seq04 are the sequences for StapleLength02, ..."""
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


"""Do not make any changes below here!"""
"""------------------------------------------------------------------------------------"""


staples = (StapleLength01, StapleLength02, StapleLength03, StapleLength04)
stapleVersions = ((1.1, 1.2), (2.1, 2.2), (3.1, 3.2), (4.1, 4.2))
stapleVersionColors = ((0xb72525, 0xDB4D4D), # red  #e88d8d  #db4ddb
					   (0xcdcd00, 0xffff1a), # yellow #e6e600 #EFFE81
					   (0x006B24, 0x67cd00), # green #003813 #8dff1a
					   (0x1a1aff, 0x6767ff)) # blue #0000cd #000081

polyTStaples = (PolyTStapleLength,)
polyTStapleVersions = ((38.1, 38.2),)
polyTStapleVersionColors = ((0xFF9900, 0xb36b00),) #orange #ffb84d #673d00

if use_random_staples:
	stapleVersionSeqs = main.generateRandomSeqs(staples, stapleVersions)
	polyTStapleVersionSeqs = main.generateRandomSeqs(polyTStaples, polyTStapleVersions)
else:
	stapleVersionSeqs = ((Seq01,
						  Seq02),
						 (Seq03,
						  Seq04),
						 (Seq05,
						  Seq06),
						 (Seq07,
						  Seq08))
	polyTStapleVersionSeqs = ((Seq09,
						 Seq10),)                  
				  


step = 7

if numberofbestones > numTests:
	print "The number of iterations for different designs has to be equal to or larger than the number of best versions. Try again."
	quit()

(permutations, allStartBases, sequences, repeatData) = main.main(cadnanoFile, step,
										staples, stapleVersions, stapleVersionSeqs,
										polyTStaples, polyTStapleVersions, polyTStapleVersionSeqs,
										numTests, minRepeatLength)
										

with open(picklefile, "w") as pickleoutputfile:
	pickle.dump([permutations, allStartBases, sequences, stapleVersionSeqs, polyTStapleVersionSeqs], pickleoutputfile)

    	
    	
    	
# saving files and calculating repetitivness    	

repetetivness = []
repetetivnessnonsorted = []

for i in range(len(permutations)):

	outputjson = outputFolder + "cadnano" + str(i) + ".json"

	sequence = sequences[i] 
	permutation = permutations[i]
	startBases = allStartBases[i]
	
	# saving sequence file	
	outputsequence = outputFolder + "sequence" + str(i) + ".txt"
	text_file = open(outputsequence, "w")
	text_file.write(sequence)
	text_file.close()
	
	sequence = list(sequence)
	
	# calculating repetitivness	
	for j in range(len(repeatData[i])):
		for k in range(len(repeatData[i][j]['indices'])):
			position = repeatData[i][j]['indices'][k]
			length = repeatData[i][j]['length']
			for l in range(length):
				sequence[position + l] = 'x'
	countrepbase = sequence.count('x')	
	sequencelength = len(sequence)
	repetetivness.append(float(countrepbase / sequencelength))
	repetetivnessnonsorted.append(float(countrepbase / sequencelength))

	# saving cadnano
	with open(cadnanoFile, "r") as cadnanoInput:
		cadnanoJson = json.load(cadnanoInput)
	allHelices = cadnanoJson["vstrands"]
	brokenjson = main.breakStaples(permutation, startBases, staples, stapleVersions, stapleVersionColors,
							 polyTStaples, polyTStapleVersions, polyTStapleVersionColors, allHelices)
	with open(outputjson, "w") as brokenjsonoutput:
			cadnanojson = {"vstrands": brokenjson}
			json.dump(cadnanojson, brokenjsonoutput)

# plot repetitiveness
repetetivness.sort()
lowestrepnumber = []
lowestrepnumberrepetitivness = []
for i in range(numberofbestones):
	lowestrepnumber.append(repetetivnessnonsorted.index(repetetivness[i]))
	lowestrepnumberrepetitivness.append(repetetivness[i])
print "Best design(s) in terms of degree of repetitivness in scaffold sequence: " + str(lowestrepnumber).strip('[]') + " - ranked from lowest to highest with the following degree(s) of repetitivness: " + str(lowestrepnumberrepetitivness).strip('[]')

plt.plot(repetetivnessnonsorted, 'ro')
plt.xlabel('Number of design', fontsize=12)
plt.ylabel('Degree of repetitiveness in scaffold sequence', fontsize=12)
plt.show()                      