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
numberofbestones = 3


"""The minimum repeat length for repetitive motifs of the scaffold sequence can be chosen here
(default is 12 and we do not recommend to change it)."""
minRepeatLength = 12


staples = (56, 63, 70, 77)
stapleVersions = ((56.1, 56.2), (63.1, 63.2), (70.1, 70.2), (77.1, 77.2))
stapleVersionColors = ((0x8F0000, 0xDB4D4D), # red
					   (0x999900, 0xE6E600), # yellow
					   (0x006B24, 0x33AD5C), # green
					   (0x660080, 0xD633FF)) # purple

polyTStaples = (38,)
polyTStapleVersions = ((38.1, 38.2),)
polyTStapleVersionColors = ((0xFF9900, 0x005CE6),) #orange and blue



use_random_staples = True

if use_random_staples:
	stapleVersionSeqs = main.generateRandomSeqs(staples, stapleVersions)
	polyTStapleVersionSeqs = main.generateRandomSeqs(polyTStaples, polyTStapleVersions)
else:
	stapleVersionSeqs = (("GAGTTTTACGTCTAGTCTCCGCTACAAATGGAGTCACGAAATAGGGCACCATCGTC",
						  "TAGCCTCTCAGTGTAGTTAAGATTTAGGAGTTCGCAACTGTGAGGACTTCGTGCGA"),
						 ("TACGTGCATTCGCTTTGGAGGCATTCTCGCTTCCAAACCATCGATGTTTACGTAGGCGCTGTT",
						  "GACAAATATCTTCTGCACAAATCCCGTCAGAGAGCCGCGTGTACTGGATTTATCGGCCGACAT"),
						 ("CATAATACTCATATGTGATGCTCGAAACTGCTGAACGGTGTTAACTGCTATGAAGACCATAAGTCATGAC",
						  "CGCGCACCCATCCGCCCTATTGAAACGGGTTGTTGCGAAGCGTAAGGAGCACAGCGAGGGGCGGGAGCGC"),
						 ("ATACAAGCAATCCACGCCGACCGGCCGATCGAAAGGACGGTCATATACCCGTATTGTCCTGTTAGTCAAACTGGGAC",
						  "ACAATCCACGGCAAATACTCCTGATGATCATATGCACGGTCTCCTTCGCTCGCAGGCCTCAACAACCGGCCATACTG"))
	polyTStapleVersionSeqs = (("AGCATACGTACCCTGATCCCAGTGTAGATATACAGAAT",
						 "AACAGCTGGCCATTGCAGGGTATGCCCATAGACGCGAA"),)                  
				  

"""Do not make any changes below here!"""
"""------------------------------------------------------------------------------------"""

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