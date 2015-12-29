from __future__ import division
from __future__ import with_statement

import copy
import json
import math
from mpl_toolkits.mplot3d.axes3d import Axes3D
import random

import external.rstr_max as SuffixArray
import matplotlib.pyplot as plt
import numpy as np
import staple_combinations as combinations
import random_staple_sequence_generator as random_staple_sequences


"""This script supervises the staple breaking (it calls staple_combinations.py) and 
finally analyses the repetitiveness of the final scaffold sequences. Therefore it uses 
an external script (see external folder (sequence was obtained 
from [https://code.google.com/p/py-rstr-max/](https://code.google.com/p/py-rstr-max/)))."""


def main(cadnanoFile, step, staples, stapleVersions, stapleVersionSeqs,
         polyTStaples, polyTStapleVersions, polyTStapleVersionSeqs, numTests, minRepeatLength):

##### IMPORT CADNANO SCAFFOLD/STAPLE ROUTING
    with open(cadnanoFile) as data_file:    
        unbrokenCadnanoJson = json.load(data_file)
    #end with
    
    allHelices = unbrokenCadnanoJson["vstrands"]
    
#### FIND STRANDS
    strandsInfo = strandLocationsAndLengths(allHelices)
    strandLengths = strandsInfo["lengths"]  

##### PERMUTE EACH STRAND
    permutations = []
    allStartBases = []
    sequences = []
    repeatData = []
    
    successfulTests = 0
    failedTests = 0
    
    while successfulTests < numTests:
        permutation = combinations.randomSolution(strandLengths, staples, stapleVersions,
                                                  polyTStaples, polyTStapleVersions)
        randomStartBases = randomStartLocations(strandsInfo, allHelices, step)

##### FIND SCAFFOLD SEQUENCE FOR EACH PERMUTATION
        scaffoldSequence = makeScaffoldSequence(permutation, randomStartBases,
                                                stapleVersions, stapleVersionSeqs,
                                                polyTStapleVersions, polyTStapleVersionSeqs, allHelices)
        
##### SAVE PERMUTATION INFO
        if random_staple_sequences.checkScaffoldCorrectness(scaffoldSequence):
            successfulTests += 1
            permutations.append(permutation)
            allStartBases.append(randomStartBases)
            sequences.append(scaffoldSequence)
            repeatData.append(findDistances(scaffoldSequence, minRepeatLength, allHelices, step))
        else:
            failedTests += 1
        #end if
	
    #end for
    
    
    

##### SAVE/OUTPUT DATA
    return (permutations, allStartBases, sequences, repeatData)
#end def


'''
Traverses all of the staple strands and outputs their lengths and their beginning positions.
For broken strands, the beginning position is the 5' end.
For looped strands, the beginning position is chosen at random at a point that is suitable for the 5' end of a break
'''
def strandLocationsAndLengths(allHelices):
    ## This function walks along the staple strand which contains the given base pair. It finds the
    ## length, going backwards as well if necessary (if staple is not loop). 
    def walkOneStrand(startHelix, startPos, basesSeen, allHelices):
        
        # will be used to track which bases have already been traversed
        
        currHelix = startHelix
        currPos = startPos
        length = 0
        
        # walk forwards along the scaffold till either you hit the end or end up where you started.
        # if you hit the end, go backwards to find the beginning. Otherwise it's a loop; find a
        # crossover (so you know a good breaking point). 
        # When done, you have your start point and your strand length.
        while True:
            if (basesSeen[currHelix][currPos] == True):
                # is loop and we're back at the beginning.
                # Find a good breaking point and start loop right after it
                while True:
                    (nextHelix, nextPos) = walkForward(1, currHelix, currPos, "stap", allHelices)
                    if (nextHelix != currHelix):
                        # nextHelix occurs right after a crossover
                        (begHelix, begPos) = (nextHelix, nextPos)
                        break
                    else:
                        (currHelix, currPos) = (nextHelix, nextPos)
                    #end if
                #end while
                    
                break
            elif (currHelix == -1):
                # is not loop, we go back to where we started, walk backwards, and add to length
                # till we reach the beginning
                (currHelix, currPos) = (startHelix, startPos)
                while True:
                    (prevHelix, prevPos) = walkBackward(1, currHelix, currPos, "stap", allHelices)
                    if (prevHelix != -1):
                        (currHelix, currPos) = (prevHelix, prevPos)
                        basesSeen[currHelix][currPos] = True
                        length += 1
                    else:
                        break
                    #end if
                #end while
                (begHelix, begPos) = (currHelix, currPos)
                break
            else:
                length += 1 
                basesSeen[currHelix][currPos] = True
                
                # move to next position
                (currHelix, currPos) = walkForward(1, currHelix, currPos, "stap", allHelices)
            #end if 
        #end while
        
        return ((begHelix, begPos), length)
    #end def
    
    startBases = []
    lengths = []
    
    seenBase = [[False] * len(helix["stap"]) for helix in allHelices]
    
    # go through every scaffold base; if it has a staple there and you haven't already traversed it,
    # do so.
    for helixNum in range(len(allHelices)):
        helix = allHelices[helixNum]
        for pos in range(len(helix["stap"])):
            
            if ((seenBase[helixNum][pos] == False) and (helix["stap"][pos] != [-1, -1, -1, -1])):
                (startBase, length) = walkOneStrand(helixNum, pos, seenBase, allHelices)
                startBases.append(startBase)
                lengths.append(length)
            #end if
        #end for
    #end for
    
    return {"startBases":startBases, "lengths":lengths}
#end def


def randomStartLocations(strandsInfo, allHelices, step):
    newStartBases = []
    
    startBases = strandsInfo["startBases"]
    strandLengths = strandsInfo["lengths"]
    
    for i in range(len(startBases)):
        startBase = startBases[i]
        prevBase = walkBackward(1, startBase[0], startBase[1], "stap", allHelices)
        
        if (prevBase == (-1, -1)):
            # is not a loop, don't change start base
            newStartBases.append(startBase)
        else:
            stepsForward = random.randrange(0, strandLengths[i], step)
            newBase = walkForward(stepsForward, startBase[0], startBase[1], "stap", allHelices)
            newStartBases.append(newBase)
        #end if
        
    return newStartBases
#end def
        
              
def makeScaffoldSequence(permutation, startBases, stapleVersions, stapleVersionSeqs,
                         polyTStapleVersions, polyTStapleVersionSeqs, allHelices):

    complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
    
    flatStapleVersions = flattenOneLevel(stapleVersions)
    flatStapleVersionSeqs = flattenOneLevel(stapleVersionSeqs)
    
    flatPolyTVersions = flattenOneLevel(polyTStapleVersions)
    flatPolyTVersionSeqs = flattenOneLevel(polyTStapleVersionSeqs)
    
    scafSeqs = [["?"] * len(helix["stap"]) for helix in allHelices]
    
    # copy staple sequences to scaffold
    for strand in range(len(permutation)):
        staplePermutation = permutation[strand]
        startBase = startBases[strand]
        
        currHelix = startBase[0]
        currPos = startBase[1]
        
        for stapleVersion in staplePermutation:
            # make it work for normal staples or poly-T staples
            if stapleVersion in flatStapleVersions:
                stapleVersionIndex = flatStapleVersions.index(stapleVersion)
                stapleSeqList = flatStapleVersionSeqs[stapleVersionIndex]
            else:
                stapleVersionIndex = flatPolyTVersions.index(stapleVersion)
                stapleSeqList = flatPolyTVersionSeqs[stapleVersionIndex]
            #end if
            
            for base in stapleSeqList:
                scafSeqs[currHelix][currPos] = complement[base]
                (currHelix, currPos) = walkForward(1, currHelix, currPos, "stap", allHelices)
            #end for
        #end for
    #end for
    
    #walk scaffold from start and put sequences in order
    (currHelix, currPos) = findBegOfScaffold(allHelices)
    
    scafSeqList = []
    
    scafSeqList.append(scafSeqs[currHelix][currPos])
    (currHelix, currPos) = walkForward(1, currHelix, currPos, "scaf", allHelices)
    
    while (currHelix!= -1):
        scafSeqList.append(scafSeqs[currHelix][currPos])
        (currHelix, currPos) = walkForward(1, currHelix, currPos, "scaf", allHelices)
    #end while
    
    scafSeq = ''.join(scafSeqList)

    return scafSeq
#end def


def breakStaples(permutation, startBases, staples, stapleVersions, stapleVersionColors,
                         polyTStaples, polyTStapleVersions, polyTStapleVersionColors, allHelices):
    brokenHelices = copy.deepcopy(allHelices)
    #remove color
    for helix in brokenHelices:
        helix["stap_colors"] = []
    #end for
    
    flatStapleVersions = flattenOneLevel(stapleVersions)
    flatStapleVersionLengths = flattenOneLevel([[staples[i]] * len(stapleVersions[i]) for i in range(len(staples))])
    flatStapleVersionColors = flattenOneLevel(stapleVersionColors)
    
    flatPolyTVersions = flattenOneLevel(polyTStapleVersions)
    flatPolyTVersionLengths = flattenOneLevel([[polyTStaples[i]] * len(polyTStapleVersions[i]) for i in range(len(polyTStaples))])
    flatPolyTVersionColors = flattenOneLevel(polyTStapleVersionColors)
    
    for strand in range(len(permutation)):
        staplePermutation = permutation[strand]
        startBase = startBases[strand]
        
        currHelix = startBase[0]
        currPos = startBase[1]
        
        for stapleVersion in staplePermutation:
            # make it work for normal staples or poly-T staples
            if stapleVersion in flatStapleVersions:
                stapleVersionIndex = flatStapleVersions.index(stapleVersion)
                stapleLength = flatStapleVersionLengths[stapleVersionIndex]
                stapleColor = flatStapleVersionColors[stapleVersionIndex]
            else:
                stapleVersionIndex = flatPolyTVersions.index(stapleVersion)
                stapleLength = flatPolyTVersionLengths[stapleVersionIndex]
                stapleColor = flatPolyTVersionColors[stapleVersionIndex]
            #end if
            
            #break beginning of staple
            brokenHelices[currHelix]["stap"][currPos][0:2] = [-1, -1]
            
            #add color
            brokenHelices[currHelix]["stap_colors"].append([currPos, stapleColor])
            
            #walk to end of staple
            (currHelix, currPos) = walkForward(stapleLength - 1, currHelix, currPos, "stap", allHelices)
                
            # break end of staple
            brokenHelices[currHelix]["stap"][currPos][2:4] = [-1, -1]
            
            # walk forward to start at beginning of next staple
            (currHelix, currPos) = walkForward(1, currHelix, currPos, "stap", allHelices)
        #end for
    
    #end for
    
    return brokenHelices
#end def


def generateRandomSeqs(staples, stapleVersions):
    stapleSeqs = []
    for stapleIndex in range(len(staples)):
        stapleVersionSeqs = []
        for stapleVersion in stapleVersions[stapleIndex]:
            randomSeq = random_staple_sequences.main(staples[stapleIndex])
            stapleVersionSeqs.append(randomSeq)
        #end for
        stapleSeqs.append(stapleVersionSeqs)
    #end for
    return stapleSeqs
#end def


'''
Finds all repeat sequences longer than a given length in a given string.

@return: a list of lists; for each repeated sequence, contains the tuple (length of sequence,
number of times repeated)
'''
def findRepeatSeqs(seq, minRepeatLength):
    repeatData = []
    suffixArrayBuilder = SuffixArray.Rstr_max()
    suffixArrayBuilder.add_str(seq)
    suffixArray = suffixArrayBuilder.go()
    for (dummy1, numRepeats), (repeatLength, dummy2) in suffixArray.iteritems():
        if (repeatLength >= minRepeatLength):
            repeatData.append((repeatLength, numRepeats),)
        #end if
    #end for
    return repeatData
#end def


def findDistances(seq, minRepeatLength, allHelices, step):
    
    suffixArray = SuffixArray.Rstr_max()
    suffixArray.add_str(seq)
    suffixArrayList = suffixArray.go()
    
    repeatSeqs = []
    
    for (uniqueSeqEndIndexInStr, numRepeats), (repeatLength, startIndexInArray) in suffixArrayList.iteritems():
        if (repeatLength >= minRepeatLength):
            repeatData = {"length":repeatLength, "indices":[], "positions":[], "distances":[]}
            
            for arrayIndex in range(startIndexInArray, startIndexInArray + numRepeats):
                offset_global = suffixArray.res[arrayIndex]
                startIndexInStr = suffixArray.idxPos[offset_global]
                
                # measure position from middle (or just past middle, for an even-length repeat)
                middleIndex = startIndexInStr + (repeatLength // 2)
                
                pos3D = find3DPosOfScaffoldBase(middleIndex, allHelices, step)
                
                # add info to repeatData
                repeatData["indices"].append(startIndexInStr)
                repeatData["positions"].append(pos3D)
                
                iDistancesList = []
                i = arrayIndex - startIndexInArray
                for j in range(i):
                    ijDistance = find3DDistance(pos3D, repeatData["positions"][j])
                    repeatData["distances"][j].append(ijDistance)
                    iDistancesList.append(ijDistance)
                #end for
                iDistancesList.append(0)
                repeatData["distances"].append(iDistancesList)
            #end for
            
            repeatSeqs.append(repeatData)
        #end if
    #end for
    return repeatSeqs      
#end def

'''
Returns a tuple representing the x, y, and z positions of the base at the specified index in the
scaffold sequence.
Note that y is measured from the top of the caDNAno helix diagram. X is measured from the far left,
and z from position 0 on the helix. 
'''
def find3DPosOfScaffoldBase(seqIndex, allHelices, step):
    (currHelix, currPos) = findBegOfScaffold(allHelices)
    (baseHelix, basePos) = walkForward(seqIndex, currHelix, currPos, "scaf", allHelices)
    
    row = allHelices[baseHelix]["row"]
    col = allHelices[baseHelix]["col"]
    
    return find3DPosOfBase(row, col, basePos, step)


def find3DPosOfBase(row, col, base, step):   
    bpWidth = 0.33 # nm
    helixDiam = 2 # nm
    
    if (step == 8): # is rectangular
        x = helixDiam * col
        y = helixDiam * row
    elif (step == 7): # is hexagonal
        # helix width * (3/2) * (sqrt(3)/2) * 2 / 3
        # = helix width * sqrt(3) / 2
        # = col width
        hexColWidth = helixDiam * math.sqrt(3) / 2
        
        x = hexColWidth * col
        
        helicesBefore = row
        if (col // 2 == 0): # is even column
            emptyHelicesBefore = (row + 1) // 2
        elif (col // 2 == 1): # is odd column
            emptyHelicesBefore = (row // 2) + 0.5
        #end if
        
        y = helixDiam * (row + helicesBefore)
    # end if
    
    z = bpWidth * base
    
    return (x, y, z)
#end def


def find3DDistance(pos1, pos2):
    (x1, y1, z1) = pos1
    (x2, y2, z2) = pos2
    (dx, dy, dz) = (x2-x1, y2-y1, z2-z1)
    return math.sqrt(dx**2 + dy**2 + dz**2)
#end def


def findMaxDistance(allHelices, step):
    # this just draws a box around the caDNAno diagram and measures distance from opposite corners, TODO improve
    maxRow = -1
    minRow = float("inf")
    maxCol = -1
    minCol = float("inf")
    minBase = 0
    maxBase = len(allHelices[0]["scaf"])
    for helixNum in range(len(allHelices)):
        row = allHelices[helixNum]["row"]
        col = allHelices[helixNum]["col"]
        if row > maxRow:
            maxRow = row
        if row < minRow:
            minRow = row
        if col > maxCol:
            maxCol = col
        if col < minCol:
            minCol = col
        #end if
    #end for
        
    minPos = find3DPosOfBase(minRow, minCol, minBase, step)
    maxPos = find3DPosOfBase(maxRow, maxCol, maxBase, step)
    
    return find3DDistance(minPos, maxPos)
#end def
        

def saveAsJson(allHelices, filename):
    toSave = {}
    toSave["vstrands"] = allHelices
    with open(filename, 'w') as outFile:
        json.dump(toSave, outFile)
    #end with
#end def


'''
Utility functions to move forward and backward along a strand.
'''
def walkBackward(numSteps, currHelix, currPos, scafOrStap, allHelices):
    for i in range(numSteps):
        adjacentBases = allHelices[currHelix][scafOrStap][currPos]
        currHelix = adjacentBases[0]
        currPos = adjacentBases[1]
    return (currHelix, currPos)
#end def


def walkForward(numSteps, currHelix, currPos, scafOrStap, allHelices):
    for i in range(numSteps):
        adjacentBases = allHelices[currHelix][scafOrStap][currPos]
        currHelix = adjacentBases[2]
        currPos = adjacentBases[3]
    return (currHelix, currPos)
#end def


def findBegOfScaffold(allHelices):
    for helixIndex in range(len(allHelices)):
        helix  = allHelices[helixIndex]
        for i in range(len(helix["scaf"])):
            scafBase = walkForward(1, helixIndex, i, "scaf", allHelices)
            if (scafBase != (-1, -1)):
                break
            #end if
        #end for
        if (scafBase != (-1, -1)):
            break
        #end if
    #end for
    
    scafStartBase = findBegOfStrand(scafBase[0], scafBase[1], "scaf", allHelices)
    
    return scafStartBase
#end def


def findBegOfStrand(startHelix, startPos, ScafOrStap, allHelices):
    (currHelix, currPos) = (startHelix, startPos)
    (nextHelix, nextPos) = walkBackward(1, currHelix, currPos, ScafOrStap, allHelices)
    while True:
        if (nextHelix == -1):
            return (currHelix, currPos)
        elif ((nextHelix, nextPos) == (startHelix, startPos)):
            return None
        else:
            (currHelix, currPos) = (nextHelix, nextPos)
            (nextHelix, nextPos) = walkBackward(1, currHelix, currPos, ScafOrStap, allHelices)
        #end if
    #end while
#end def


def findLengthOfSubStrand(startBase, endBase, scafOrStap, allHelices):
    currBase = startBase
    length = 1
    while not (currBase == endBase):
        currBase = walkForward(1, currBase[0], currBase[1], scafOrStap, allHelices)
        length += 1
    #end while
    return length
#end def
    
'''
Utility function to remove (one level of) nesting in lists while maintaining order.
'''
def flattenOneLevel(origList):
    return [item for sublist in origList for item in sublist]
#end def