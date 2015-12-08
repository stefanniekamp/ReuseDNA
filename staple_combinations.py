from __future__ import division
import bisect
import itertools
import math
import random
import timeit

'''
General order of things:

"Strand": a staple strand that has not yet been broken - it can be a loop or have endpoints.

"Staple": a strand that has been broken into one of the allowed lengths.

"Staple Version": a staple that has been assigned to one of the specific sequences given by the user;
    also, the unique numeric identifier for a specific sequence.

"Partition": a list of staple lengths which adds up to the length of a strand;
    one of (usually) many ways in which a particular strand can be broken up.
    Always sorted by staple length.

"Partition Version": a partition in which each of the staple lengths has been assigned to a
    specific sequence of that length, in the form of a list of numeric sequence identifiers.
    Always in sorted order (first by length, then by version).
    
External Variables:

"strandLengths": a list of the lengths of all the strands in the given design.
    
"staples": list of possible staple lengths. Currently passed by user.
    Always sorted.

"stapleVersions": two-dimensional list; each inner list corresponds to a staple length, and is a
    list of all versions of that length.
    Always sorted, outer by length, inner by version.
    
"polyTStaples/StapleVersions": same as staples/stapleVersions, but reserved for the strands
    with poly-Ts. If the design is done as expected, they should be all the same length and shorter
    than the normal staples.
'''


'''
Given a list of all the staple strand lengths in a design and available staples, randomly chooses
staples to break up the design.

@return: list of lists; for every strand in strandLengths, a list of staple versions, in the order
they should be broken, which add up to that strand length.
'''
def randomSolution(strandLengths, staples, stapleVersions, polyTStaples, polyTStapleVersions):
    # improve randomness
    random.seed(timeit.default_timer())
    
    solutionList = [[] for i in range(len(strandLengths))]
    for strandLength in set(strandLengths):
        if (strandLength < min(staples)):
            # is a poly-T strand
            solutionsForLength = makeStapleList(strandLengths.count(strandLength), strandLength,
                                            polyTStaples, polyTStapleVersions)
        else:
            # is a normal strand
            solutionsForLength = makeStapleList(strandLengths.count(strandLength), strandLength,
                                            staples, stapleVersions)
        #end if
        
        currentSolution = 0
        for i in range(len(strandLengths)):
            if (strandLengths[i] == strandLength):
                solutionList[i] = solutionsForLength[currentSolution]
                currentSolution+=1
            #end if
        #end for
    #end for
    
    return solutionList
#end def


'''
Given the length of a strand, returns a list of (multiple, if desired) random breakings (partitions)
into available staple versions.

NOTE ON RANDOMNESS OF SAMPLING: The random partition and version is chosen at random from all
possible versions of all possible partitions; however, rather than also choosing from all possible
orderings of all versions, the version is simply chosen and then shuffled.
This is because for looped strands, we later randomly rotate the chosen solution so we're not always
breaking at the top left corner of a loop. Rotating would result in duplicates for some
partition versions. We could solve this by finding all possible *necklaces* rather than
orderings; however, there is (as far as I can tell) no way to find the number of fixed-content
necklaces for a given partition version without explicitly calculating all of them. That would 
require more memory than my computer has. If you've read this far, you may like this paper on
fixed-content necklace generation:
http://www.sciencedirect.com/science/article/pii/S0304397503000495

@param numToChoose: Number of different permutations to generate.
@param strandLength: Number of base pairs in the loop.
@param staples: Sorted list of possible staple lengths.
@param stapleVersions: List of lists, one for each staple; each list contains unique numerical
identifiers for the different possible staples of that length - they can be floats, must be sorted.
@param step: Length between possible crossovers (7 for hexagonal structures, 8 for rectangular)

@return: List containing the requested number of solutions for the loop, each solution containing a
random staple partition, random staple versions, in random order.
'''
def makeStapleList(numToChoose, strandLength, staples, stapleVersions):
    
    partitions = staplePartitions(strandLength, staples)
    if len(partitions) < 1:
        print "Could not partition staple with length {}".format(strandLength)
    #end if
    partitionVersionCounts = partitionVersionCountsForStrand(partitions, staples, stapleVersions)
    totalSolutions = partitionVersionCounts[-1]
    chosenSolutions = []
    
    for i in range(numToChoose):
        partitionVersionToChoose = random.randrange(0, totalSolutions)
        partitionIndex = bisect.bisect_right(partitionVersionCounts, partitionVersionToChoose)
        
        solutionsBeforePartition = partitionVersionCounts[partitionIndex - 1] \
                                   if partitionIndex > 0 else 0
        partitionVersionNum = partitionVersionToChoose - solutionsBeforePartition
        partitionVersions = staplePartitionVersions(partitions[partitionIndex],
                                                    staples, stapleVersions)
        
        shuffledVersion = partitionVersions[partitionVersionNum][:]
        random.shuffle(shuffledVersion)
         
        chosenSolutions.append(shuffledVersion)
    # end for
    
    return chosenSolutions
#end def


'''
Finds all of the different ways a strand can be broken into staples of the given lengths.

@param strandLength: Length of the strand to be broken.
@param staples: List of the possible staple lengths Must be sorted and cannot contain duplicates.

@return: List of lists; each inner list is a sorted list of staple lengths which add up to loopSize,
and the outer list is sorted in lexicographical order.
'''
def staplePartitions(strandLength, staples):


    allSubsets = [] # will be list of sorted partitions
    smallest = staples[0]
    
    # (base case) check whether the sum can be completed with just the smallest integer
    # (if it can, no point going further)
    if (smallest == strandLength):
        allSubsets.append([smallest])
    elif smallest < strandLength:
        # find all subsets that contain at least one more copy of the smallest integer
        # (add integer to every set, subtract it from sum)
        smallerSets = staplePartitions(strandLength - smallest, staples)
        for smallerSet in smallerSets:
            smallerSet.insert(0, smallest)
        allSubsets.extend(smallerSets)
        # find all subsets that don't (remove smallest integer from set, then recurse)
        if (len(staples) > 1):
            allSubsets.extend(staplePartitions(strandLength, staples[1:]))
        # end if
    #end if
    
    return allSubsets
    
#end def


'''
Finds all possible partition versions for the given partition.

@param staplePartition: Sorted list of staple lengths.
@param staples: List of all possible staple lengths.
@param stapleVersions: List of lists, one for each staple; each list contains unique identifiers for
the different possible staples of that length.

@return: list of every possible partition version; each inner
list is sorted and the outer list is in lexicographical order.
'''
def staplePartitionVersions(staplePartition, staples, stapleVersions):
    versionCombinationsForStaple = []
    
    # finding every possible combination of staple versions for one staple at a time, given the
    # count of that staple in the partition, then saving in versionCombinationsForStaple
    for i in range(len(staples)):
        combIter = itertools.combinations_with_replacement(\
                stapleVersions[i], staplePartition.count(staples[i]))
        versionCombinationsForStaple.append([x for x in combIter])
    # end for
    
    # combining version sets so we have one set for each staple
    partitionVersions = []
    versionCombinationIndices = [0 for i in range(len(staples))]
    
    def convertVersionComboToList(versionCombinationIndices, versionCombinationsForStaple):        
        return [item for i in range(len(staples)) for item in versionCombinationsForStaple[i][versionCombinationIndices[i]]]
    # end def
    
    partitionVersions.append(convertVersionComboToList(versionCombinationIndices, versionCombinationsForStaple))
    
    stapleIndex = len(staples) - 1
    while (stapleIndex > -1):
        #increment versionCombinationIndices
        versionCombinationIndices[stapleIndex] += 1
        if (versionCombinationIndices[stapleIndex] >= len(versionCombinationsForStaple[stapleIndex])):
            # if you went too high, set digit to 0, start over at end of versionCombinationIndices
            versionCombinationIndices[stapleIndex] = 0
            stapleIndex = stapleIndex - 1
        else:
            # save combination of versions denoted by versionCombinationIndices
            partitionVersions.append(convertVersionComboToList(versionCombinationIndices, versionCombinationsForStaple))
            stapleIndex = len(staples) - 1
        #end if
    #end while

    return partitionVersions
# end def


'''
Gives the possible number of combinations of staple versions, for all possible partitions.

@return: list of number of combinations *so far* for each partition - each element is the number of 
versions for the corresponding partition, plus the number of versions for all the partitions 
preceding it in the list. The last number in the list is thus the total number of combinations 
of the staple versions which add up to the loop length.
'''
def partitionVersionCountsForStrand(partitionsForStrand, staples, stapleVersions):
    
    partitionVersionNums = []
    partitionVersionsSoFar = 0
    for partition in partitionsForStrand:
        partitionVersionsSoFar += numVersionsForPartition(partition, staples, stapleVersions)
        partitionVersionNums.append(partitionVersionsSoFar)
        
    return partitionVersionNums
# end def


'''
Finds the number of ways to assign staple versions to staples for a given partition.
'''
def numVersionsForPartition(partition, staples, stapleVersions):
    numVersions = 1
    
    def numCombinationsWithReplacement(n, k):
        n1 = n + k - 1
        return int(math.factorial(n1) / (math.factorial(k) * math.factorial(n1 - k)))
    # end def
    
    for stapleIndex in range(len(staples)):
        n = len(stapleVersions[stapleIndex])
        k = partition.count(staples[stapleIndex])
        numVersions *= numCombinationsWithReplacement(n, k)
    #end for
    
    return numVersions
#end def