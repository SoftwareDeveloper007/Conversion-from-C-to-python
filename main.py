import numpy as np
import sys


# ***********************************************************************************#
#   U(w,l,n,p) gives the total number of reads uncorrectable with witLength = w     #
#   w = witLength, l = readLength, n = numberOfReads, p = error                     #
# ***********************************************************************************#

def U(witLength, readLength, numberOfReads):
    sum = 0
    temp = 0
    error = 0.01

    fw = np.zeros((readLength + 1, readLength))

    for j in range(1, np.min(witLength, readLength) + 1):
        fw[1][j] = j

    for j in range(witLength + 1, np.min(2 * witLength - 1, readLength) + 1):
        fw[1][j] = 2 * witLength - j

    for j in range(2 * witLength, readLength + 1):
        fw[1][j] = 0

    for i in range(2, readLength + 1):
        fw[i][1] = 0

        for j in range(2, np.min(witLength + 1, readLength) + 1):
            sum += fw[i - 1][j - 1]
            fw[i][j] = int(sum)

        for j in range(witLength + 2, readLength + 1):
            sum -= fw[i - 1][j - witLength - 1]
            sum += fw[i - 1][j - 1]
            fw[i][j] = int(sum)

    sum = 0
    temp = 2
    k = 1

    while k <= readLength:

        temp = float(fw[k][readLength]) * error ** k * (1 - error) ** (readLength - k) * float(numberOfReads)
        sum += temp
        k += 1
        if temp > 0 and temp < 1:
            k = readLength + 1

    return sum


# ***********************************************************************************#
#   D(w,l,n,L,p) = expected number of destructible reads, that is, correct reads    #
#   that are turned wrong because they contain va with v containing errors,         #
#   a correct but v appearing elsewhere as correct vb                               #
# ***********************************************************************************#

def D(witLength, readLength, numberOfReads, genomeLength):
    # ** v has some erros: (1 - (1 - p) ^ w)
    # ** a is correct: (1 - p)
    # ** v appears elsewhere: (1 - (1 - 1 / 4 ^ w) ^ L)
    # ** followed by b <> a: 3 / 4

    error = 0.01
    # ** prob 1 wit gone wrong
    qw = (1 - (1 - error) ** witLength) * (1 - error) * (1 - (1 - 0.25 ** witLength) ** genomeLength) * 0.75

    expCorrReads = (1 - error) ** readLength * float(numberOfReads)

    # ** 1 - (1 - qw) ^ (l - w) = prob 1 read gone wrong

    return (1 - (1 - qw) ** (readLength - witLength)) * expCorrReads


# ***********************************************************************************#
#   compute the paremeters needed for correction                                    #
#   witLength_small = wm and witLength_large = wM                                   #
#   no larger than 28 so that it fits into 7 bytes with 2 bits per base             #
# ***********************************************************************************#

def computeWitLength(readLength, numberOfReads, genomeLength, witLength_small, witLength_large):
    # ** compute wM, wm = witness lengths
    # #** wM = min (w | D(w) < 0.0001 * expErrReads)
    # ** wm = arg min (U(w) + D(w))
    error = 0.01

    expErrReads = (1 - (1 - error) ** readLength) * float(numberOfReads)
    temp = 0
    temp1 = 0
    temp2 = 0
    witLength_large = 1

    while D(witLength_large, readLength, numberOfReads, genomeLength) > 0.0001 * expErrReads:
        witLength_large += 1

    if witLength_large > 25:
        witLength_large = 25

    witLength_small = witLength_large

    temp = U(witLength_small, readLength, numberOfReads) + D(witLength_small, readLength, numberOfReads, genomeLength)

    temp1 = U(witLength_small - 1, readLength, numberOfReads) + D(witLength_small - 1, readLength, numberOfReads,
                                                                  genomeLength)

    temp2 = U(witLength_small + 1, readLength, numberOfReads) + D(witLength_small + 1, readLength, numberOfReads,
                                                                  genomeLength)

    if temp1 < temp:
        while temp1 < temp:
            witLength_small -= 1
            temp = U(witLength_small, readLength, numberOfReads) + D(witLength_small, readLength, numberOfReads,
                                                                     genomeLength)
            temp1 = U(witLength_small - 1, readLength, numberOfReads) + D(witLength_small - 1, readLength,
                                                                          numberOfReads,
                                                                          genomeLength)

    else:
        while temp2 < temp:
            witLength_small += 1
            temp = U(witLength_small, readLength, numberOfReads) + D(witLength_small, readLength, numberOfReads,
                                                                     genomeLength)
            temp2 = U(witLength_small + 1, readLength, numberOfReads) + D(witLength_small + 1, readLength,
                                                                          numberOfReads,
                                                                          genomeLength)

    if witLength_small > 24:
        witLength_small = 24

    if witLength_small >= witLength_large:
        witLength_large = witLength_small + 1


# ***********************************************************************************#
#   compute Threshold = T                                                            #
#   counts T or larger are correct, smaller are wrong                                #
# ***********************************************************************************#

def computeThreshold(readLength, numberOfReads, genomeLength, witLength, threshold):
    # compute T(done for each iteration)
    # T = min(k | Wc(k) > We(k) + 1)
    error = 0.01
    qc = (readLength - witLength) * (1 - error) ** (witLength + 1) / float(genomeLength)
    qe = (readLength - witLength) * (1 - error) ** witLength * error / float(3 * genomeLength)
    # while (dbinom(T, numberOfReads, qc) <= dbinom(T, numberOfReads, qe)) {
    threshold = 1

    while (1 - error) ** threshold * (1 - qc) ** (numberOfReads - threshold) <= (error / 3) ** threshold * (1 - qe) ** (
        numberOfReads - threshold):
        threshold += 1

    threshold += 2

# ***********************************************************************************#
#   Read input from the file provided by the user.  The reads are stored in a        #
#   character array and then converted to binary and stored in the binReads array.   #
# ***********************************************************************************#

def readInput(binReads, *argv, fileFormatFlag, totalBases, numberOfReads, maxReadLength, goodRead, readLocation,
              currentMemory, peakMemory):
    '''Create and initialize variables. '''


    baseCounter = 0
    readsSize = 0
    chunkSize = 0
    binReadSize = 0
    offset = 0
    k = 0
    i = 0
    j = 0
    totalNs = 0
    badReads = 0
    binReadsLocation = 0
    previousPeak = 0
    currReadLength = 0
    arraySize = 5000000
    oldArraySize = 0
    chooseLetter = 0

    try:
        inReadfile = open(argv[1], 'r')
    except:
        print("ERROR: Unable to open text file. ")
        sys.exit(1)

    readsSize = len(inReadfile.read())
    chunkSize = int(0.1*readsSize)

    '''Store the reads file into a character array. '''
    charReads = np.zeros((1, chunkSize))


