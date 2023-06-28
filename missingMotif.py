'''
File: missingMotif.py
executable: missingMotif.py
Purpose: find missing genomic sequences in .fna files
Input/Output: STDIN/STDOUT

Options:

--minMotif  minimum motif size to evaluate (int)
--maxMotif maximum motif size to evaluate (int)
--cutoff  Z-score cutoff (int)

Ex: python missingMotif.py --minMotif 3 --maxMotif 8 --cutoff 0 < Zm4-genomic.fna > output.txt

Author: Quin Lamothe
Due: 10/12/2020
Created 10/6/2020

'''

import sys
import argparse

class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
                header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class CommandLine:
    """
    Allows use of command line for program inputs.
    Input: command line rendered variables.
    Output: values assigned to variables.
    """

    def __init__(self, inOpts = None):

        self.parser = argparse.ArgumentParser()

        self.parser.add_argument('-m', '--minMotif', type=int,
                                 action='store', default=3, help='min motif length')
        self.parser.add_argument('-x', '--maxMotif', type=int,
                                 action='store', default=8, help='max motif length')
        self.parser.add_argument('-c', '--cutoff', type=float,
                                 action='store', default=0.0, help='cutoff for max z value')
        self.parser.add_argument('-k', '--kScoring', type=bool,
                                 action='store', default=False, help='k values')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    """
    Used to signal a Usage error, evoking a usage
    statement and eventual exit when raised.
    """

    def __init__(self, msg):
        self.msg = msg

class workDoer:
    """
    This class is designed to preform the work necessary to develop information pertaining
    to the expectations, z and k values of specific kmers in a given dataset.

    Functions:

    revComp(): Returns the reverse complement of a given kmer.

    countDict(): Develops a dictionary holding counts of both kmer and reverse kmer.

    expected(): Develops expectations and z scores for each given kmer.

    kScores(): This function is designed to develop k scores for each given kmer.

    printer(): This function organizes and prints the calculations as relate to each given kmer.

    """

    def __init__(self, minMotif, maxMotif, cutoff, kScoring):

        self.minMotif = minMotif
        self.maxMotif = maxMotif
        self.cutoff = cutoff
        self.kScoring = kScoring

#        if kScoring == True:
#            self.kScoring = True
#        else:
#            self.kScoring = False

        self.rSeqDict = {}
        self.kmerList = []
        self.zScoreList = []
        self.eDict = {}
        self.pDict = {}
        self.cDict = {}
        self.cDictList = []
        self.zDict = {}
        self.printDictList = []
        self.finalTupList = []
        self.kDict = {}
        self.kSizeDict = {}
        self.normDict = {}
        self.printDict = {}
        self.n = 0

    def revComp(self, kmer):
        """
        Function designed to return the reverse complement of each underrepresented sequence.

        input: sequence
        output: reverse complement
        """

        rKmer = kmer.lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G').replace(' ', '')[::-1]

        self.rSeqDict[kmer] = rKmer

        return rKmer

    def countDict(self, seq):
        """
        Function which joins the sequence and reverse complements,
        places them into a dictionary as keys, and then attaches the number
        of occurances.

        --- reverse complements -- kmer may have reverse complement that isnt the kmer -- report as single count
        kmer------>ATCG, irrelevant ------> TAGC,  valid ----> reverse CGAT
        take count of kmer and reverse complement summed when calculating null models
        --- palindromes, disregard
        ---reverse complements palindromes  ----- dont add to both, just count once

        input: sequence fragments
        Output: a list of dictionaries containing frequencies of sequence fragments of given lengths
        with associated count developed by recognition
        """

        dec = self.maxMotif

        #using n as the total length of the genome
        acceptableBases = ['A', 'G', 'C', 'T', 'N']

        for base in seq:
            if base in acceptableBases:
                self.n += len(base)

        #while loop allows decrementation to capture kmers 1 - max
        while dec > 0:

            # 1st forloop iterates through all indexes of sequence
            for seqPosition in range(0, len(seq) - dec - 1):
                # kmers are taken by string slicing
                kmer = seq[seqPosition: seqPosition + dec]

                # should use only the count of the acceptable bases for the total length of the genomic sequence
                #for base in kmer:
                #    if base not in acceptableBases:
                #        continue

                # adding to the counts related to each kmer's presence
                if kmer in self.cDict.keys():
                    self.cDict[kmer] += 1
                else:
                    self.cDict[kmer] = 1

                #reverse complement palendromes should not be counted with the palendrome
                if self.revComp(kmer) in self.cDict.keys():
                    self.cDict[self.revComp(kmer)] += 1
                    self.cDict[kmer] += 1

                    if self.revComp(kmer)[::-1] in self.cDict.keys():
                        self.cDict[self.revComp(kmer)] -= 1
                        self.cDict[kmer] -= 1
                else:
                    self.cDict[self.revComp(kmer)] = 1

            dec -= 1

    def expected(self):
        """
        Develops a float value for the expectation placed on encountering a particular kmer.

        kmers 'to the left' * kmers 'to the right' will be divided by the 'middle kmers'

        Input:dictionary of counts
        Output: values for expectations of appearance of given kmer
        """

        for kmer in self.cDict.keys():
            #list of expectations E(K) developed per kmer --
            # EXAMPLE: (c(k1,k2) * c(k2,k3)) / c(k2) = E(K3)

            if len(kmer) > self.minMotif - 1:

                #self.printDict[kmer] = []

                # 1 find count of kmer with first position missing
                suffixPosMissing = self.cDict[kmer[1:]] + self.cDict[self.revComp(kmer)[1:]]

                # 2 find count of kmer with last position missing
                prefixPosMissing = self.cDict[kmer[:-1]] + self.cDict[self.revComp(kmer)[:-1]]

                # 3 find count of kmer with first and last position missing
                prefixSuffix = self.cDict[kmer[1:- 1]] + self.cDict[self.revComp(kmer)[1:-1]]

                E = ((suffixPosMissing * prefixPosMissing) / prefixSuffix)

                E = E * 2

                p = E / self.n

                sd = pow(E * (1 - p), .5)

                z = (self.cDict[kmer] - E) / sd

                self.printDict[kmer] = [kmer, self.revComp(kmer), self.cDict[kmer], E, z]

    def kScores(self):
        """
        This function is designed to develop k scores used to rank z scores.

        Input: expectations, counts of kmers
        Output: k scores for each kmer
        """

        dec = self.maxMotif
        normList = []
        meankDict = {}
        sdkDict = {}
        sumNormDict = {}
        sumNormSquaredDict = {}

        #normalize k scores
        for i in range(1, self.maxMotif + 1):
            self.kSizeDict[i] = self.n / i

        while dec != self.minMotif:

            for kmer in self.eDict.keys():
                # developing the norm for each kmer by dividing the E(K) for each kmer by the
                # total number of kmers of the same length

                #developing a normal value for kmers of a given size to be summed
                normList.append(self.kSizeDict[dec] / self.eDict[kmer])

                # forming a sum of the normalized values -- there should be one of these for kmers of each given size
            sumNormDict[dec] = sum(normList)

            # square of that sum
            sumNormSquaredDict[dec] = pow(sumNormDict[dec], 2)

            #developing a mean k value for kmers of a given size
            meankDict[dec] = sumNormDict[dec] / self.kSizeDict[dec]

            #foming a standard deviation for kmers of a specific size
            sdkDict[dec] = pow((sumNormSquaredDict[dec] / self.kSizeDict[dec]) - ((pow(sumNormDict[dec], 2)) / pow(self.kSizeDict[dec], 2)), .5)

            # calculate z scores over the normalized score found for
            # each kmer within each group size
            for kmer in self.cDict.keys():

                if skdDict[dec] != 0:
                    k = (self.cDict[kmer] - meankDict[dec]) / sdkDict[dec]
                    self.kDict[kmer] = k
                    self.printDict[kmer].append(k)

            dec -= 1

    def printer(self):
        """
        Function designed to handle the formatting and printing of the data developed.

        Input: kmers, reverse complement of kmers, counts of appearances, expected values (z scores),
        probabilities of occurances.

        Output: formatted ^^
        """
        print("N: ", self.n - 8)

        list3 = []
        list4 = []
        list5 = []
        list6 = []
        list7 = []
        list8 = []

        dec = self.maxMotif

        # the lists of kmer and reverse complement are sorted alphabetically -- then rearranged
        # checking if there are duplicates in the alphabetically ordered kmers,
        # as reverse complements will be ordered the same, causing duplicate items
        # these items will be happened to a list

        # flip kmer pairs to be alphabetically ordered

        #sortList = []
        for key in self.printDict:
            #make list of lists -- item[0] is kmer and item[1] rev complement
            sortedAlph = sorted([key, self.printDict[key][1]])
            self.printDict[sortedAlph[0]][0] = sortedAlph[0]
            self.printDict[sortedAlph[0]][1] = sortedAlph[1]
            #sortList.append(sortedAlph[0])

        #print(item)
        """
        #print(sortList)
        dupCounterDict = {}
        dupList = []
        for item in sortList:
            duplicateCheck = sortList.count(item)
            if duplicateCheck > 1:
                dupCounterDict[item] = duplicateCheck
                dupList.append(item)

        #for dup in dupList:
        #    if dup in self.printDict.keys():
        #        self.printDict[dup][3] *= dupCounterDict[dup]

        dupRemoverDict = {}
        for key in self.printDict.keys():
            if key not in dupList:
                dupRemoverDict[key] = self.printDict[key]

        self.printDict = dupRemoverDict
        """

        for key in self.printDict:
            if len(key) == 8:
                list8.append(self.printDict[key])

            if len(key) == 7:
                list7.append(self.printDict[key])

            if len(key) == 6:
                list6.append(self.printDict[key])

            if len(key) == 5:
                list5.append(self.printDict[key])

            if len(key) == 4:
                list4.append(self.printDict[key])

            if len(key) == 3:
                list3.append(self.printDict[key])

        list8.sort(key=lambda x: x[4], reverse=False)
        list7.sort(key=lambda x: x[4], reverse=False)
        list6.sort(key=lambda x: x[4], reverse=False)
        list5.sort(key=lambda x: x[4], reverse=False)
        list4.sort(key=lambda x: x[4], reverse=False)
        list3.sort(key=lambda x: x[4], reverse=False)
    
        finalPrintDict = {8: list8, 7: list7, 6: list6, 5: list5, 4: list4, 3: list3}

        if self.kScoring == True:
            for key in finalPrintDict.keys():
                for item in finalPrintDict[key]:
                    if item[4] < self.cutoff:
                        print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}\t{5:0.3f}'.format(item[0], item[1], item[2], item[3],
                                                                               item[4], item[5]))
        else:
            for key in finalPrintDict.keys():
                for item in finalPrintDict[key]:
                    if item[4] < self.cutoff:
                        print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(item[0], item[1], item[2], item[3], item[4]))


def main(myCommandLine=None):

    if myCommandLine == None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(myCommandLine)
    try:
        print(myCommandLine.args)

    except Usage as err:
        print(err.msg)

    myReader = FastAreader()

    #linker allows handling of sys.stdin
    linker = workDoer(myCommandLine.args.minMotif, myCommandLine.args.maxMotif, myCommandLine.args.cutoff,
                      myCommandLine.args.kScoring)

    #feeds sequences from fasta file into function
    for head, seq in myReader.readFasta():
        linker.countDict(seq)
    linker.expected()

    if linker.kScoring == True:
        print("sequence ", "reverse ", " count ", " expect ", " Z score", " k score")
        linker.kScores()
        linker.printer()
    else:
        print("sequence ", "reverse ", " count ", " expect ", " Z score")
        linker.printer()

if __name__ == '__main__':
    main()
