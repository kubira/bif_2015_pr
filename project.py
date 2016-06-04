#!/usr/bin/env python3

############################################################
#                                                          #
# Filename: project.py                                     #
# Author:   Bc. Radim Kubis, xkubis03                      #
# Date:     14.04.2015                                     #
#                                                          #
#  BIF Project 2014/2015 - Analysis of human genome GRCh3  #
#                                                          #
############################################################


### MODULE(S) IMPORT #######################################

# System-specific parameters and symbols
import sys


### CLASSES FOR SCRIPT #####################################

# Class for chromosome range
class ChromosomeRange:
    # Constructor for chromosome range
    #     chromosome - chromosome ID (String)
    #     start      - start position of range (Integer)
    #     stop       - stop position of range (Integer)
    def __init__(self, chromosome, start, stop):
        # Set chromosome ID of range
        self.chromosome = chromosome
        # Set start positon of range
        self.start = start
        # Set stop position of range
        self.stop = stop

    # Function for length of range
    #     return - length of range
    def length(self):
        # Length is difference of stop and start + 1
        # Return length of range
        return ((self.stop - self.start) + 1)

# Class for list of chromosomes ranges
class ChromosomesRangesList:
    # Constructor for list of chromosomes ranges
    def __init__(self):
        # Initialization of empty dictionary
        self.rangeList = {}

    # Function to add ChromosomeRange object to list
    #     r - ChromosomeRange object to add (ChromosomeRange)
    def add(self, r):
        # If chromosome ID is not in dictionary keys
        if r.chromosome not in self.rangeList:
            # Create new item for chromosome ID in dictionary
            self.rangeList[r.chromosome] = []
        # Append ChromosomeRange object to list of chromosome ranges
        self.rangeList[r.chromosome].append(r)

    # Function to reduce overlapping chromosome ranges
    def reduceList(self):
        # For every chromosome ID in dictionary
        for chrId in self.rangeList.keys():
            # Set index of range with maximum stop position to 0 (first range)
            maxStopId = 0
            # Create new list for reduced ranges of chromosome
            newList = []
            # For every range of chromosome (sorted ascending by start position)
            for r in sorted(self.rangeList[chrId], key=lambda x: x.start):
                # If actual range is overlapping range on maxStopId index
                if r.start <= self.rangeList[chrId][maxStopId].stop:
                    # Merge ranges (select max of choices):
                    #     actual range is inside - stop position does not change
                    #     actual range stops after max - set new stop position
                    self.rangeList[chrId][maxStopId].stop = max(
                        r.stop, self.rangeList[chrId][maxStopId].stop
                    )
                # If actual range is not overlapping previous range
                else:
                    # Add previous range to new list of reduced ranges
                    newList.append(self.rangeList[chrId][maxStopId])
                    # Set index of range with maximum stop position to actual
                    maxStopId = self.rangeList[chrId].index(r)
            # Append last (reduced) range to new list of ranges
            newList.append(self.rangeList[chrId][maxStopId])
            # Set new reduced list of ranges to chromosome
            self.rangeList[chrId] = newList

    # Function to get total length of ranges by chromosome ID
    #     chrId - chromosome ID (String)
    #
    #     return - total length of ranges by chromosome ID (String)
    def lengthByChromosomeId(self, chrId):
        # Total length of ranges by chromosome ID is sum of ranges lengths
        # Return sum of ranges lengths by chromosome ID
        return sum([r.length() for r in self.rangeList[chrId]])

    # Function to get total length of ranges for all chromosomes IDs
    #     return - total length of ranges for all chromosomes IDs
    def totalSize(self):
       # Total length of ranges in list is sum of lengths by chromosomes IDs
       # Return sum of all ranges in list
       return sum([self.lengthByChromosomeId(k) for k in self.rangeList.keys()])


### CONSTANTS FOR SCRIPT ###################################

# Exit without error
SUCCESS = 0
# Exit with some error
FAILURE = 1
# Chromosomes to analyze
CHROMOSOMES = [
     '1',  '2',  '3',  '4',  '5',  '6',
     '7',  '8',  '9', '10', '11', '12',
    '13', '14', '15', '16', '17', '18',
    '19', '20', '21', '22',  'X',  'Y',
]
# Lengths of chromosomes (1-22, X, Y)
CHROMOSOMES_LENGTHS = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
    145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
    101991189,  90338345,  83257441,  80373285,  58617616,  64444167,  46709983,
     50818468, 156040895,  57227415,
]
# Index of chromosome id
CHROMOSOME_ID = 0
# Index of feature (gene, start_codon, ...)
FEATURE = 2
# Index of start position
START = 3
# Index of stop position
STOP = 4
# Index of chromosome attributes
CHROMOSOME_ATTRS = 8
# Index of count
COUNT = 0
# Index of ranges list
LIST = 1
# Index of key
KEY = 0
# Index of value
VALUE = 1
# Number of expected arguments
EXPECTED = 2
# Groups of genes
CODING = [
    'protein_coding',
    'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene',
    'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene',
]
SMALL = [
    'snRNA', 'rRNA', 'snoRNA', 'miRNA', 'misc_RNA',
]
LONG = [
    'lincRNA', 'non_coding', 'processed_transcript', 'antisense',
    '3prime_overlapping_ncrna', 'sense_intronic', 'sense_overlapping',
    'known_ncrna',
]
PSEUDO = [
    'unitary_pseudogene', 'IG_C_pseudogene', 'translated_processed_pseudogene',
    'polymorphic_pseudogene', 'TR_J_pseudogene', 'IG_J_pseudogene',
    'TR_V_pseudogene', 'IG_V_pseudogene', 'pseudogene',
    'unprocessed_pseudogene', 'transcribed_unprocessed_pseudogene',
    'translated_unprocessed_pseudogene', 'transcribed_processed_pseudogene',
    'processed_pseudogene', 'transcribed_unitary_pseudogene',
]
PROTEIN = [
    'Coding transcripts', 'CDS',
]


### FUNCTIONS FOR SCRIPT ###################################

# Function for usage print
#     scriptName - name of running script (String)
def printUsage(scriptName):
    # Print usage to stderr
    print("Usage:\t%s <input_file>" % scriptName, file=sys.stderr)

# Function for error print
#     text - text of error message (String)
def printE(text):
    # Print error to stderr
    print("ERROR: %s" % text, file=sys.stderr)

# Function for warning print
#     text - text of warning message (String)
def printW(text):
    # Print warning to stderr
    print("WARNING: %s" % text, file=sys.stderr)


### VARIABLES FOR SCRIPT ###################################

# Number of command line arguments
numberOfInputArgs = len(sys.argv)
# Object for input file
inputFile = None
# Output table with statistics
outputTable = {}
# Total length of chromosomes 1-22, X, Y
genomeLength = float(sum(CHROMOSOMES_LENGTHS))


### PARSING ARGUMENTS ######################################

# If script is run without arguments
if numberOfInputArgs == 1:
    # Print error message
    printE("No input file defined.")
    # Print usage of script
    printUsage(sys.argv[0])
    # Exit with error
    exit(FAILURE)
# If script has at least one argument (input file name)
else:
    # If script has more than expected arguments
    if numberOfInputArgs > EXPECTED:
        # Print warning about ignoring extra argument(s)
        printW("Ignoring extra input argument(s): %s" % sys.argv[EXPECTED:])
    # Try to open input file
    try:
        # Opening input file
        inputFile = open(sys.argv[1], 'r')
    # Catch I/O error
    except IOError as e:
        # Print error message
        printE("%s: '%s'" % (e.strerror, sys.argv[1]))
        # Exit with error
        exit(FAILURE)


### PROCESSING INPUT FILE ##################################

# Read lines of input file
for line in inputFile:
    # Strip and split line
    splt = line.strip(';\n').split('\t')
    # Check if line is valid chromosome and gene/CDS feature to analyze
    if splt[CHROMOSOME_ID] in CHROMOSOMES and splt[FEATURE] in ['gene', 'CDS', 'transcript']:
        ### PROCESS LINE ###################################

        # Split attributes from line
        attributes = splt[CHROMOSOME_ATTRS].split(';')
        # Dictionary for attributes
        attrDict = {}
        # Create dictionary of attributes
        for attr in attributes:
            # Split attribute to key and value
            attrSplt = attr.strip().split(' ')
            # Create item in dictionary
            attrDict[attrSplt[KEY]] = attrSplt[VALUE].strip("\"")
        # If line has gene_biotype attribute
        if 'gene_biotype' in attrDict:
            # If feature is gene
            if splt[FEATURE] == 'gene':
                # If gene_biotype does not exist in outputTable statistics
                if attrDict['gene_biotype'] not in outputTable:
                    # Create new item with initial values:
                    #     count = 0
                    #     ranges = new ChromosomesRangesList()
                    outputTable[attrDict['gene_biotype']] = [
                        0,
                        ChromosomesRangesList(),
                    ]
                # Increment value of count
                outputTable[attrDict['gene_biotype']][COUNT] += 1
                # Add range to list
                outputTable[attrDict['gene_biotype']][LIST].add(
                    ChromosomeRange(
                        splt[CHROMOSOME_ID],
                        int(splt[START]),
                        int(splt[STOP])
                    )
                )
            # Else if gene_biotype is protein_coding
            elif attrDict['gene_biotype'] == 'protein_coding':
                # If feature is CDS
                if splt[FEATURE] == 'CDS':
                    # If CDS does not exist in outputTable statistics
                    if 'CDS' not in outputTable:
                        # Create new item for CDS with initial values:
                        #     count = 0
                        #     ranges = new ChromosomesRangesList()
                        outputTable['CDS'] = [
                            0,
                            ChromosomesRangesList(),
                        ]
                    # Increment value of CDS count
                    outputTable['CDS'][COUNT] += 1
                    # Add range to list
                    outputTable['CDS'][LIST].add(
                        ChromosomeRange(
                            splt[CHROMOSOME_ID],
                            int(splt[START]),
                            int(splt[STOP])
                        )
                    )
                # Else if feature is transcript
                elif splt[FEATURE] == 'transcript':
                    # If Coding transcripts does not exist in outputTable statistics
                    if 'Coding transcripts' not in outputTable:
                        # Create new item for Coding transcripts with initial values:
                        #     count = 0
                        #     ranges = new ChromosomesRangesList()
                        outputTable['Coding transcripts'] = [
                            0,
                            ChromosomesRangesList(),
                        ]
                    # Increment value of Coding transcripts count
                    outputTable['Coding transcripts'][COUNT] += 1
                    # Add range to list
                    outputTable['Coding transcripts'][LIST].add(
                        ChromosomeRange(
                            splt[CHROMOSOME_ID],
                            int(splt[START]),
                            int(splt[STOP])
                        )
                    )
        # Else if line does not have gene_biotype attribute
        else:
            # Print warning about ingnoring line without gene_biotype attribute
            printW("Ignoring line without gene_biotype attribute: '%s'" % line)

# Close input file
inputFile.close()


### PRINTING OUTPUT TABLE ##################################

### CODING GENES ###########################################

# Initialization of variables for summary values
count = size = gcov = 0
# Header print
print("Coding genes;Count;Size [bp];G.cov. [%]")
# Values print
for name in CODING:
    # Reduce ranges of gene_biotype by name from CODING
    outputTable[name][LIST].reduceList()
    # Save size of actual name
    nameSize = outputTable[name][LIST].totalSize()
    # Save percent for actual name
    namePercent = (nameSize / genomeLength) * 100
    # Print values
    print("%(name)s;%(count)d;%(size)d;%(percent).6f" % {
            'name': name,
            'count': outputTable[name][COUNT],
            'size': nameSize,
            'percent': namePercent,
        }
    )
    # Summary count
    count += outputTable[name][COUNT]
    # Summary size
    size += nameSize
    # Summary percent
    gcov += namePercent
# Summary print
print("Summary;%(count)d;%(size)d;%(gcov).6f" % {
        'count': count, 'size': size, 'gcov': gcov,
    }
)
# Empty line print
print(";;;")


### SMALL NON-CODING GENES #################################

# Initialization of variables for summary values
count = size = gcov = 0
# Header print
print("Small non-coding genes;Count;Size [bp];G.cov. [%]")
# Values print
for name in SMALL:
    # Reduce ranges of gene_biotype by name from SMALL
    outputTable[name][LIST].reduceList()
    # Save size of actual name
    nameSize = outputTable[name][LIST].totalSize()
    # Save percent for actual name
    namePercent = (nameSize / genomeLength) * 100
    # Print values
    print("%(name)s;%(count)d;%(size)d;%(percent).6f" % {
            'name': name,
            'count': outputTable[name][COUNT],
            'size': nameSize,
            'percent': namePercent,
        }
    )
    # Summary count
    count += outputTable[name][COUNT]
    # Summary size
    size += nameSize
    # Summary percent
    gcov += namePercent
# Summary print
print("Summary;%(count)d;%(size)d;%(gcov).6f" % {
        'count': count, 'size': size, 'gcov': gcov,
    }
)
# Empty line print
print(";;;")


### LONG NON-CODING GENES ##################################

# Initialization of variables for summary values
count = size = gcov = 0
# Header print
print("Long non-coding genes;Count;Size [bp];G.cov. [%]")
# Values print
for name in LONG:
    # Reduce ranges of gene_biotype by name from LONG
    outputTable[name][LIST].reduceList()
    # Save size of actual name
    nameSize = outputTable[name][LIST].totalSize()
    # Save percent for actual name
    namePercent = (nameSize / genomeLength) * 100
    # Print values
    print("%(name)s;%(count)d;%(size)d;%(percent).6f" % {
            'name': name,
            'count': outputTable[name][COUNT],
            'size': nameSize,
            'percent': namePercent,
        }
    )
    # Summary count
    count += outputTable[name][COUNT]
    # Summary size
    size += nameSize
    # Summary percent
    gcov += namePercent
# Summary print
print("Summary;%(count)d;%(size)d;%(gcov).6f" % {
        'count': count, 'size': size, 'gcov': gcov,
    }
)
# Empty line print
print(";;;")


### PSEUDO GENES ###########################################

# Initialization of variables for summary values
count = size = gcov = 0
# Header print
print("Pseudo genes;Count;Size [bp];G.cov. [%]")
# Values print
for name in PSEUDO:
    # Reduce ranges of gene_biotype by name from PSEUDO
    outputTable[name][LIST].reduceList()
    # Save size of actual name
    nameSize = outputTable[name][LIST].totalSize()
    # Save percent for actual name
    namePercent = (nameSize / genomeLength) * 100
    # Print values
    print("%(name)s;%(count)d;%(size)d;%(percent).6f" % {
            'name': name,
            'count': outputTable[name][COUNT],
            'size': nameSize,
            'percent': namePercent,
        }
    )
    # Summary count
    count += outputTable[name][COUNT]
    # Summary size
    size += nameSize
    # Summary percent
    gcov += namePercent
# Summary print
print("Summary;%(count)d;%(size)d;%(gcov).6f" % {
        'count': count, 'size': size, 'gcov': gcov,
    }
)
# Empty line print
print(";;;")


### PROTEIN CODING GENES ###################################

# Header print
print("Protein Coding Genes;Count;Size [bp];G.cov. [%]")
# Values print
for name in PROTEIN:
    # Reduce ranges of gene_biotype by name from PSEUDO
    outputTable[name][LIST].reduceList()
    # Save size of actual name
    nameSize = outputTable[name][LIST].totalSize()
    # Save percent for actual name
    namePercent = (nameSize / genomeLength) * 100
    # Print values
    print("%(name)s;%(count)d;%(size)d;%(percent).6f" % {
            'name': name,
            'count': outputTable[name][COUNT],
            'size': nameSize,
            'percent': namePercent,
        }
    )


### END OF SCRIPT ##########################################

# Exit without error
exit(SUCCESS)

