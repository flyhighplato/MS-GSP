'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import logging
import re
import sys

#### Initialization utilities

# Set up logging
def initLogger():  
    logging.basicConfig( level=logging.DEBUG,
                         format='%(asctime)s %(levelname)s %(message)s',
                         filename='msgsp.log',
                         filemode='w' )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)s %(levelname)s: %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

# Loads data file into a database of sequences, where each sequence is a series list of transactions
def loadData( db, fileName ):
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        seqStr = re.split( r"<{|}{|}>", line )
        seqLst = []
        for transStr in seqStr:
            if ( transStr != '' ):
                transLst = []
                transItemsLst = transStr.split( ", " ) 
                for transItem in transItemsLst:
                    transLst.append( int(transItem) )
                seqLst.append( transLst )
        db.append( seqLst )

# Loads parameters for minimum item support and support difference constraint from data file
def loadParams( map, fileName):
    sdc = sys.float_info.max
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        param = re.findall(r"\d+\.*\d*",line)
        # Make sure we're dealing with the two types of parameters that we can handle
        assert( ( len( param ) == 2 and line.startswith("MIS") ) or ( len( param ) == 1 and line.startswith( "SDC" ) ) )
        # Specifying minimum support for an item
        if(line.startswith("MIS")):
            map[ int(param[0]) ] = float(param[1]) # mapping item id -> MIS
        # Specifying support difference constraint
        else:
            sdc = float(param[0])
    return sdc

# Sorts transactions in parameter sequence database by user supplied MIS values
def sortData( seqDB, misMap):
    for seq in seqDB:
        for trans in seq:
            trans.sort(key=lambda x:misMap[x])

#### Sequence utilities

# Returns True if seqA contains seqB, False otherwise
def seqContains( seqSup, seqSub ):
    if ( len( seqSub) == 0 ):
        return True
    idxSub = 0
    for idxSup in range( 0, len( seqSup ) ):
        if ( set( seqSub[ idxSub ] ).issubset( set(seqSup[ idxSup ] ) ) ):
            idxSub += 1
            if ( idxSub == len( seqSub ) ):
                return True
    return False

# Structure for representing a sequence object
class Sequence:
    # Constructor
    def __init__(self, seq=[], support=-1.0):
        self.seq = seq
        self.support = support
    
    # String representation
    def __repr__(self):
        return str(self.seq) + ":" + str(self.support)
 
    # Computes and caches support for this sequence
    def cacheSupport(self, seqDB):
        support = 0.0
        for seq in seqDB:
            if ( seqContains( seq, self.seq ) ):
                support += 1.0
        support /= float(len(seqDB))
        
    # Returns support from [0.0 to 1.0] of this sequence - will assert if support is invalid (not been cached)
    def getSupport(self):
        assert(( self.support >= 0.0 ) and ( self.support <= 1.0))
        return self.support
           

#### GSP Algorithm

# Initial pass for MS-GSP
def initPass( seqDB, misMap, L, F ):
    # Objectives:
    # 1. Find item 'I' with lowest MIS that has support greater than or equal to MIS(I)
    # 2. Output L : set of items that have support greater than MIS(I)
    # 3. Output F : set of 1-sequences such that support <{J}> is greater than or equal to MIS(J), F is a subset of L
    
    # Count unique items in sequence database
    itemCountsMap = {}
    for seq in seqDB:
        uniqueItemsInSeq = set()
        for trans in seq:
            for item in trans:
                uniqueItemsInSeq.add(item)
        for item in uniqueItemsInSeq:
            if item not in itemCountsMap:
                itemCountsMap[item]  = 1.0
            else:
                itemCountsMap[item] += 1.0
                    
    # Make possible 1-sequences
    for key in itemCountsMap.keys():
        support = itemCountsMap[key] / float( len( seqDB ) )
        assert( ( 0.0 <= support ) and ( 1.0 >= support ) )
        L.append( Sequence( [key], support ) )
    
    # Sort possible 1-sequences by MIS
    L.sort( key=lambda x:misMap[ x.seq[0] ] )
    
    # Determine first item with support >= MIS(item) - this item will have the lowest satisfied MIS
    toItemId = lambda oneSeq : oneSeq.seq[ 0 ]
    idxItemWithLowestSatisfiedMIS = len( L )
    for idx in range( 0, len( L ) ):
        if ( L[idx].getSupport() >= misMap[ toItemId( L[idx] ) ] ):
            idxItemWithLowestSatisfiedMIS = idx
            break
    
    # Trim L of all items with support lower than item with lowest satisfied MIS
    # Note use of [:] to modify list in place rather than allocate a new list
    minGlobalSatisfiedMis = misMap[ toItemId( L[idxItemWithLowestSatisfiedMIS] ) ]
    satisfiesMinGlobalSatisfiedMis = lambda seq : (seq.getSupport() >= minGlobalSatisfiedMis)
    L[:] = [ seq for seq in L if satisfiesMinGlobalSatisfiedMis( seq ) ]
    
    # Determine frequent 1-sequences
    for seq in L:
        if ( seq.getSupport() >= misMap[ toItemId( seq ) ] ):
            F.append( seq )

    logging.getLogger("InitPass").info("L = " + str(L))
    logging.getLogger("InitPass").info("F = " + str(F))

# Main body of MS-GSP
def MSGSPMain():       
    seqDB = []                           # The sequence database
    misMap = {}                          # A map of form item id -> minimum item support
    dataPath = "../../../Data/data.txt"  # The path to the input data @TODO: Read from arguments!
    paramPath = "../../../Data/para.txt" # The path to the MIS parameter data @TODO: Read from arguments!

    # Initialize logging
    initLogger()
    
    # Load sequence data from file
    loadData( seqDB, dataPath )
    logging.getLogger("MSGSPMain").info("Loaded seqDB: " + str(seqDB))
    
    # Initialize MIS values and support difference constraint threshold @TODO: support SDC!
    sdc = loadParams(misMap,paramPath)
    logging.getLogger("MSGSPMain").info("sdc: " + str(sdc))
    
    # Sort data according to MIS values
    sortData( seqDB, misMap )
    logging.getLogger("MSGSPMain").info("Sorted seqDB: " + str(seqDB))
    
    # Generate all frequent 1-sequences
    L=[]
    F=[]
    initPass( seqDB, misMap, L, F ) 
    
#### Application entry point

if __name__ == '__main__':
    MSGSPMain()
