'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import logging
import math
import re
import sys

####

# Structure for easier passing around of "global" parameters
class Context:
    # Constructor
    def __init__(self):
        self.seqDB = []                # The sequence database
        self.misMap = {}               # A map of form item id -> minimum item support
        self.supportMap = {}           # A map of form item id -> actual support
        self.sdc = sys.float_info.max  # The maximum support difference constraint allowed between two sequences      
        
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
def loadParams( map, fileName ):
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
def sortData( seqDB, misMap ):
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
    def __init__( self, seq=[], mis=-1.0, count=-1.0, support=-1.0 ):
        self.seq = seq
        self.mis = mis
        self.count = count
        self.support = support
    
    # String representation
    def __repr__(self):
        strSeq = "<"
        for trans in self.seq:
            strSeq += "{"
            if ( len( trans ) >= 1 ):
                strSeq += str( trans[0] )
            for idx in range ( 1, len( trans ) ):
                strSeq += "," + str( trans[ idx ] )
            strSeq += "}"
        strSeq += "> Count: " + str( int(self.getCount()) )
        return strSeq
    
    # Returns MIS value assumed to be greater than or equal to 0.0
    def getMis(self):
        assert( self.mis >= 0.0 )
        return self.mis
   
    # Computes and caches support for this sequence
    # Also, of note: the optimal way to avoid reads from the database is to compute the sequence counts for a single record in the database
    # (i.e . for db for seqs instead of for seqs for db as is the case here)
    def cacheSupport(self, seqDB):
        self.count = 0.0
        for seq in seqDB:
            if ( seqContains( seq, self.seq ) ):
                self.count += 1.0
        self.support = self.count / float(len(seqDB))
        
    # Returns support from [0.0 to 1.0] of this sequence - will assert if support is invalid (not been cached)
    def getSupport(self):
        assert( ( self.support >= 0.0 ) and ( self.support <= 1.0 ) )
        return self.support
    
    # Returns count for this sequence - will assert if count is invalid (not been cached)
    def getCount(self):
        assert( self.count >= 0.0 )
        return self.count
    
    # Returns item id of first item in sequence
    def getFirstItemId(self):
        assert( (len( self.seq ) > 0) and (len( self.seq[0] ) > 0 ) )
        return self.seq[0][0]
        
    # Returns item id of last item in sequence
    def getLastItemId(self):
        assert( (len( self.seq ) > 0) and (len( self.seq[0] ) > 0 ) )
        return self.seq[-1][-1]

    # Returns True if item with parameter id contains lowest MIS and MIS is unique within the sequence, False otherwise
    def itemHasUniqueMinMis(self, itemId, misMap):
        assert( itemId in misMap )
        itemMis = misMap[ itemId ]
        itemCount = 0
        for trans in self.seq:
            for transItemId in trans:
                if ( transItemId != itemId ):
                    assert( itemId in misMap )
                    if ( itemMis >= misMap[ transItemId ] ):
                        return False
                else:
                    itemCount += 1
                    if ( itemCount > 1 ):
                        return False
        assert( itemCount > 0 )
        return True;
    
    # Returns True if first item contains lowest MIS and MIS is unique within sequence, False otherwise
    def firstItemHasUniqueMinMis(self, misMap):
        return self.itemHasUniqueMinMis( self.getFirstItemId(), misMap )
    
    # Returns True if last item contains lowest MIS and MIS is unique within sequence, False otherwise
    def lastItemHasUniqueMinMis(self, misMap):
        return self.itemHasUniqueMinMis( self.getLastItemId(), misMap )
    
    # Returns True if removing first element from this sequence and last element from parameter sequence results in same sequence, False otherwise
    def canJoin(self, seqObj):
        # Bounds checking
        assert( (len( self.seq ) > 0) and (len( self.seq[0] ) > 0 ) )
        assert( (len( seqObj.seq ) > 0) and (len( seqObj.seq[0] ) > 0 ) )
        # Remove first element from this sequence (copying sequences because k should be small)
        tmpDelFirst = self.seq[:]
        del tmpDelFirst[0][0]
        assert( tmpDelFirst != self.seq )
        # Remove last element from parameter sequence
        tmpDelLast = seqObj.seq[:]
        del tmpDelLast[-1][-1]
        assert( tmpDelLast != seqObj.seq )
        # Return True if both sequences are the same, False otherwise
        return (tmpDelFirst == tmpDelLast)
    
    # Returns new raw sequence with the result of joining this sequence with parameter sequence
    def join(self, seqObj):
        # Bounds checking
        assert( (len( self.seq ) > 0) and (len( self.seq[0] ) > 0 ) )
        assert( (len( seqObj.seq ) > 0) and (len( seqObj.seq[0] ) > 0 ) )
        # Get item id to merge into this sequence
        mergeItemId = seqObj.getLastItemId()
        joinedRawSeq = self.seq[:]
        # Determine if last item should be appended to last transaction or if a new transaction should be created
        if ( len( seqObj.seq[-1] ) > 1 ):
            joinedRawSeq[-1].append( mergeItemId ) 
        else:
            joinedRawSeq.append( [mergeItemId] )
        return joinedRawSeq

# Returns True if parameter raw sequence is not found within list, False otherwise
def isUniqueSeqWithinList( aList, aRawSeq ):
    for seqObj in aList:
        if seqObj.seq == aRawSeq:
            return False
    return True

# Appends the raw sequence as a sequence object and caches the support of the sequence
# Using this utility function because we can't overload constructors in Python
def appendSeqAndCacheSupport( aList, aRawSeq, anMis, seqDB ):
    assert( isUniqueSeqWithinList( aList, aRawSeq ) )
    aList.append( Sequence( aRawSeq, anMis ) )
    aList[ -1 ].cacheSupport( seqDB )

# Extracts all sequences from C which have a support greater than or equal to their MIS and stores them in F   
def extractAllSeqsWhichSatisfyTheirMis( F, C ):
    F[:] = [ seq for seq in C if ( seq.getSupport() >= seq.getMis() ) ]

#### GSP Algorithm

# Initial pass for MS-GSP
def initPass( L, F, ctx ):
    # Objectives:
    # 1. Find item 'I' with lowest MIS that has support greater than or equal to MIS(I)
    # 2. Output L : set of items that have support greater than MIS(I)
    # 3. Output F : set of 1-sequences such that support <{J}> is greater than or equal to MIS(J), note: F is a subset of L
    
    # Count unique items in sequence database
    itemCountsMap = {}
    for seq in ctx.seqDB:
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
    L[:] = []
    for key in itemCountsMap.keys():
        count = itemCountsMap[ key ]
        support = count / float( len( ctx.seqDB ) )
        assert( ( 0.0 <= support ) and ( 1.0 >= support ) )
        L.append( Sequence( [[key]], ctx.misMap[ key ], count, support ) )
        ctx.supportMap[key] = support # Not sure if we need this!
    
    # Sort possible 1-sequences by MIS
    L.sort( key=lambda seqObj : seqObj.getMis() )
    
    # Determine first item with support >= MIS(item) - this item will have the lowest satisfied MIS
    idxItemWithLowestSatisfiedMIS = len( L )
    for idx in range( 0, len( L ) ):
        if ( L[idx].getSupport() >= L[idx].getMis() ):
            idxItemWithLowestSatisfiedMIS = idx
            break
    
    # Trim L of all items with support lower than item with lowest satisfied MIS
    # Note use of [:] to modify list in place rather than allocate a new list
    minGlobalSatisfiedMis = L[idxItemWithLowestSatisfiedMIS].getMis()
    L[:] = [ seq for seq in L if (seq.getSupport() >= minGlobalSatisfiedMis) ]
    
    # Determine frequent 1-sequences
    extractAllSeqsWhichSatisfyTheirMis( F, L )
    
    logging.getLogger("InitPass").info("L = " + str(L))
    logging.getLogger("InitPass").info("F = " + str(F))

# Determines candidate 2-sequences
def level2CandidateGen( C, L, ctx ):
    C[:] = []
    for idxL in range(0, len(L)):
        seqL = L[ idxL ]
        # Assert we're working with sequences of length 1
        assert( ( len( seqL.seq ) == 1 ) and ( len( seqL.seq[0] ) == 1 ) )
        # See if 'l' satisfies it's own MIS
        if ( seqL.getSupport() >= seqL.getMis() ):
            lId = seqL.getFirstItemId()
            # BEGIN UNSUPPORTED CODE BLOCK
            # Also create 2-tuple <{l}{l}> - it could exist!
            appendSeqAndCacheSupport( C, [ [ lId ], [ lId ] ], seqL.getMis(), ctx.seqDB )
            # END UNSUPPORTED CODE BLOCK
            # Create 2-tuples with all sequences 'h' where MIS(h) >= MIS(l)  
            for idxH in range( idxL+1, len(L) ):
                seqH = L[ idxH ]
                # Assert sequences were pre-sorted by MIS
                assert( seqH.getMis() >= seqL.getMis() )
                # Only create tuple if sup(h) is greater than MIS(l)
                if ( ( seqH.getSupport() >= seqL.getMis() ) and ( math.fabs( seqH.getSupport() - seqL.getSupport() ) <= ctx.sdc ) ):
                    # Join seqL and seqH to create both <{l,h}> and <{l,h}>
                    hId = seqH.getFirstItemId()
                    assert ( hId != lId ) # assert these items are unique!
                    appendSeqAndCacheSupport( C, [ [ lId, hId ] ], seqL.getMis(), ctx.seqDB )
                    appendSeqAndCacheSupport( C, [ [ lId ], [ hId ] ], seqL.getMis(), ctx.seqDB )
                    # BEGIN UNSUPPORTED CODE BLOCK
                    # Also create 2-tuple <{h},{l}> - else how could we get it?
                    appendSeqAndCacheSupport( C, [ [ hId ], [ lId ] ], seqL.getMis(), ctx.seqDB )
                    # END UNSUPPORTED CODE BLOCK

# Determines candidate k-sequences where k is not 2
def MSCandidateGenSPM( C, Fprev, ctx ):
    # Join step: create candidate sequences by joining Fk-1 with Fk-1
    # NOTE: seq1 joins seq2 and seq2 joins with seq1 iff seq1 = <abab...ab> and seq2 = <baba..ba>
    for idx1, seq1 in enumerate( Fprev ):
        if ( seq1.firstItemHasUniqueMinMis( ctx.misMap ) ):
            # @TODO:
            return
        elif ( seq1.lastItemHasUniqueMinMis( ctx.misMap ) ):
            # @TODO:
            return
        else:
            # BEGIN UNSUPPORTED CODE BLOCK
            # Attempt to join current sequence with itself, this is possible
            if ( seq1.canJoin( seq1 ) ):
                appendSeqAndCacheSupport( C, seq1.join(seq1), seq1.getMis(), ctx.seqDB )
            # END UNSUPPORTED CODE BLOCK
            
            for idx2 in range( idx1+1, len( Fprev ) ):
                seq2 = Fprev[ idx2 ]
                if ( seq1.canJoin( seq2 ) ):
                    appendSeqAndCacheSupport( C, seq1.join(seq2), min( seq1.getMis(), ctx.misMap[ seq2.getLastItemId() ] ), ctx.seqDB )
                # BEGIN UNSUPPORTED CODE BLOCK
                if ( seq2.canJoin( seq1 ) ):
                    appendSeqAndCacheSupport( C, seq2.join(seq2), min( seq2.getMis(), ctx.misMap[ seq1.getLastItemId() ] ), ctx.seqDB )
                # END UNSUPPORTED CODE BLOCK

def printFreqSeqs( FHist ):
    for idxK in range( 0, len( FHist ) ):
        print ( "The number of length ", (idxK+1), " sequential patterns is ", len( FHist[idxK] ) )
        for seq in FHist[ idxK ]:
            print( seq )
        print()
  
# Main body of MS-GSP
def MSGSPMain():       
    ctx = Context()                      # Structure for storing "global" parameters
    dataPath = "../../../Data/data.txt"  # The path to the input data @TODO: Read from arguments!
    paramPath = "../../../Data/para.txt" # The path to the MIS parameter data @TODO: Read from arguments!

    # Initialize logging
    initLogger()
    
    # Load sequence data from file
    loadData( ctx.seqDB, dataPath )
    logging.getLogger("MSGSPMain").info("Loaded seqDB: " + str(ctx.seqDB))
    
    # Initialize MIS values and support difference constraint threshold @TODO: support SDC!
    ctx.sdc = loadParams(ctx.misMap, paramPath)
    logging.getLogger("MSGSPMain").info("sdc: " + str(ctx.sdc))
    
    # Sort data according to MIS values
    sortData( ctx.seqDB, ctx.misMap )
    logging.getLogger("MSGSPMain").info("Sorted seqDB: " + str(ctx.seqDB))
    
    # Generate all frequent 1-sequences
    CHist = [[]] # used for generating candidate 2-sequences
    FHist = [[]] # the set of frequent 1-sequences
    initPass( CHist[0], FHist[0], ctx )
    logging.getLogger("MSGSPMain").info("Frequent 1-sequences: " + str(FHist[0])) 
    
    # Generate candidate 2-sequences
    CHist.append([])
    level2CandidateGen( CHist[1], CHist[0], ctx )
    logging.getLogger("MSGSPMain").info("Candidate 2-sequences: " + str(CHist[1]))
    
    # Obtain all frequent 2-sequences
    FHist.append([])
    extractAllSeqsWhichSatisfyTheirMis( FHist[-1], CHist[-1] )
    logging.getLogger("MSGSPMain").info("Frequent 2-sequences: " + str(FHist[1]))
    
    # The max length of frequent k-sequences to obtain
    maxK = 3 # @TODO: Read from arguments/load from file
 
    # Generate remaining k-sequences   
    for idxK in range( 2, maxK ):
        assert( len( CHist) == idxK )
        assert( len( FHist) == idxK )
        CHist.append([])
        # Generate candidate k-sequences
        MSCandidateGenSPM( CHist[-1], FHist[-1], ctx )
        FHist.append([])
        extractAllSeqsWhichSatisfyTheirMis( FHist[-1], CHist[-1] )
        
    # Output frequent sequences
    printFreqSeqs( FHist )
    
#### Application entry point

if __name__ == '__main__':
    MSGSPMain()
