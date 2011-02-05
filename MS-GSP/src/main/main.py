'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import copy
import logging
import math
import re
import sys

####

# Structure for easier passing around of "global" parameters
class Context:
    # Constructor
    def __init__(self):
        self.rawSeqDB = []             # The sequence database
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
def loadData( rawSeqDB, fileName ):
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        strRawSeq = re.split( r"<{|}{|}>", line )
        lstRawSeq = []
        for strTrans in strRawSeq:
            if ( strTrans != '' ):
                lstTrans = []
                lstTransItems = strTrans.split( ", " ) 
                for transItem in lstTransItems:
                    lstTrans.append( int(transItem) )
                lstRawSeq.append( lstTrans )
        rawSeqDB.append( lstRawSeq )

# Loads parameters for minimum item support and support difference constraint from data file
def loadParams( misMap, fileName ):
    sdc = sys.float_info.max
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        param = re.findall(r"\d+\.*\d*",line)
        # Make sure we're dealing with the two types of parameters that we can handle
        assert( ( len( param ) == 2 and line.startswith("MIS") ) or ( len( param ) == 1 and line.startswith( "SDC" ) ) )
        # Specifying minimum support for an item
        if(line.startswith("MIS")):
            misMap[ int(param[0]) ] = float(param[1]) # mapping item id -> MIS
        # Specifying support difference constraint
        else:
            sdc = float(param[0])
    return sdc

# Sorts transactions in parameter sequence database by user supplied MIS values
def sortData( rawSeqDB, misMap ):
    for rawSeq in rawSeqDB:
        for trans in rawSeq:
            trans.sort(key=lambda x:misMap[x])

#### Sequence utilities

# Returns True if seqA contains seqB, False otherwise
def rawSeqContains( rawSeqSup, rawSeqSub ):
    if ( len( rawSeqSub ) == 0 ):
        return True
    idxSub = 0
    for idxSup in range( 0, len( rawSeqSup ) ):
        if ( set( rawSeqSub[ idxSub ] ).issubset( set( rawSeqSup[ idxSup ] ) ) ):
            idxSub += 1
            if ( idxSub == len( rawSeqSub ) ):
                return True
    return False

# Structure for representing a sequence object
class Sequence:
    # Constructor
    def __init__( self, rawSeq=[], mis=-1.0, count=-1.0, support=-1.0 ):
        self.rawSeq = rawSeq
        self.mis = mis
        self.count = count
        self.support = support
    
    # String representation
    def __repr__(self):
        strSeq = "<"
        for trans in self.getRawSeq():
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
    
    # Accessor for internal raw sequence
    def getRawSeq(self):
        return self.rawSeq
    
    # Returns total number of items in the sequence
    def length(self):
        length=0
        for trans in self.getRawSeq():
            length += len( trans )
        return length
    
    # Returns total number of transaction in the sequence
    def size(self):
        return len(self.getRawSeq())
   
    # Computes and caches support for this sequence
    # Also, of note: the optimal way to avoid reads from the database is to compute the sequence counts for a single record in the database
    # (i.e . for db for seqs instead of for seqs for db as is the case here)
    def cacheSupport(self, rawSeqDB):
        self.count = 0.0
        for rawSeq in rawSeqDB:
            if ( rawSeqContains( rawSeq, self.getRawSeq() ) ):
                self.count += 1.0
        self.support = self.count / float(len(rawSeqDB))
        
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
        assert( (len( self.getRawSeq() ) > 0) and (len( self.getRawSeq()[0] ) > 0 ) )
        return self.getRawSeq()[0][0]
        
    # Returns item id of last item in sequence
    def getLastItemId(self):
        assert( ( len( self.getRawSeq() ) > 0) and ( len( self.getRawSeq()[0] ) > 0 ) )
        return self.getRawSeq()[-1][-1]

    # Returns True if item with parameter id contains lowest MIS and MIS is unique within the sequence, False otherwise
    def itemHasUniqueMinMis(self, itemId, misMap):
        assert( itemId in misMap )
        itemMis = misMap[ itemId ]
        itemCount = 0
        for trans in self.getRawSeq():
            for transItemId in trans:
                if ( transItemId != itemId ):
                    assert( itemId in misMap )
                    if ( itemMis >= misMap[ transItemId ] ):
                        return False
                else:
                    itemCount += 1
                    if ( itemCount > 1 ):
                        return False
        assert( itemCount == 1 )
        return True;
    
    # Returns True if first item contains lowest MIS and MIS is unique within sequence, False otherwise
    def firstItemHasUniqueMinMis(self, misMap):
        return self.itemHasUniqueMinMis( self.getFirstItemId(), misMap )
    
    # Returns True if last item contains lowest MIS and MIS is unique within sequence, False otherwise
    def lastItemHasUniqueMinMis(self, misMap):
        return self.itemHasUniqueMinMis( self.getLastItemId(), misMap )
    
    def getItemAtIdx(self, idxItemToGet):
        idxCurItem = 0
        if ( idxItemToGet < 0 ):
            idxItemToGet = self.length() + idxItemToGet
        assert( 0 <= idxItemToGet < self.length() )
        for trans in self.getRawSeq():
            for idxItemInTrans in range( len( trans ) ):
                if ( idxCurItem == idxItemToGet ):
                    return trans[ idxItemInTrans ]
                idxCurItem += 1
        assert( False ) # Should never reach here
    
    # Returns new raw sequence with item at idx missing
    def getRawSeqWithoutItemAtIdx(self, idxItemToDel):
        outRawSeq = copy.deepcopy( self.getRawSeq() )
        idxCurItem = 0
        if ( idxItemToDel < 0 ):
            idxItemToDel = self.length() + idxItemToDel
        assert( 0 <= idxItemToDel < self.length() )
        for idxTrans, trans in enumerate( outRawSeq ):
            for idxItemInTrans in range( len(trans) ):
                if ( idxCurItem == idxItemToDel ):
                    del trans[ idxItemInTrans ]
                    if ( len( trans ) == 0 ):
                        del outRawSeq[ idxTrans ]
                    assert( outRawSeq != self.getRawSeq() )
                    return outRawSeq
                idxCurItem += 1
        assert( False ) # Should never reach here!  
    
    # Returns True if removing first element from this sequence and last element from parameter sequence results in same sequence, False otherwise
    def canJoin(self, seqObj):
        # Bounds checking
        assert( ( len( self.getRawSeq() ) > 0 ) and ( len( self.getRawSeq()[0] ) > 0 ) )
        assert( (len( seqObj.getRawSeq() ) > 0) and ( len( seqObj.getRawSeq()[0] ) > 0 ) )
        return ( self.getRawSeqWithoutItemAtIdx( 0 ) == seqObj.getRawSeqWithoutItemAtIdx( -1 ) )
                
    # Returns new raw sequence with the result of joining this sequence with parameter sequence
    def join(self, seqObj, misMap):
        # Bounds checking
        assert( ( len( self.getRawSeq() ) > 0 ) and ( len( self.getRawSeq()[0] ) > 0 ) )
        assert( ( len( seqObj.getRawSeq() ) > 0 ) and (len( seqObj.getRawSeq()[0] ) > 0 ) )
        # Get item id to merge into this sequence
        itemToMerge = seqObj.getLastItemId()
        joinedRawSeq = copy.deepcopy(self.getRawSeq())
        # Determine if last item should be appended to last transaction or if a new transaction should be created
        if ( len( seqObj.getRawSeq()[-1] ) > 1 ):
            assert ( itemToMerge not in joinedRawSeq[-1] )
            joinedRawSeq[-1].append( itemToMerge )
            # sort by MIS
            joinedRawSeq[-1].sort( key = lambda itemId : misMap[ itemId ] )
        else:
            joinedRawSeq.append( [itemToMerge] )
        return joinedRawSeq

# Returns True if parameter raw sequence is not found within parameter sequence object list, False otherwise
def isUniqueRawSeqWithinList( lstSeqObjs, rawSeq ):
    for seqObj in lstSeqObjs:
        if seqObj.getRawSeq() == rawSeq:
            return False
    return True

# Appends the raw sequence as a sequence object and caches the support of the sequence
# Using this utility function because we can't overload constructors in Python
def appendSeqObjAndCacheSupport( lstSeqObjs, rawSeq, mis, rawSeqDB ):
    assert( isUniqueRawSeqWithinList( lstSeqObjs, rawSeq ) )
    lstSeqObjs.append( Sequence( rawSeq, mis ) )
    lstSeqObjs[ -1 ].cacheSupport( rawSeqDB )

# Extracts all sequences from C which have a support greater than or equal to their MIS and stores them in F   
def extractAllSeqObjsWhichSatisfyTheirMis( F, C ):
    F[:] = [ seqObj for seqObj in C if ( seqObj.getSupport() >= seqObj.getMis() ) ]

#### GSP Algorithm

# Initial pass for MS-GSP
def initPass( L, F, ctx ):
    # Objectives:
    # 1. Find item 'I' with lowest MIS that has support greater than or equal to MIS(I)
    # 2. Output L : set of items that have support greater than MIS(I)
    # 3. Output F : set of 1-sequences such that support <{J}> is greater than or equal to MIS(J), note: F is a subset of L
    
    # Count unique items in sequence database
    itemCountsMap = {}
    for rawSeq in ctx.rawSeqDB:
        uniqueItemsInRawSeq = set()
        for trans in rawSeq:
            for item in trans:
                uniqueItemsInRawSeq.add(item)
        for item in uniqueItemsInRawSeq:
            if item not in itemCountsMap:
                itemCountsMap[item]  = 1.0
            else:
                itemCountsMap[item] += 1.0
                        
    # Make possible 1-sequences
    L[:] = []
    for key in itemCountsMap.keys():
        count = itemCountsMap[ key ]
        support = count / float( len( ctx.rawSeqDB ) )
        assert( ( 0.0 <= support ) and ( 1.0 >= support ) )
        L.append( Sequence( [[key]], ctx.misMap[ key ], count, support ) )
        ctx.supportMap[key] = support # Cache item supports
    
    # Sort possible 1-sequences by MIS
    L.sort( key = lambda seqObj : seqObj.getMis() )
    
    # Determine first item with support >= MIS(item) - this item will have the lowest satisfied MIS
    idxItemWithLowestSatisfiedMIS = len( L )
    for idx in range( 0, len( L ) ):
        if ( L[idx].getSupport() >= L[idx].getMis() ):
            idxItemWithLowestSatisfiedMIS = idx
            break
    
    # Trim L of all items with support lower than item with lowest satisfied MIS
    # Note use of [:] to modify list in place rather than allocate a new list
    minGlobalSatisfiedMis = L[idxItemWithLowestSatisfiedMIS].getMis()
    L[:] = [ seqObj for seqObj in L if (seqObj.getSupport() >= minGlobalSatisfiedMis) ]
    
    # Determine frequent 1-sequences
    extractAllSeqObjsWhichSatisfyTheirMis( F, L )

# Determines candidate 2-sequences
def level2CandidateGen( C, L, ctx ):
    C[:] = []
    for idxL in range( len(L) ):
        seqObjL = L[ idxL ]
        # Assert we're working with sequences of length 1
        assert( ( len( seqObjL.getRawSeq() ) == 1 ) and ( len( seqObjL.getRawSeq()[0] ) == 1 ) )
        # See if 'l' satisfies it's own MIS
        if ( seqObjL.getSupport() >= seqObjL.getMis() ):
            lId = seqObjL.getFirstItemId()
            # Create 2-tuple <{l}{l}> - it could exist!
            appendSeqObjAndCacheSupport( C, [ [ lId ], [ lId ] ], seqObjL.getMis(), ctx.rawSeqDB )
            # Create 2-tuples with all sequences 'h' where MIS(h) >= MIS(l)  
            for idxH in range( idxL+1, len(L) ):
                seqObjH = L[ idxH ]
                # Assert sequences were pre-sorted by MIS
                assert( seqObjH.getMis() >= seqObjL.getMis() )
                # Only create tuple if sup(h) is greater than MIS(l)
                if ( ( seqObjH.getSupport() >= seqObjL.getMis() ) and ( math.fabs( seqObjH.getSupport() - seqObjL.getSupport() ) <= ctx.sdc ) ):
                    # Join seqL and seqH to create both <{l,h}> and <{l,h}>
                    hId = seqObjH.getFirstItemId()
                    assert ( hId != lId ) # assert these items are unique!
                    appendSeqObjAndCacheSupport( C, [ [ lId, hId ] ], seqObjL.getMis(), ctx.rawSeqDB )
                    appendSeqObjAndCacheSupport( C, [ [ lId ], [ hId ] ], seqObjL.getMis(), ctx.rawSeqDB )
                    # Also create 2-tuple <{h},{l}> - else how could we get it?
                    appendSeqObjAndCacheSupport( C, [ [ hId ], [ lId ] ], seqObjL.getMis(), ctx.rawSeqDB )
    
def MSCandidateGenSPM_conditionalJoinWhenFirstItemHasUniqueMinMis( seqObj1, seqObj2, C, ctx ):
    assert( (ctx.misMap[seqObj2.getLastItemId()]>ctx.misMap[seqObj1.getFirstItemId()]) and ( seqObj1.getMis() == ctx.misMap[seqObj1.getFirstItemId()] ) )
    # Enforce sdc - early out if support difference exceeds threshold
    if ( math.fabs( ctx.supportMap[ seqObj1.getItemAtIdx( 1 ) ] - ctx.supportMap[ seqObj2.getItemAtIdx( -1 ) ] ) > ctx.sdc ):
        return
    if seqObj1.getRawSeqWithoutItemAtIdx(1)==seqObj2.getRawSeqWithoutItemAtIdx(-1):
        if ( len( seqObj2.getRawSeq()[-1] ) == 1 ):
            c1 = copy.deepcopy( seqObj1.getRawSeq() )
            c1.append([seqObj2.getLastItemId()])
            appendSeqObjAndCacheSupport( C, c1, seqObj1.getMis(), ctx.rawSeqDB )
                
            if( (seqObj1.size()==2) and (seqObj1.length()==2) and (ctx.misMap[seqObj2.getLastItemId()]>ctx.misMap[seqObj1.getLastItemId()]) ):
                c2 = copy.deepcopy(seqObj1.getRawSeq())
                c2[-1].append( seqObj2.getLastItemId() )
                appendSeqObjAndCacheSupport( C, c2, seqObj1.getMis(), ctx.rawSeqDB )
                            
        elif (seqObj1.length()>2) or ((seqObj1.length()==2 and seqObj1.size()==1) and (ctx.misMap[seqObj2.getLastItemId()]>ctx.misMap[seqObj1.getLastItemId()])):
            c2 = copy.deepcopy(seqObj1.getRawSeq())
            c2[-1].append( seqObj2.getLastItemId() )
            c2[-1].sort( key = lambda itemId: ctx.misMap[ itemId ] ) # Maintain lexographic order by MIS value
            appendSeqObjAndCacheSupport( C, c2, seqObj1.getMis(), ctx.rawSeqDB )

def MSCandidateGenSPM_conditionalJoinWhenLastItemHasUniqueMinMis( seqObj1, seqObj2, C, ctx ):
    assert( (ctx.misMap[seqObj2.getLastItemId()]<ctx.misMap[seqObj1.getFirstItemId()]) and ( seqObj2.getMis() == ctx.misMap[seqObj2.getLastItemId()] ) )
    # Enforce sdc - early out if support difference exceeds threshold
    if ( math.fabs( ctx.supportMap[ seqObj1.getItemAtIdx( 0 ) ] - ctx.supportMap[ seqObj2.getItemAtIdx( -2 ) ] ) > ctx.sdc ):
        return
    if ( seqObj1.getRawSeqWithoutItemAtIdx(0) == seqObj2.getRawSeqWithoutItemAtIdx(-2) ):
        if ( len( seqObj1.getRawSeq()[0] ) == 1 ):
            c1 = [[seqObj1.getFirstItemId()]]
            c1.extend( copy.deepcopy(seqObj2.getRawSeq()) )
            appendSeqObjAndCacheSupport( C, c1, seqObj2.getMis(), ctx.rawSeqDB )
                
            if( (seqObj2.size()==2) and (seqObj2.length()==2) and (ctx.misMap[seqObj1.getFirstItemId()]<ctx.misMap[seqObj2.getFirstItemId()]) ):
                c2 = copy.deepcopy(seqObj2.seq)
                c2[0].insert( 0, seqObj1.getFirstItemId() )
                appendSeqObjAndCacheSupport( C, c2, seqObj2.getMis(), ctx.rawSeqDB )

        elif (seqObj2.length()>2) or ((seqObj2.length()==1 and seqObj2.size()==2) and (ctx.misMap[seqObj2.getFirstItemId()]>ctx.misMap[seqObj1.getFirstItemId()])):
            c2 = copy.deepcopy(seqObj2.seq)
            c2[0].insert( 0, seqObj1.getFirstItemId() )
            c2[0].sort( key = lambda itemId: ctx.misMap[ itemId ] ) # Maintain lexographic order by MIS value
            appendSeqObjAndCacheSupport( C, c2, seqObj2.getMis(), ctx.rawSeqDB ) 

# Extracts all k-sequences in C such that all k-1 subsequences are frequent based on FPrev (i.e. Fk-1)
def MSCandidateGenSPM_prune( C, FPrev, misMap ):
    CPruned = []
    for candidateSeqObj in C:
        bAreAllSubsFreq = True
        for idxItemToDel in range( candidateSeqObj.length() ):
            # @TODO: Handle case where multiple objects with same min MIS exist in candidate sequence
            if ( candidateSeqObj.getMis() != misMap[ candidateSeqObj.getItemAtIdx( idxItemToDel ) ] ):
                rawCandidateSeqObj = candidateSeqObj.getRawSeqWithoutItemAtIdx( idxItemToDel )
                bSubExists = False
                for frequentSeqObj in FPrev:
                    if ( frequentSeqObj.getRawSeq() == rawCandidateSeqObj ):
                        bSubExists = True
                        break
            
                if not bSubExists:
                    bAreAllSubsFreq = False
                    break
        
        if bAreAllSubsFreq:
            CPruned.append( candidateSeqObj )
    return CPruned

# Determines candidate k-sequences where k is not 2
def MSCandidateGenSPM( C, FPrev, ctx ):
    # Join step: create candidate sequences by joining Fk-1 with Fk-1
    # NOTE: seqObj1 joins seqObj2 and seqObj2 joins with seqObj1 iff seqObj1 = <abab...ab> and seqObj2 = <baba..ba>
    for seqObj1 in FPrev:
        bSeqObj1_FirstItemHasUniqueMinMis = seqObj1.firstItemHasUniqueMinMis( ctx.misMap )
        for seqObj2 in FPrev:
            if ( bSeqObj1_FirstItemHasUniqueMinMis and (ctx.misMap[seqObj2.getLastItemId()] > ctx.misMap[seqObj1.getFirstItemId()]) ):
                MSCandidateGenSPM_conditionalJoinWhenFirstItemHasUniqueMinMis( seqObj1, seqObj2, C, ctx )        
            elif ( seqObj2.lastItemHasUniqueMinMis( ctx.misMap ) and (ctx.misMap[seqObj2.getLastItemId()] < ctx.misMap[seqObj1.getFirstItemId()]) ):
                MSCandidateGenSPM_conditionalJoinWhenLastItemHasUniqueMinMis( seqObj1, seqObj2, C, ctx )  
            elif ( seqObj1.canJoin( seqObj2 ) and ( ctx.sdc > math.fabs( ctx.supportMap[ seqObj1.getItemAtIdx(0) ] - ctx.supportMap[ seqObj2.getItemAtIdx(-1) ] ) ) ):
                appendSeqObjAndCacheSupport( C, seqObj1.join(seqObj2, ctx.misMap), min( seqObj1.getMis(), ctx.misMap[ seqObj2.getLastItemId() ] ), ctx.rawSeqDB )
    # Prune any candidate sets if all their k-1 subsets are not frequent (with the exception of the subset missing the item with the lowest mis)
    C[:] = MSCandidateGenSPM_prune( C, FPrev, ctx.misMap )

def printFreqSeqObjs( FHist ):
    for idxK in range( len( FHist ) ):
        print ( "The number of length ", (idxK+1), " sequential patterns is ", len( FHist[idxK] ) )
        for seqObj in FHist[ idxK ]:
            print( seqObj )
        print()
  
# Main body of MS-GSP
def MSGSPMain():       
    ctx = Context()                      # Structure for storing "global" parameters
    dataPath = "../../../Data/data.txt"  # The path to the input data @TODO: Read from arguments!
    paramPath = "../../../Data/para.txt" # The path to the MIS parameter data @TODO: Read from arguments!

    # Initialize logging
    initLogger()
    
    # Load sequence data from file
    loadData( ctx.rawSeqDB, dataPath )
    logging.getLogger("MSGSPMain").info("Loaded seqDB: " + str(ctx.rawSeqDB))
    
    # Initialize MIS values and support difference constraint threshold @TODO: support SDC!
    ctx.sdc = loadParams(ctx.misMap, paramPath)
    logging.getLogger("MSGSPMain").info("sdc: " + str(ctx.sdc))
    
    # Sort data according to MIS values
    sortData( ctx.rawSeqDB, ctx.misMap )
    logging.getLogger("MSGSPMain").info("Sorted seqDB: " + str(ctx.rawSeqDB))
    
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
    extractAllSeqObjsWhichSatisfyTheirMis( FHist[-1], CHist[-1] )
    logging.getLogger("MSGSPMain").info("Frequent 2-sequences: " + str(FHist[1]))
    
    # The max length of frequent k-sequences to obtain
    maxK = 10 # @TODO: Read from arguments/load from file
 
    # Generate remaining k-sequences   
    for idxK in range( 2, maxK ):
        assert( len( CHist) == idxK )
        assert( len( FHist) == idxK )
        CHist.append([])
        # Generate candidate k-sequences
        MSCandidateGenSPM( CHist[-1], FHist[-1], ctx )
        logging.getLogger("MSGSPMain").info("Candidate " + str(idxK+1) + "-sequences: " + str(CHist[-1]))
        FHist.append([])
        extractAllSeqObjsWhichSatisfyTheirMis( FHist[-1], CHist[-1] )
        logging.getLogger("MSGSPMain").info("Frequent " + str(idxK+1) + "-sequences: " + str(FHist[-1]))
    
    # Output frequent sequences
    printFreqSeqObjs( FHist )
    
#### Application entry point

if __name__ == '__main__':
    MSGSPMain()
