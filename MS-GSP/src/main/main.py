'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import copy
import logging
import math


    
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


# Returns True if the support difference between items at parameter indices of two sequences satisfies the sdc, False otherwise
def satisfiesSDC( seqObj1, idxItem1, seqObj2, idxItem2, ctx ):
    return ( math.fabs( seqObj1.getSupportForItemAtIdx( idxItem1, ctx.supportMap ) - seqObj2.getSupportForItemAtIdx( idxItem2, ctx.supportMap ) ) <= ctx.sdc )

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
                if((seqObjL.getFirstItemId()==13 and seqObjH.getFirstItemId()==2) or (seqObjL.getFirstItemId()==2 and seqObjH.getFirstItemId()==13)):
                    print("FOUND")
                # Assert sequences were pre-sorted by MIS
                assert( seqObjH.getMis() >= seqObjL.getMis() )
                # Only create tuple if sup(h) is greater than MIS(l)
                if ( ( seqObjH.getSupport() >= seqObjL.getMis() ) and satisfiesSDC( seqObjL, 0, seqObjH, 0, ctx ) ):
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
    if not satisfiesSDC( seqObj1, 1, seqObj2, -1, ctx ):
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
    if not satisfiesSDC( seqObj1, 0, seqObj2, -2, ctx ):
        return
    if ( seqObj1.getRawSeqWithoutItemAtIdx(0) == seqObj2.getRawSeqWithoutItemAtIdx(-2) ):
        if ( len( seqObj1.getRawSeq()[0] ) == 1 ):
            c1 = [[seqObj1.getFirstItemId()]]
            c1.extend( copy.deepcopy(seqObj2.getRawSeq()) )
            appendSeqObjAndCacheSupport( C, c1, seqObj2.getMis(), ctx.rawSeqDB )
                
            if( (seqObj2.size()==2) and (seqObj2.length()==2) and (ctx.misMap[seqObj1.getFirstItemId()]<ctx.misMap[seqObj2.getFirstItemId()]) ):
                c2 = copy.deepcopy(seqObj2.getRawSeq())
                c2[0].insert( 0, seqObj1.getFirstItemId() )
                appendSeqObjAndCacheSupport( C, c2, seqObj2.getMis(), ctx.rawSeqDB )

        elif (seqObj2.length()>2) or ((seqObj2.length()==1 and seqObj2.size()==2) and (ctx.misMap[seqObj2.getFirstItemId()]>ctx.misMap[seqObj1.getFirstItemId()])):
            c2 = copy.deepcopy(seqObj2.getRawSeq())
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
            elif ( seqObj1.canJoin( seqObj2 ) and satisfiesSDC( seqObj1, 0, seqObj2, -1, ctx ) ):
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
def MSGSPMain(maxK = 10):                             # Structure for storing "global" parameters
    dataPath = "../../../Data/data.txt"  # The path to the input data @TODO: Read from arguments!
    paramPath = "../../../Data/para.txt" # The path to the MIS parameter data @TODO: Read from arguments!

    
    
    ctx = Context(dataPath,paramPath)
    
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
    #printFreqSeqObjs( FHist )
    return FHist
    
#### Application entry point

if __name__ == '__main__':
    from Sequence import Sequence,isUniqueRawSeqWithinList
    from Context import Context
    
    # Initialize logging
    initLogger()
    FHist=MSGSPMain(2)
    printFreqSeqObjs( FHist )
else:
    from main.Sequence import Sequence,isUniqueRawSeqWithinList
    from main.Context import Context
