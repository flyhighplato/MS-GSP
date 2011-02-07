'''
Created on Feb 5, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import copy

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
    
    # Returns support for item at parameter idx
    def getSupportForItemAtIdx(self, idx, supportMap ):
        assert( self.getItemAtIdx( idx ) in supportMap )
        return supportMap[ self.getItemAtIdx( idx ) ]
    
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

# Returns True if parameter raw sequence is not found within parameter sequence object list, False otherwise
def isUniqueRawSeqWithinList( lstSeqObjs, rawSeq ):
    for seqObj in lstSeqObjs:
        if seqObj.getRawSeq() == rawSeq:
            return False
    return True