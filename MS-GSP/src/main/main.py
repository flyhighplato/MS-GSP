'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import re
import sys

#### Global variables

# The sequence database
seqDB = []

# A map of form item id -> minimum item support
misMap = {}

# The path to the input data
dataPath = "../../../Data/data.txt"

# The path to the MIS parameter data
paramPath = "../../../Data/para.txt"

# Support difference constraint
sdc = 0.0

# Global minimum MIS value
globalMinMis = 0.0;

#### Initialization utilities

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
    minMis=sys.float_info.max
    
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        param = re.findall(r"\d+\.*\d*",line)
        
        # Make sure we're dealing with the two types of parameters that we can handle
        assert( ( len( param ) == 2 and line.startswith("MIS") ) or ( len( param ) == 1 and line.startswith( "SDC" ) ) )
        
        # Specifying minimum support for an item
        if(line.startswith("MIS")):
            map[int(param[0])] = float(param[1])
            if(minMis>float(param[1])):
                minMis=float(param[1])
        # Specifying support difference constraint
        else:
            sdc = float(param[0])

    return minMis

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
    def __init__(self, seq=[], count=0):
        self.seq = seq
        self.count = count
    
    # String representation
    def __repr__(self):
        return str(self.seq) + ":" + str(self.count)
    
    # Returns True if this sequence contains parameter sequence, False otherwise
    def contains(self, seq):
        return seqContains( self.seq, seq )

#### GSP Algorithm

# Initial pass for MS-GSP
def initPass( seqDB, misMap, globalMinMis, L, F ):
    T = {}
    #Count unique items in support
    for seq in seqDB:
        S = set()
        for trans in seq:
            for item in trans:
                S.add(item)
        for item in S:
            if item not in T:
                T[item]=1
            else:
                T[item]+=1
                    
    #Make possible 1-sequences
    for key in T.keys():
        if(T[key]/len(seqDB)>=globalMinMis):
            L.add(key)
    
    #Make frequent 1-sequences
    for key in T.keys():
        if(T[key]/len(seqDB)>=misMap[key]):
            F.append([Sequence(key,T[key])])

# Main body of MS-GSP
def MSGSPMain():
    loadData( seqDB, dataPath )
    globalMinMis = loadParams(misMap,paramPath)
    sortData( seqDB, misMap )
    
    L=set()
    F=[]
    initPass( seqDB, misMap, globalMinMis, L, F ) 
    
#### Application entry point

if __name__ == '__main__':
    MSGSPMain()
