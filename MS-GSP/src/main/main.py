'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import re

#### Global variables

# The sequence database
seqDB = []

# The parameters
misMap = {}

# The path to the input data
dataPath = "../../../Data/data.txt"

# The path to the MIS parameter data
paramPath = "../../../Data/para.txt"

# Support difference constraint
sdc = 0.0

#### I/O

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

def loadParam( map, fileName):
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        param = re.findall(r"\d+\.*\d*",line)
        
        #Make sure we're dealing with the two types of parameters that we can handle
        assert( ( len( param ) == 2 and line.startswith("MIS") ) or ( len( param ) == 1 and line.startswith( "SDC" ) ) )
        
        #Specifying minimum support for an item
        if(line.startswith("MIS")):
            map[int(param[0])] = float(param[1])
        #Specifying support difference constraint
        else:
            sdc = float(param[0])
        
if __name__ == '__main__':
    loadData( seqDB, dataPath )
    loadParam(misMap,paramPath)
    
    print( seqDB )
    print( misMap ) 
    print( sdc )
