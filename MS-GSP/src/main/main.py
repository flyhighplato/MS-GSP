'''
Created on Jan 26, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import re

#### Global variables

# The sequence database
seqDB = []

# The path to the input data
dataPath = "../../../Data/data.txt"

#### I/O

def loadData( db, fileName ):
    print( "opening: ", fileName )
    FILE = open( fileName, "r" )
    for line in FILE:
        line = line.rstrip('\n')
        seqStr = re.split( r"<{|}{|}>", line )
        print( line )
        seqLst = []
        for transStr in seqStr:
            if ( transStr != '' ):
                transLst = []
                transItemsLst = transStr.split( ", " ) 
                for transItem in transItemsLst:
                    transLst.append( int(transItem) )
                seqLst.append( transLst )
        db.append( seqLst )
        
if __name__ == '__main__':
    loadData( seqDB, dataPath )
    print ( "-------------")
    print ( seqDB )
