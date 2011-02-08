'''
Created on Feb 5, 2011

@author: alanperezrathke
@author: garyturovsky
@author: alanperezrathke
'''

import logging
import re
import sys

# Structure for easier passing around of "global" parameters
class Context:
    # Constructor
    def __init__(self, dataPath, paramPath):
        logging.getLogger("Context").info("Creating new context")
        self.rawSeqDB = []             # The sequence database
        self.misMap = {}               # A map of form item id -> minimum item support
        self.supportMap = {}           # A map of form item id -> actual support
        self.sdc = sys.float_info.max  # The maximum support difference constraint allowed between two sequences      
        
        loadData( self.rawSeqDB, dataPath )
        logging.getLogger("Context").info("Loaded seqDB: " + str(self.rawSeqDB))
        
        self.sdc = loadParams( self.misMap, paramPath)
        logging.getLogger("Context").info("sdc: " + str(self.sdc))
        
        sortData(self.rawSeqDB, self.misMap)
        logging.getLogger("Context").info("Sorted seqDB: " + str(self.rawSeqDB))

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