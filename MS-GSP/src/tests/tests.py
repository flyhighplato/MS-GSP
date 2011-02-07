'''
Created on Feb 5, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import unittest
import copy
import sys
import main

from main.Context import Context
from main.Sequence import rawSeqContains
from main.main import MSGSPMain

class TestMSGSPOutput(unittest.TestCase):

    def setUp(self):
        dataPath = "../../../Data/data.txt"  # The path to the input data @TODO: Read from arguments!
        paramPath = "../../../Data/para.txt" # The path to the MIS parameter data @TODO: Read from arguments!

        self.ctx = Context(dataPath,paramPath)
        self.maxK=3
    
    def getSupport(self,item):
        support=0
        for rawSeq in self.ctx.rawSeqDB:
            if(rawSeqContains(rawSeq,[[item]])):
                support+=1
        return support
    
    # Output discrepancies between the brute force output and the algorithm
    def reportDiscrepancies(self,nextSeqs,FHist,k):
        nextSeqsCopy=copy.deepcopy(nextSeqs)
        print("\r\n -- " + str(k+1) + "-sequences --")
        countMissing=0
        countFound=0
        
        # Check for missing sequences
        for idx,seq in enumerate(copy.deepcopy(nextSeqsCopy)):
            freqCount=0
            containingSeqs=[]
            for rawSeq in self.ctx.rawSeqDB:
                if(rawSeqContains(rawSeq,seq)):
                    containingSeqs.append(rawSeq)
                    freqCount+=1
            actualSupport = freqCount/len(self.ctx.rawSeqDB)
            minMIS=sys.float_info.max
            maxSupport=sys.float_info.min
            minSupport=sys.float_info.max
            
            
            for trans in seq:
                for item in trans:
                    if(self.ctx.misMap[item]<minMIS):
                        minMIS=self.ctx.misMap[item]
                    sup=self.getSupport(item)
                    if(sup>maxSupport):
                        maxSupport=sup
                    if(sup<minSupport):
                        minSupport=sup
                    
            maxSupport=maxSupport/len(self.ctx.rawSeqDB)
            minSupport=minSupport/len(self.ctx.rawSeqDB)
                        
            if(minMIS<=actualSupport and maxSupport-minSupport<=self.ctx.sdc):
                
                bIsInFHist=False
                for fseq in FHist[k]:
                    if(seq==fseq.getRawSeq()):
                        bIsInFHist=True
                        break
                    
                if(not bIsInFHist):
                    countMissing+=1
                    print("MISSING:",seq,":\r\n\tMIS ==> ",minMIS,"<",actualSupport,",\r\n\tSDC ==> (",maxSupport,"-",minSupport,"=",(maxSupport-minSupport),") <",self.ctx.sdc," from: ")
                    for idx,contSeqs in enumerate(containingSeqs):
                        print("\t\t",str(idx+1) + ".",contSeqs)
                else:
                    countFound+=1
            else:
                    
                nextSeqsCopy.remove(seq)
        
        countIncorrect=0
        
        # Check for incorrect sequences
        for fseq in FHist[k]:
            bIsInNextSeqs=False
            for seq in nextSeqsCopy:
                if(seq==fseq.getRawSeq()):
                    bIsInNextSeqs=True
                    break
            if(not bIsInNextSeqs):
                freqCount=0
                containingSeqs=[]
                for rawSeq in self.ctx.rawSeqDB:
                    if(rawSeqContains(rawSeq,fseq.getRawSeq())):
                        containingSeqs.append(rawSeq)
                        freqCount+=1
                actualSupport = freqCount/len(self.ctx.rawSeqDB)
                minMIS=sys.float_info.max
                maxSupport=sys.float_info.min
                minSupport=sys.float_info.max
                
                for trans in fseq.getRawSeq():
                    for item in trans:
                        if(self.ctx.misMap[item]<minMIS):
                            minMIS=self.ctx.misMap[item]
                        sup=self.getSupport(item)
                    if(sup>maxSupport):
                        maxSupport=sup
                    if(sup<minSupport):
                        minSupport=sup
                
                maxSupport=maxSupport/len(self.ctx.rawSeqDB)
                minSupport=minSupport/len(self.ctx.rawSeqDB)
                countIncorrect+=1
                print("INCORRECT?:",fseq.getRawSeq(),":\r\n\tMIS ==> ",minMIS,"<",actualSupport,",\r\n\tSDC ==> (",maxSupport,"-",minSupport,"=",(maxSupport-minSupport),") <",self.ctx.sdc," found in: ")
                for idx,contSeqs in enumerate(containingSeqs):
                        print("\t\t",str(idx+1) + ".",contSeqs)
                
        print(countFound," items correct")
        print(countIncorrect," items incorrect")
        print(countMissing," items missing")
        self.assertTrue(countMissing==0 and countIncorrect==0,"There are missing and/or incorrect items!")
        
    def test_BruteForce(self):
        # Get history of frequent sequences from the algorithm
        FHist=MSGSPMain(self.maxK)

        # Grab all unique items in our sequence DB
        uniqueItems=set()
        for rawSeq in self.ctx.rawSeqDB:
            for trans in rawSeq:
                for item in trans:
                    uniqueItems.add(item)
        
        # Build 1-sequence collection
        nextSeqs=[]
        for item in uniqueItems:
            nextSeqs.append( [[item]] )
        
        self.reportDiscrepancies(nextSeqs,FHist,0)
        
        # Build k-sequence collections and test each one
        for k in range(1,self.maxK):
            prevSeqs=copy.deepcopy(nextSeqs)
            nextSeqs=[]
            for seq in prevSeqs:         
                for item in uniqueItems:
                    if(item not in seq[-1]):
                        c=copy.deepcopy(seq)
                        c[-1].append(item)
                        c[-1].sort(key=lambda x: self.ctx.misMap[x])
                        nextSeqs.append(c)
                    
                    c=copy.deepcopy(seq)
                    c.append([item])
                    nextSeqs.append(c)
   
            self.reportDiscrepancies(nextSeqs,FHist,k)

if __name__ == '__main__':
    unittest.main()