'''
Created on Feb 5, 2011

@author: garyturovsky
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
        self.maxK=2
    
    #Output discrepancies between the brute force output and the algorithm
    def reportDiscrepancies(self,nextSeqs,FHist,k):
        print("\r\n -- " + str(k+1) + "-sequences --")
        countMissing=0
        countFound=0
        
        #Check for missing sequences
        for idx,seq in enumerate(copy.deepcopy(nextSeqs)):
            freqCount=0
            containingSeqs=[]
            for rawSeq in self.ctx.rawSeqDB:
                if(rawSeqContains(rawSeq,seq)):
                    containingSeqs.append(rawSeq)
                    freqCount+=1
            actualSupport = freqCount/len(self.ctx.rawSeqDB)
            minMIS=sys.float_info.max
            maxMIS=sys.float_info.min
            
            for trans in seq:
                for item in trans:
                    if(self.ctx.misMap[item]<minMIS):
                        minMIS=self.ctx.misMap[item]
                    if(self.ctx.misMap[item]>maxMIS):
                        maxMIS=self.ctx.misMap[item]
               
            if(minMIS<=actualSupport and maxMIS-minMIS<self.ctx.sdc):
                bIsInFHist=False
                for fseq in FHist[k]:
                    if(seq==fseq.getRawSeq()):
                        bIsInFHist=True
                        break

                if(not bIsInFHist):
                    countMissing+=1
                    print("MISSING:",seq,":\r\n\tMIS ==> ",minMIS,"<",actualSupport,",\r\n\tSDC ==> (",maxMIS,"-",minMIS,"=",(maxMIS-minMIS),") <",self.ctx.sdc," from: ")
                    for idx,contSeqs in enumerate(containingSeqs):
                        print("\t\t",idx+1,".",contSeqs)
                else:
                    countFound+=1
            else:
                nextSeqs.remove(seq)
        
        #Check for incorrect sequences
        for fseq in FHist[k]:
            bIsInNextSeqs=False
            for seq in nextSeqs:
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
                maxMIS=sys.float_info.min
                
                for trans in fseq.getRawSeq():
                    for item in trans:
                        if(self.ctx.misMap[item]<minMIS):
                            minMIS=self.ctx.misMap[item]
                        if(self.ctx.misMap[item]>maxMIS):
                            maxMIS=self.ctx.misMap[item]
                print("INCORRECT:",fseq.getRawSeq(),":\r\n\tMIS ==> ",minMIS,"<",actualSupport,",\r\n\tSDC ==> (",maxMIS,"-",minMIS,"=",(maxMIS-minMIS),") <",self.ctx.sdc," found in: ")
                for idx,contSeqs in enumerate(containingSeqs):
                        print("\t\t",idx+1,".",contSeqs)
                
        print(countFound," items correct")
        print(len(FHist[k])-countFound," items incorrect")
        print(countMissing," items missing")
        self.assertTrue(countMissing==0 and len(FHist[k])-countFound==0,"There are missing and/or incorrect items!")
        
        
    def test_BruteForce(self):
        #Get history of frequent sequences from the algorithm
        FHist=MSGSPMain(self.maxK)

        #Grab all unique items in our sequence DB
        uniqueItems=set()
        for rawSeq in self.ctx.rawSeqDB:
            for trans in rawSeq:
                for item in trans:
                    uniqueItems.add(item)
        
        #Build 1-sequence collection
        nextSeqs=[]
        for item in uniqueItems:
            nextSeqs.append( [[item]] )
        
        self.reportDiscrepancies(nextSeqs,FHist,0)
        
        #Build k-sequence collections and test each one
        for k in range(1,self.maxK):
            prevSeqs=copy.deepcopy(nextSeqs)
            nextSeqs=[]
            for seq in prevSeqs:
                for item in uniqueItems:
                    if(item not in seq[-1]):
                        c=copy.deepcopy(seq)
                        c[-1].append(item)
                        nextSeqs.append(c)
                    
                    c=copy.deepcopy(seq)
                    c.append([item])
                    nextSeqs.append(c)
   
            self.reportDiscrepancies(nextSeqs,FHist,k)

if __name__ == '__main__':
    unittest.main()