'''
Created on Feb 5, 2011

@author: garyturovsky
@author: alanperezrathke
'''

import copy
import random
import sys
import unittest

from main.Context import Context
from main.Sequence import rawSeqContains
from main.main import MSGSPMain

class TestMSGSPOutput(unittest.TestCase):
    def generateInputs(self):
        dataFile = open(self.dataPath, 'w')
        paramFile = open(self.paramPath, 'w')
        
        allSeqs=[]
        uniqueItems={}
        for i in range(0,100):
            seq=[]
            seqLength = random.randint(1,self.maxSeqLength)
            while seqLength>0:
                transLen = random.randint(1,seqLength)
                trans=[]
                for j in range(0,transLen):
                    seqLength-=1
                    item = random.randint(1,self.maxUniqueItems)
                    trans.append(item)
                    if item not in uniqueItems:
                        uniqueItems[item]=random.random()
                seq.append(trans)
            allSeqs.append(seq)
        
        for seq in allSeqs:
            s = "<"
            for trans in seq:
                s+="{"
                for item in trans[0:-2]:
                    s+=str(item) + ", "
                s+=str(trans[-1]) + "}"
            s+=">\r\n"
            dataFile.write(s)
            
        for key in uniqueItems.keys():
            s = "MIS( " + str(key) + " ) = " + str(uniqueItems[key]) + "\r\n"
            paramFile.write(s)
        
        paramFile.write("SDC = " + str(random.random()) + "\r\n")
        
         
    def setUp(self):
        self.dataPath = "../../../Data/UnitTestData.txt"  # The path to the input data @TODO: Read from arguments!
        self.paramPath = "../../../Data/UnitTestPara.txt" # The path to the MIS parameter data @TODO: Read from arguments!
        self.maxUniqueItems=30
        self.maxSeqLength=30
        self.maxK=3
        
        self.generateInputs()
        print("Input files generated...")
        self.ctx = Context(self.dataPath,self.paramPath)
        print("Context created...")
        
    
    #Returns amount of sequences in support of item
    def getSupport(self,item):
        support=0
        for rawSeq in self.ctx.rawSeqDB:
            if(rawSeqContains(rawSeq,[[item]])):
                support+=1
        return support
    
    def getSeqSupport(self,seq):
        containingSeqs=[]
        for rawSeq in self.ctx.rawSeqDB:
            if(rawSeqContains(rawSeq,seq)):
                containingSeqs.append(rawSeq)
        return containingSeqs
    
    def getSeqStats(self,seq):
        class SeqStats:
            minMIS=sys.float_info.max
            maxSupport=sys.float_info.min
            minSupport=sys.float_info.max
            containingSeqs=[]
            freqCount=0
            actualSupport=0.0
        
        stats = SeqStats()
        
        stats.containingSeqs = self.getSeqSupport(seq)
        stats.freqCount=len(stats.containingSeqs) 
        stats.actualSupport = stats.freqCount/len(self.ctx.rawSeqDB)

        for trans in seq:
            for item in trans:
                if(self.ctx.misMap[item]<stats.minMIS):
                    stats.minMIS=self.ctx.misMap[item]
                sup=self.getSupport(item)
                if(sup>stats.maxSupport):
                    stats.maxSupport=sup
                if(sup<stats.minSupport):
                    stats.minSupport=sup
                
        stats.maxSupport=stats.maxSupport/len(self.ctx.rawSeqDB)
        stats.minSupport=stats.minSupport/len(self.ctx.rawSeqDB)
        
        return stats
        
    # Output discrepancies between the brute force output and the algorithm
    def reportDiscrepancies(self,nextSeqs,FHist,k):
        nextSeqsCopy=copy.deepcopy(nextSeqs)
        print("\r\n -- " + str(k+1) + "-sequences --")
        countMissing=0
        countFound=0
        
        #Check for missing sequences in the MSGSP output
        for idx,seq in enumerate(copy.deepcopy(nextSeqsCopy)):
                    
            stats=self.getSeqStats(seq)
       
            if(stats.minMIS<=stats.actualSupport and stats.maxSupport-stats.minSupport<=self.ctx.sdc):
                
                bIsInFHist=False
                for fseq in FHist[k]:
                    if(seq==fseq.getRawSeq()):
                        bIsInFHist=True
                        break
                    
                if(not bIsInFHist):
                    countMissing+=1
                    print("MISSING:",seq,":\r\n\tMIS ==> ",stats.minMIS,"<",stats.actualSupport,",\r\n\tSDC ==> (",stats.maxSupport,"-",stats.minSupport,"=",(stats.maxSupport-stats.minSupport),") <",self.ctx.sdc," from: ")
                    for idx,contSeqs in enumerate(stats.containingSeqs):
                        print("\t\t",str(idx+1) + ".",contSeqs)
                        
                else:
                    countFound+=1
            else:
                    
                nextSeqsCopy.remove(seq)
        
        #Check for incorrect sequences
        countIncorrect=0

        for fseq in FHist[k]:
            
            #Check if it exists in collection of sequences from the brute force approach
            bIsInNextSeqs=False
            for seq in nextSeqsCopy:
                if(seq==fseq.getRawSeq()):
                    bIsInNextSeqs=True
                    break
            if(not bIsInNextSeqs):
                
                #Calculate minimum MIS, maximum/minimum support for items in the sequence
                stats=self.getSeqStats(fseq.getRawSeq())
       
                #This sequence wasn't in the brute force approach so it might be incorrect (or the brute force algorithm is incorrect)
                countIncorrect+=1
                print("INCORRECT?:",fseq.getRawSeq(),":\r\n\tMIS ==> ",stats.minMIS,"<",stats.actualSupport,",\r\n\tSDC ==> (",stats.maxSupport,"-",stats.minSupport,"=",(stats.maxSupport-stats.minSupport),") <",self.ctx.sdc," found in: ")
                for idx,contSeqs in enumerate(stats.containingSeqs):
                        print("\t\t",str(idx+1) + ".",contSeqs)
        
        #Print summary
        print(countFound," items correct")
        print(countIncorrect," items incorrect")
        print(countMissing," items missing")
        self.assertTrue(countMissing==0 and countIncorrect==0,"There are missing and/or incorrect items!")
        
    def test_BruteForce(self):
    
        # Get history of frequent sequences from the algorithm
        FHist=MSGSPMain(self.maxK,self.dataPath,self.paramPath)
        
        print("Finished running MSGSPMain. Beginning verification...")

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
                    
                    stats=self.getSeqStats(c)
       
                    if(stats.minMIS<=stats.actualSupport and stats.maxSupport-stats.minSupport<=self.ctx.sdc):
                        nextSeqs.append(c)
   
            self.reportDiscrepancies(nextSeqs,FHist,k)

if __name__ == '__main__':
    unittest.main()