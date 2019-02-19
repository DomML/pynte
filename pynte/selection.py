################################################################################
#    Contain class that create protein domains to be used for representations
#        Interactions between residues in a same domain will no be shown
#        Interactions involving residues not in a domain will not be shown too
################################################################################
from operator import __getitem__
from warnings import catch_warnings
from Bio.PDB.Polypeptide import *

def residue_filter(anyMap, resNameFilter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """
    """
    for resI in list(anyMap):
#         print "##", three_to_one(resI.get_resname())
        if three_to_one(resI.get_resname()) not in resNameFilter:
            del anyMap[resI]
        else:
            for resJ in list(anyMap[resI]):
                if three_to_one(resJ.get_resname()) not in resNameFilter:
                    del anyMap[resI][resJ]

    return anyMap

class Domain(object):
    """
    """
    def __init__(self, Chain, firstOfDomain=-1, lastOfDomain=-1):
        """
        """
        pass
        
        self.chain=Chain
        self.domainBoundaries={}
        self.domainBoundariesInverted = {}
        self.firstOfChain=-1
        self.lastOfChain=-1
        self.domains = {}
        
        self.firstOfChain, self.lastOfChain = \
            [i.get_id()[1] for i in self.chain.get_residues()] \
            [0::len(self.chain)-1]
    
        for i in self.chain.get_residues():
            self.domains[i] = -1
            
        if firstOfDomain != -1:
            self.addDomain(firstOfDomain, lastOfDomain)
    
    def __contains__(self, item):
        """
        """
        pass
        
#         print self.domains.keys()
        return item in self.domains.keys()
    
    def __getitem__(self, item):
        """
        """
        pass
    
        return self.domains[item]

    def addDomain(self, firstOfDomain=-1, lastOfDomain=-1):
        """
        """
        pass
        
        if (lastOfDomain == -1 or
            lastOfDomain > self.lastOfChain):
            lastOfDomain = self.lastOfChain
        
        if firstOfDomain == -1:
            firstOfDomain = self.firstOfChain
            
        if (firstOfDomain > lastOfDomain or
            firstOfDomain > self.lastOfChain):
            return
        
        
        if len(self.domainBoundaries) > 0:
            myKey = max(self.domainBoundaries.keys()) + 1
            self.domainBoundaries[myKey] = (firstOfDomain, lastOfDomain)
        else:
            myKey = 0
            self.domainBoundaries[myKey] = (firstOfDomain, lastOfDomain)

#         print self.domains
        for i in xrange(firstOfDomain, lastOfDomain + 1):
            self.domains[self.chain[i]] = myKey

        # Invert boundaries, to get "owner"'s domain of a residue
        for i in self.domainBoundaries.keys():
            for j in xrange(self.domainBoundaries[i][0], 
                            self.domainBoundaries[i][1] + 1):
                self.domainBoundariesInverted[j] = i
        
#         print self.domainBoundaries
        return len(self.domainBoundaries)-1

    def delDomain(self, domainID):
        """
        """
        pass
        
        del self.domainBoundaries[domainID]
        
    
    def mergeDomain(self, domainID_1, domainID_2):
        """
        """
        pass
        
        if not (domainID_1 in self.domains.values() and
                domainID_2 in self.domains.values()):
            return
        
        if not (domainID_1 != domainID_2):
            return
        
        # Change ID
        myKey = max(self.domainBoundaries.keys()) + 1
        for i in self.domains.keys():
            if (self.domains[i] == domainID_1 or
                self.domains[i] == domainID_2):
                self.domains[i] = myKey
        
        self.domainBoundaries[myKey] = [self.domainBoundaries[domainID_1], 
                                        self.domainBoundaries[domainID_2]]
        del self.domainBoundaries[domainID_1]
        del self.domainBoundaries[domainID_2]
        


    def filterInteraction(self, myMap):
        """
        """
        pass

        newMap = {}
    
        try:
            for resI in myMap:
                # If this residue have to be shown
                if self[resI] != -1:
                    for resJ in myMap[resI]:
                        # If this residue have to be shown
                        if self[resJ] != -1:
                            # If residues are in differents groups
                            if self[resI] != self[resJ]:
                                if not resI in newMap:
                                    newMap[resI] = {}
                                if not resJ in newMap:
                                    newMap[resJ] = {}

                                if not resI in newMap[resJ]:
                                    newMap[resJ][resI] = {}
                                if not resJ in newMap[resI]:
                                    newMap[resI][resJ] = {}

                                newMap[resI][resJ] = myMap[resI][resJ]
                                newMap[resJ][resI] = myMap[resJ][resI]


        except :
            pass

        return newMap