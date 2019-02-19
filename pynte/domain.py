################################################################################
#    domain.py :
#        * Define domains in protein sequence to be used in representation,
#        to filter interaction (remove intra-domain interactions during
#        representation process).
################################################################################
from Cython.Compiler.Nodes import PassStatNode

class Domain(object):
    """
    """
    pass

    def __init__(self, distanceMap):
        """
        """
        pass
        
        self.distanceMap = distanceMap
        self.resDomainDict = {}
    
    
    def __getattr__(self, attr):
        return self.resDomainDict[attr]
    def __setattr__(self, attr, value):
        self.resDomainDict[attr] = value
        
    def add_domain(self, id, chain, start = -1, stop = -1):
        """
        """
        pass
        
        
    def del_domain(self, id):
        """
        """
        pass
        
        
    def merge_domain(self, id1, id2):
        """
        """
        pass