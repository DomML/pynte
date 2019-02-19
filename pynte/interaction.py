# from copy import deepcopy as dcp
from numpy import matrix
import numpy as np
# from scipy.spatial.distance import euclidean
from math import pi
import distanceMatrix as dmtx
from numpy import sqrt

def dcp(myMap):
    """
    'Soft' copy, named dcp because I'm lazy.
    """
    pass

    newMap = {}
    for res1 in myMap:
        newMap[res1] = {}
        for res2 in myMap[res1]:
            newMap[res1][res2] = {}
            for at1 in myMap[res1][res2]:
                newMap[res1][res2][at1] = {}
                for at2 in myMap[res1][res2][at1]:
                    pass
                    newMap[res1][res2][at1][at2] = myMap[res1][res2][at1][at2]
    
    return newMap

# np.set_printoptions(precision=2, linewidth=500)

def intraHydrophobicInteraction(distMat, cutOffDist = 4.2, cutoffCount = 5):
    """
    """
    sqCutOff = cutOffDist * cutOffDist
    notHybRes = []
    notHybAt = []
    
#     print len(distMat)
#     print sorted([i.get_id()[1] for i in distMat.keys()])
#     print distMat
    tmp = {}
    for i in distMat:
        pass
        for j in distMat[i]:
            tmp[j] = 0
#     print len(tmp)
#     print sorted([i.get_id()[1] for i in tmp.keys()])
    
    hybMat = dcp(distMat)#{}
    
#     for resI in distMat:

    
    # List non-hydrophobic residues
    for res in hybMat:
        
        hasHyb = False
        for at in res:
            if at.get_occupancy() > 0:
                hasHyb = True
                continue

        if hasHyb == False:
            notHybRes.append(res)
  
      
    # Remove residues that does not contain hydrophobics atoms
    for i in notHybRes:
        del hybMat[i]
    
    for i in hybMat.keys():
        for j in notHybRes:
            if hybMat[i].has_key(j):
                del hybMat[i][j]
  
    # Remove non-hydrophobics atoms, thus interaction involving at least
    # one non-hydrophobic atom
    for res in hybMat:
        for at in res:
            if at.get_occupancy() <= 0:
                notHybAt.append(at)
           
        for res2 in hybMat[res]:
            for at in res2:
                if at.get_occupancy() <= 0:
                    notHybAt.append(at)
    
    for i in hybMat:
        for j in hybMat[i]:
                
            for at in notHybAt:
                if at in hybMat[i][j]:
                    del hybMat[i][j][at]
                        
            for at in hybMat[i][j]:
                for at2 in notHybAt:
                    if at2 in hybMat[i][j][at]:
                        del hybMat[i][j][at][at2]
 
    # Adjust atomic distance
    atomsPairToDel = []
    for i in hybMat:
        for j in hybMat[i]:
            for at1 in hybMat[i][j]:
                for at2 in hybMat[i][j][at1]:
                    if hybMat[i][j][at1][at2] > sqCutOff:
                        atomsPairToDel.append((at1, at2))
       
    for i in atomsPairToDel:
        res1 = i[0].get_parent()
        res2 = i[1].get_parent()
        del hybMat[res1][res2][i[0]][i[1]]
 
    # Remove interactions between residues involving less than X interactions
    # By default, 6 (between CH2 and CH3)
    for resI in hybMat:
        resToDel = []
        for resJ in hybMat[resI]:
            count = 0
            for atI in hybMat[resI][resJ]:
                for atJ in hybMat[resI][resJ][atI]:
                    count +=1
     
            if count < cutoffCount:
                resToDel.append(resJ)
              
        for i in resToDel:
            del hybMat[resI][i]


    return hybMat

def intraHydrogenBond(distMat, cutOffDist = 4.0, HhH_angle = 40):
    """
    """
    HBMat = dcp(distMat)
    newHBMat = {}
    HBat = {}
    intraDmx = dmtx.intraResidueDist(HBMat.keys())

    sqCutOff = cutOffDist ** 2
    heavyAt = ["N", "O", "F"]
    heavyPairs = {}
    HydrogenHeavyPairs = {}
    resDic = {}
    
    # List of heavy atoms
    heavyDic = {}
    hydrogenDic = {}
    for res in HBMat:
        for atInd, at in enumerate(res):
            if at.element in heavyAt:
                heavyDic[at] = atInd
                resDic[at.get_parent()] = None
            elif at.element == "H":
                hydrogenDic[at] = atInd
#     print len(heavyDic)


    # List of Hydrogens close to an heavy atom
    for i in hydrogenDic:
        resI = i.get_parent()
        for j in heavyDic:
            resJ = j.get_parent()
            # Same residue
            if resI == resJ:
                dij = intraDmx[resI][i][j]
                if dij < 1.3**2:
#                     if not HeavyHydrogenPairs.has_key(j):
#                         HeavyHydrogenPairs[j] = []
                    HydrogenHeavyPairs[j] = i
                 
     
    # Find close enough heavy atoms, including one with bound hydrogen
    for heavy in heavyDic:
        resHeavy = heavy.get_parent()
        # Loop on heavy's bound to hydrogen
        for heavyH in HydrogenHeavyPairs:
            resHeavyH = heavyH.get_parent()
            hydrogen = HydrogenHeavyPairs[heavyH]
            
            # First check if resHeavy and resHeavyH CA's are close enough
            if HBMat[resHeavy].has_key(resHeavyH):
                # And then if heavy's are close too
                if HBMat[resHeavy][resHeavyH].has_key(heavy):
                    if HBMat[resHeavy][resHeavyH][heavy].has_key(heavyH):
                        if HBMat[resHeavy][resHeavyH][heavy][heavyH] < sqCutOff:
                            # And if not bind heavy and hydrogene are too
                            if HBMat[resHeavy][resHeavyH][heavy][hydrogen] < sqCutOff:
                                pass
                            
                                # Compute angle heavy--Hydrogen-heavy
                                vec_heavyH_hyd = hydrogen.get_coord() - heavyH.get_coord()
                                vec_hyd_heavy = hydrogen.get_coord() - heavy.get_coord()
                                
                                # Vector lenght
                                len_heavyH_hyd = intraDmx[resHeavyH][heavyH][hydrogen]
                                len_hyd_heavy = HBMat[resHeavy][resHeavyH][heavy][hydrogen]
                                # Becareful, lenght are squared in distance matrix !
                                len_heavyH_hyd = sqrt(len_heavyH_hyd)
                                len_hyd_heavy = sqrt(len_hyd_heavy)

                                # Dot product (OK)
                                dotP = np.dot(vec_heavyH_hyd, vec_hyd_heavy)
                                angle = dotP / (len_heavyH_hyd * len_hyd_heavy)
                                if (angle < -1):
                                	angle = -0.999999
                               	elif (angle > 1):
                               		angle = 0.999999
#                                 print angle
                                angle = np.arccos(angle)
                                angle = angle * 180 / pi
                                
                                if angle < HhH_angle:
                                    if not resHeavy in newHBMat:
                                        newHBMat[resHeavy] = {}
                                    if not resHeavyH in newHBMat[resHeavy]:
                                        newHBMat[resHeavy][resHeavyH] = {}
                                    if not heavy in newHBMat[resHeavy][resHeavyH]:
                                        newHBMat[resHeavy][resHeavyH][heavy] = {}
                                    newHBMat[resHeavy][resHeavyH][heavy][heavyH] = \
                                            HBMat[resHeavy][resHeavyH][heavy][heavyH]
    return newHBMat
    
def intraSaltBridge(distMat, cutOffDist = 4.0, include_histidine = False):
    """
    { Positive : { Negative : { at+ : { at- : value }}}}
    """
    sqCutOff = cutOffDist ** 2
    sbMat = dcp(distMat)
    # Electro +/- atoms and corresponding AA
    positive_Res = {"ARG" : None, "LYS" : None}
    positiveAtoms = {"NH1":None, "NH2":None, "NE":None, "NZ":None}
    if include_histidine:
        positive_Res["HIS"] = None
        positiveAtoms = {"NH1":None, "NH2":None, "NE":None, "NZ":None,
                         "CG":None, "ND1":None, "CD2":None, "CE1":None, "NE2":None}
 
    negative_Res = {"ASP" : None, "GLU" : None}
    negativeAtoms = {"OD1":None, "OD2":None, "OE1":None, "OE2":None}
     
    resToAtom = dict(positive_Res.items()
                     + negative_Res.items())


    # Delete useless residues
    resToDel = []
    for res in sbMat:
        if res.get_resname() not in resToAtom:
            resToDel.append(res)

    for i in resToDel:
        del sbMat[i]
     
    for i in sbMat.keys():
        for j in resToDel:
            if sbMat[i].has_key(j):
                del sbMat[i][j]


    # Select positive residue as key and negative as value
    notPosResToDel = []
    notNegResToDel = {}

    # Delete **NON** positive Residues
    for res in sbMat:
        if res.get_resname() not in positive_Res:
            notPosResToDel.append(res)
     
    for i in notPosResToDel:
        del sbMat[i]
         
    # Delete **NON** negative Residues
    for res in sbMat:
        for res2 in sbMat[res]:
            if res2.get_resname() not in negative_Res:
                notNegResToDel[res2] = None
     
    for res in sbMat:
        for resToDel in notNegResToDel.keys():
            if sbMat[res].has_key(resToDel):
                del sbMat[res][resToDel]

    # Delete Non-positive atoms
    emptyPosToDel = []
    for i in sbMat:
        for j in sbMat[i]:
            notPosAt = []
            pass
            
            for at1 in sbMat[i][j]:
                if at1.get_id() not in positiveAtoms:
                    notPosAt.append(at1)
 
            for at1 in notPosAt:
                del sbMat[i][j][at1]
                if len(sbMat[i][j]) == 0:
                    emptyPosToDel.append(i)
    
    for i in set(emptyPosToDel):
        del sbMat[i]

    # Delete Non-negative atoms
    for i in sbMat:
        emptyNegToDel = []
        for j in sbMat[i]:
            for at1 in sbMat[i][j]:
                notNegAt = []
                pass

                for at2 in sbMat[i][j][at1]:
                    if at2.get_id() not in negativeAtoms:
                        notNegAt.append(at2)

                for at2 in notNegAt:
                    del sbMat[i][j][at1][at2]
                    if len(sbMat[i][j][at1]) == 0:
                        emptyNegToDel.append(j)


        for j in set(emptyNegToDel):
            del sbMat[i][j]
    
    # Adjust atomic distance
    atomsPairToDel = []
    for i in sbMat:
        for j in sbMat[i]:
            for at1 in sbMat[i][j]:
                for at2 in sbMat[i][j][at1]:
                    if sbMat[i][j][at1][at2] > sqCutOff:
                        atomsPairToDel.append((at1, at2))
     
    for i in atomsPairToDel:
        res1 = i[0].get_parent()
        res2 = i[1].get_parent()
        del sbMat[res1][res2][i[0]][i[1]]

#     cleanDicMat(sbMat)
    return sbMat