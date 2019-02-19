import numpy
from numpy import matrix
from sys import exit


def atomSquareDistance(at1, at2):
    """
    """
    diff = at1.coord - at2.coord
    return numpy.dot(diff, diff)

def resAtomDistance(res1, res2, sqCutOff):
    """
    """

    inContact = False
#     mat = [[-1] * res2len for i in xrange(res1len)]
    mat = {}
    for i in res1:
        for j in res2:
#             mat[i][j] = -1
            sqDist =  atomSquareDistance(i, j)
            if sqDist <= sqCutOff:
                inContact = True
                if i not in mat:
                    mat[i] = {}
                mat[i][j] = sqDist
#             else:
#                 mat[ind][jnd] = -1


#     return matrix(mat), inContact
#     print res1, res2, mat, inContact
    return mat, inContact

def translateDic(dic):
    """
    """
    tDic = {}
    
    for i in dic:
        for j in dic[i]:
            if j not in tDic:
                tDic[j] = {}
            tDic[j][i] = dic[i][j]
    return tDic

def intraMolecularDistanceMatrix(chain, cutOff= 4.0, CAdist = 20.0):
    """
    CAdist = preprossessing, to sprrdup computation.
    """
    squaredCAdist = CAdist ** 2
    sqCutOff = cutOff**2
#     print "#"
#     print cutOff, sqCutOff
#     print squaredCAdist

    resresDic = {}

    for ind, i in enumerate(chain):
        for jnd, j in enumerate(chain):
            if ind < jnd:
                # Compute the CA-CA distance
                squareCacaDist = atomSquareDistance(i["CA"], j["CA"])

                # If CA are close enough
                if squareCacaDist <= squaredCAdist:
                    # Slowing, to improve
                    atDistMat, inContact = resAtomDistance(i, j, sqCutOff)
    
                    if inContact == True:
                        if not i in resresDic.keys():
                            resresDic[i] = {}
                        if not j in resresDic.keys():
                            resresDic[j] = {}

                        resresDic[i][j] = {}#atDistMat
                        resresDic[j][i] = {}#translateDic(atDistMat)
                        
                        # Copy distance atom by atom
                        for atI in atDistMat:
                            for atJ in atDistMat[atI]:
                                if atI not in resresDic[i][j]:
                                    resresDic[i][j][atI] = {}
                                if atJ not in resresDic[j][i]:
                                    resresDic[j][i][atJ] = {}

                                resresDic[i][j][atI][atJ] = atDistMat[atI][atJ]
                                resresDic[j][i][atJ][atI] = atDistMat[atI][atJ]
                                
#             if j.get_id()[1] == 453 and i.get_id()[1] == 448:
#                 tmpi = i
#                 tmpj = j
#                 print j, i, squareCacaDist
#                 print resresDic[i][j]
#                 print "##"
#                 print resresDic[j][i]
#     print resresDic[tmpi][tmpj]

    return resresDic

def intraResidueDist(chain):
    """
    """
    resMat = {}
    
    for res in chain:
#         resMat[res] = [[-1] * len(res) for i in xrange(len(res))]
        resMat[res] = {}
        
        for ind, i in enumerate(res):
            if i not in resMat:
                resMat[res][i] = {}
 
            for jnd, j  in enumerate(res):
#                 if j not in resMat:
#                     resMat[res][j] = {}
 
#                 if ind < jnd:
                    
                resMat[res][i][j] = i-j
#                     resMat[res][j][i] = resMat[res][i][j]
#                     print resMat[res][i][j], resMat[res][j][i]
#                     print resMat[res][i][j]
#             print "--", resMat[res][i]
#         print "+", resMat[res]
    
#     print resMat
    return resMat