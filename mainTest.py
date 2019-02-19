import Bio.PDB as PDB
import pynte.distanceMatrix as dm
from pynte import interaction
from pynte import representation as rp
from pynte.selection import *

from numpy import matrix

import pickle, cPickle

import warnings

def main():
    warnings.filterwarnings('ignore')
    
    myPDBparser = PDB.PDBParser(PERMISSIVE=1)
    # R11-15
    myStructure = myPDBparser.get_structure("R11-15", "../Structures/R11-15.SAXSmodel.MHP.pdb")
    pdbPath = "../Structures/R11-15.SAXSmodel.MHP.pdb"
    # 3UUN
#     myStructure = myPDBparser.get_structure("3UUN", "../Structures/3UUN_A.MHP.pdb")
#     pdbPath = "../Structures/3UUN_A.MHP.pdb"
#     # R16-17_nNOS
#     myStructure = myPDBparser.get_structure("R16-17_nNOS", "../Structures/R16-17_nNOS_model1.pdb")
#     pdbPath = "../Structures/R16-17_nNOS_model1.pdb"

    
    # for chain in myStructure.get_chains():
    #     print chain
    # 
    # prev = None
    # for atom in myStructure.get_atoms():
    #     print atom.get_id()
    #     print atom.get_serial_number()
    #     print atom.get_parent().get_full_id()
    # #     print atom.get_coord()
    # #     print atom.get_vector()
    # #     print atom.is_disordered()
    #     if prev != None:
    #         print atom - prev
    #     prev = atom
    # 
    # residue = chain[20]

    
    for chain in myStructure.get_chains():
        print chain
    
        if False:
            distMat = dm.intraMolecularDistanceMatrix(chain, cutOff=5.0)
            cPickle.dump(distMat, open("myData.pi", "w"))
        else:
            distMat = cPickle.load(open("myData.pi", "r"))

        myDomain = Domain(distMat.keys()[0].get_parent())
        myDomain.addDomain(450, 499)
        myDomain.addDomain(500)
        myDomain2 = Domain(distMat.keys()[0].get_parent())
        myDomain2.addDomain(1, 199)
        myDomain2.addDomain(200, 260)
#         myDomain.mergeDomain(0,1)
        filteredMap = filterInteraction(distMat, myDomain)
        filteredMap2 = filterInteraction(distMat, myDomain2)
        
#         print filteredMap2
        for i in filteredMap2:
            print i.get_parent().get_parent().get_parent()
            print len(list(i.get_parent().get_residues()))
            io = PDB.PDBIO()
            io.set_structure(i.get_parent().get_parent().get_parent())
            return
        
        return

        ## + 3,5s
        hybMat = interaction.intraHydrophobicInteraction(distMat)
        hybMat2 = interaction.intraHydrophobicInteraction(filteredMap)
        hybMat3 = interaction.intraHydrophobicInteraction(filteredMap2)
#         ## + 1,2s
        sbMat = interaction.intraSaltBridge(distMat)
#         ## + 4,1s
        hbMat = interaction.intraHydrogenBond(distMat)

#     vmdrep = rp.VMD_rep(hbMat, "red", "../Structures/R11-15.SAXSmodel.MHP.pdb", "../data/R11_15.tcl")
#     vmdrep.add_map(sbMat, "blue")
#     open("tmp2.txt", "w").write(vmdrep.get_vmdData())
#     open("tmp2.tcl", "w").write(vmdrep.get_vmdScript())

        rrep = rp.R_map(distMat, pdbPath, "")
#         rrep.add_map(sbMat)
        rrep.set_Rdata()
        rrep.set_Rscript()
        data = rrep.get_Rdata()
        script = rrep.get_Rscript()
#         
        open("data.Rdata", "w").write(data)
        open("script.r", "w").write(script)

#         crep = rp.circos_rep(hybMat2, "blue", pdbPath, "../data/R11_15")
#         crep.add_map(hybMat3, "orange")
#         crep.add_map(sbMat, "green")
#         crep.add_map(hbMat, "blue")
#         crep.set_data()
#         crep.set_conf()
#         open("Cdata.data", "w").write(crep.get_data())
#         open("Cdata.link", "w").write(crep.get_link())
#         open("Cdata.seq", "w").write(crep.get_seq())
#         open("Cdata.conf", "w").write(crep.get_conf())


if __name__ == "__main__":
    main()