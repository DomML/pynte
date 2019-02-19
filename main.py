#!/usr/local/bin/python
################################################################################
# * PyNte-CLI (command line interface)
# Dominique MIAS-LUCQUIN @ IGDR
# dominique.mias-lucquin@univ-rennes1.fr
# 
# A software to study interaction in protein structures
# 
# Licence of use : ???
# 
################################################################################
import argparse
from gooey import Gooey

import pynte.distanceMatrix as pdm
from pynte import representation as rp
from pynte import interaction
from pynte.selection import *
from Bio import PDB
import cPickle
from urlparse import ParseResult


# @Gooey
def main():
    """
    """
    parser = argparse.ArgumentParser(description='PyNte-CLI (command line interface)')
    parser.add_argument('-p', '--pdb', action="store", type=str)
    parser.add_argument('-y', '--hydrophobicity', action="store", type=float, default=-1.0, help="If a paltinum's PDB is given with '-p', compute hydrophobic interactions with the given cutoff.")
    parser.add_argument('-b', '--hydrogen-bond', action="store", default=-1.0, type=float)
    parser.add_argument('-s', '--salt-bridge', action="store", default=-1.0, type=float)
    parser.add_argument('-c', '--cut-off', action="store", default=50.0, type=float, nargs=1, help="General cut-off distance for distance matrix computation.")
    parser.add_argument('-l', '--load-map', action="store", help="Load the distance matrix to avoid to compute it again and again and again.", type=file)
    parser.add_argument('-m', '--dist-map', action="store", help="Save the distance matrix to avoid to compute it again and again and again.", type=str, default=None)
    parser.add_argument('-r', '--res-filter', action="store", help="Limit interactions search to residues given into parameter. Use 1-letter code (-r ACTNG, to limit to Alanin, Cystein, ..., Glycin)", type=str, default=None)
    parser.add_argument('-o', '--output', action="store", help="Output directory.", type=str, default="./")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s alpha_1')
    parser.add_argument('-d', '--domain', action='append', help="Interaction domains, format : start:end", type=str)
    parser.add_argument('-a', '--add-info', action='store', help="File adding a layer of text information on circos (1 letter by residue)", type=str)

    parseResults = parser.parse_args()
    print parseResults
    outputFolder = parseResults.output
    if (outputFolder[-1] != "/"):
        outputFolder + "/"
#     print parseResults
#     return

    #####################################
    #++++++++++++++++++++++++++++++++++++
    # Load distance matrix to avoid to compute it again and again and again.
    if parseResults.load_map != None:
        distanceMap = cPickle.load(open("myData.pi", "r"))
    else:
        # Else, load PDB, etc.
        
        #####################################
        # General cutoff distance
        cutoff = parseResults.cut_off
        hydrophobicityDist = parseResults.hydrophobicity
        hydrogen_bondDist = parseResults.hydrogen_bond
        salt_bridgeDist = parseResults.salt_bridge
        
        zplus = 5
        cutoff = max([cutoff, hydrogen_bondDist + zplus,
                      hydrophobicityDist + zplus,
                      salt_bridgeDist + zplus])
        #####################################
    
        # Open and parse the PDB
        pdb = parseResults.pdb
        try:
            pdbID = parseResults.pdb.split('/')[-1].split('.')[0]
        except AttributeError:
            return
            pass
    
        try:
            myPDBparser = PDB.PDBParser(PERMISSIVE=1, QUIET=True)
            myStructure = myPDBparser.get_structure(pdbID, pdb)
        except IOError:
            print("Bad PDB file path.")
            print("Exit prematurely.")
            return
        
#         # Renumber residue from 1
#         for i in myStructure:
#             for j in i:
#                 mini = min([k.get_id()[1] for k in j]) - 1
#                 for k in j:
#                     k.id = (k.get_id()[0], k.get_id()[1] - mini, k.get_id()[2])
#         
#         print myStructure
        
        ######################################
        # Intra or inter-molecular interactions
        if len(list(myStructure.get_chains())) == 1:
            # Intra-molecular
            singleChain = list(myStructure[0].get_chains())[0]

            if singleChain.get_id() == " ":
                singleChain.id = "Z"

            distanceMap = pdm.intraMolecularDistanceMatrix(singleChain, cutoff, 20)
            print cutoff
            
            
            singleChainID = singleChain.get_id()
        else:
            print "Inter-molecular interactions are not currently supported."
            return
            pass
    #------------------------------------

    #####################################
    # Residue filtering
    if parseResults.res_filter!= None:
        pass
        distanceMap = residue_filter(distanceMap, parseResults.res_filter)
    
    # Domains filtering
    domainsList = None
    if (parseResults.domain != None):
        domainsList = Domain(distanceMap.keys()[0].get_parent())
        for i in parseResults.domain:
            i = map(int, i.split(":")[0:2])
            domainsList.addDomain(i[0], i[1])
        
        distanceMap = domainsList.filterInteraction(distanceMap)
    #------------------------------------

    #####################################
    # Interactions and Representation script
    vmdRep = rp.VMD_rep()
    circosRep = rp.circos_rep(parseResults.output)
    #++++++++++++++++++++++++++++++++++++
    # Hydrophobic interaction (hydrophobicityMap)
    if (hydrophobicityDist > 0.0):
        interMap = interaction.intraHydrophobicInteraction(distanceMap,
                                                           hydrophobicityDist)
             

        ####
        # R
        rRep = rp.R_map(interMap)
        rRep.add_domains(domainsList)
        # Add file info (singlaChainID)
        if parseResults.add_info != None:
            rRep.add_info(parseResults.add_info)

        rRep.set_Rdata()
        rRep.set_Rscript(outputFolder, outputFolder + "hydrophoby_Rdata.data")
        open(outputFolder + "hydrophoby_Rdata.data", "w").write(rRep.get_Rdata())
        open(outputFolder + "hydrophoby_Rscript.r", "w").write(rRep.get_Rscript())
        del rRep
        ######
        # VMD
        vmdRep.add_map(interMap, "blue")
        #########
        # Circos
        circosRep.add_map(interMap, "blue")
    #++++++++++++++++++++++++++++++++++++
    # Hydrogen bond (hydrogenBondMap)
    if (hydrogen_bondDist > 0.0):
        interMap = interaction.intraHydrogenBond(distanceMap, hydrogen_bondDist)
        ####
        # R
        rRep = rp.R_map(interMap)
        rRep.add_domains(domainsList)
        # Add file info (singlaChainID)
        if parseResults.add_info != None:
            rRep.add_info(parseResults.add_info)

        rRep.set_Rdata()
        rRep.set_Rscript(outputFolder, outputFolder + "hydrogenBond_Rdata.data")
        open(outputFolder + "hydrogenBond_Rdata.data", "w").write(rRep.get_Rdata())
        open(outputFolder + "hydrogenBond_Rscript.r", "w").write(rRep.get_Rscript())
        del rRep
        ######
        # VMD
        vmdRep.add_map(interMap, "red")
        #########
        # Circos
        circosRep.add_map(interMap, "red")
    #++++++++++++++++++++++++++++++++++++
    # Salt bridge (saltBridgeMap)
    if (salt_bridgeDist > 0.0):
        interMap = interaction.intraSaltBridge(distanceMap, salt_bridgeDist)
        ####
        # R
        rRep = rp.R_map(interMap)
        rRep.add_domains(domainsList)
        # Add file info (singlaChainID)
        if parseResults.add_info != None:
            rRep.add_info(parseResults.add_info,)

        rRep.set_Rdata()
        rRep.set_Rscript(outputFolder, outputFolder + "saltBridge_Rdata.data")
        open(outputFolder + "saltBridge_Rdata.data", "w").write(rRep.get_Rdata())
        open(outputFolder + "saltBridge_Rscript.r", "w").write(rRep.get_Rscript())
        del rRep
        ######
        # VMD
        vmdRep.add_map(interMap, "green")
        #########
        # Circos
        circosRep.add_map(interMap, "green")
    #------------------------------------
    
    ########
    # Circos : Add file info (singlaChainID)
    if parseResults.add_info != None:
        circosRep.add_info(parseResults.add_info, singleChainID)
    
    #++++++++++++++++++++++++++++++++++++
    if len(vmdRep) > 0:
        pass
        open(outputFolder + "interactions.vmd.data", "w").\
            write(vmdRep.get_vmdData())
        open(outputFolder + "interactions.tcl", "w").\
            write(vmdRep.get_vmdScript(pdb, 
                                       outputFolder + "interactions.tcl",
                                       outputFolder + "interactions.vmd.data"))
        vmdRep.get_vmdPdb("struct.pdb")
    if len(circosRep) > 0:
        pass
        circosRep.set_data()
        circosRep.add_domains(domainsList)
        circosRep.set_conf(outputFolder + "Cdata.data", 
                           outputFolder + "Cdata.link",
                           outputFolder + "Cdata.seq",
                           outputFolder + ".png")
    
        open(outputFolder + "Cdata.data", "w").write(circosRep.get_data())
        open(outputFolder + "Cdata.link", "w").write(circosRep.get_link())
        open(outputFolder + "Cdata.seq", "w").write(circosRep.get_seq())
        open(outputFolder + "Cdata.conf", "w").write(circosRep.get_conf())
    #------------------------------------

    #####################################
    #++++++++++++++++++++++++++++++++++++
    # Save distance matrix
    if parseResults.dist_map != None:
        pass
        cPickle.dump(distanceMap, open(parseResults.dist_map, "w"))
    #------------------------------------

 
if (__name__ =="__main__"):
    main()