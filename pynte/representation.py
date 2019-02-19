import sys, random
from Bio.PDB.Polypeptide import three_to_one
import Bio.PDB as PDB
import os

def randomChromColor():
    """
    """
    random.seed()
    return str(#"color=" + 
               str(random.randrange(0, 255)) + "," +
               str(random.randrange(0, 255)) + "," +
               str(random.randrange(0, 255)))


class VMD_rep(object):
    """
    Class producind data and script for all-atom interaction representation
    """
    
    def __init__(self): #, firstMap, firstMapColor, dataPath
        """
        """
        pass
        print "+++++++++++"
        print "Hello VMD !"

        self.Vscript = ""
        self.Vdata = ""

        self.mapID = 1
        self.maps = {} #{1 : (firstMap, firstMapColor)}
        self.set_vmdData()
#         self.pdbPath = pdbPath
        ##############################
#         self.dataPath = dataPath
#         if (self.dataPath[-1] != "/"):
#             self.dataPath + "/"
    
    def __del__(self):
        """
        """
        pass
        print "Bye VMD..."
        print "----------"
        
    def __str__(self):
        """
        """
        return self.get_vmdScript(), self.get_vmdData()
    
    def __len__(self):
        """
        """
        return len(self.maps.keys())
    
    def add_map(self, mapToAdd, mapToAddColor):
        """
        """
        self.maps[self.mapID] = (mapToAdd, mapToAddColor)
        self.mapID = self.mapID * 2
        self.set_vmdData()
    
    def get_vmdScript(self, PDBpath="./pdbName.pdb", tclpath="./tclName.tcl",
                      datapath="./inter.vmd.data"):
        """
        """
        pass
        # Load PDB and setup display
        self.Vscript = '# How to use : vmd -e '+ tclpath +'\n'
        self.Vscript += '# Load PDB \n'
        self.Vscript += 'mol load pdb "'+ PDBpath +'"\n'
        self.Vscript += 'mol delrep 0 0\n'
        self.Vscript += '# Display it in new cartoon\n'
        self.Vscript += 'mol rep newCartoon\n'
        self.Vscript += 'mol color colorID 8\n'
        self.Vscript += 'mol selection [atomselect top "all"]\n'
        self.Vscript += 'mol addrep 0\n'
        self.Vscript += '# Clear Window \n'
        self.Vscript += 'draw delete all\n\n'

        # Load and process data
        self.Vscript += '# Load Data \n'
        self.Vscript += 'set myFile [open "'+ datapath +'" r]\n'
        self.Vscript += '\n'
        self.Vscript += '# Remove last return line char from data file\n'
        self.Vscript += 'set data [string trimright [read $myFile] "\\n"]\n'
        self.Vscript += 'close $myFile\n'
        self.Vscript += '\n'
        self.Vscript += '# List containing ResID of interaction atoms\n'
        self.Vscript += 'set resList [list]\n\n'

        # Display line between interacting residue
        self.Vscript += '# Display dotted line between interaction atoms\n'
        self.Vscript += 'set lines [split $data "\\n"]\n'
        self.Vscript += 'foreach line $lines {\n'
        self.Vscript += '    set fields [split $line "\\t"]\n'
        self.Vscript += '    set sel1 [lindex $fields 0]\n'
        self.Vscript += '    set sel2 [lindex $fields 1]\n'
        self.Vscript += '    set color [lindex $fields 2]\n'
        self.Vscript += '    #set index "index "\n'
        self.Vscript += '    #set resId "resid"\n\n'

        self.Vscript += '    # Select atoms\n'
        self.Vscript += '    set atSel1 [atomselect top $sel1]\n'
        self.Vscript += '    set atSel2 [atomselect top $sel2]\n\n'

        self.Vscript += '    # Select atom position\n'
        self.Vscript += '    lassign [$atSel1 get {x y z}] atPos1\n'
        self.Vscript += '    lassign [$atSel2 get {x y z}] atPos2\n\n'

        self.Vscript += '    # Display line between the two atoms\n'
        self.Vscript += '    draw color $color\n'
        self.Vscript += '    draw line $atPos1 $atPos2 width 1 style dashed\n\n'

        self.Vscript += '    # Select ID of residue from atoms\n'
        self.Vscript += '    # To display them in all atom later\n'
        self.Vscript += '    lassign [$atSel1 get resid] atRes1\n'
        self.Vscript += '    lassign [$atSel2 get resid] atRes2\n'
        self.Vscript += '    lappend resList $atRes1 $atRes2\n'
        self.Vscript += '}\n\n'

        self.Vscript += '# Select unique residue in the resList\n'
        self.Vscript += 'set singRes [lsort -unique $resList]\n'
        self.Vscript += 'set atSel [concat "resid" $singRes]\n\n'

        self.Vscript += 'mol rep CPK\n'
        self.Vscript += 'mol color element\n'
        self.Vscript += 'mol selection $atSel\n'
        self.Vscript += 'mol addrep 0'
        
        return self.Vscript
    
    def set_vmdData(self):
        """
        """
        for ind in self.maps:
            myMap = self.maps[ind][0]
            color = self.maps[ind][1]
             
#             print myMap
#             print "-", ind, color, myMap
            
            for res1 in myMap:
                chain1 = str(res1.get_parent().get_id())
                resname1 = str(res1.get_resname())
                resid1 = str(res1.get_id()[1])

                for res2 in myMap[res1]:
                    chain2 = str(res2.get_parent().get_id())
                    resname2 = str(res2.get_resname())
                    resid2 = str(res2.get_id()[1])

                    for at1 in myMap[res1][res2]:
                        atID1 = str(at1.get_id())
                        
                        for at2 in myMap[res1][res2][at1]:
                            pass
                            self.Vdata += "chain " + chain1 + " and resname " + resname1 \
                                + " and resid " + resid1 + " and name " + atID1 + "\t" \
                                + "chain " + chain2 + " and resname " + resname2 \
                                + " and resid " + resid2 + " and name " + str(at2.get_id()) + "\t" \
                                + color + "\n"
        
        
    def get_vmdData(self):
        """
        """
        return self.Vdata

    def get_vmdPdb(self, PDBpath="./pdbName.pdb"):
        """
        """
        # Get the structure object
        try:
            structure = self.maps[self.maps.keys()[0]][0].keys()[0].get_parent().get_parent().get_parent()
        except IndexError as e:
            pass
#             type, value, traceback = sys.exc_info()
#             print "!! -- coucou"
#             print "!! -- ", e
#             print type, value, traceback
#             quit()
        
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(PDBpath)


class R_map(object):
    """
    """
    def __init__(self, mapR): # , dataPath = "./"
        """
        """
        pass
        print "+++++++++++++"
        print "Hello R map !"
        
        self.Rscript = ""
        self.Rdata = ""
        self.domains = None
        self.additionnalData = ""

#         self.mapID = 1
        self.map = mapR
#         self.set_RData()
#        self.pdbPath = pdbPath
        ##############################
#         self.dataPath = dataPath
#         if (self.dataPath[-1] != "/"):
#             self.dataPath + "/"
        

    def __del__(self):
        """
        """
        pass
        print "bye R map..."
        print "------------"
    
#     def __len__(self):
#         """
#         """
#         return len(self.maps.keys())

    def set_Rdata(self):
        """
        """
        pass
        mapID = self.map
#         for mapID in self.maps:
#             pass
        map, mapLen, start, stop, sequence, position = self.map_Rdata()
#             print start
#             print stop
            
        self.Rdata += "==map==\n"
        self.Rdata += map
        self.Rdata += "==mapLen==\n"
        self.Rdata += mapLen

        minMax = min([i.get_id()[1] for i in start.keys()[0]])
        self.Rdata += "==chainStart==\n"
        # If there is domains
        if self.domains != None:
            for i in self.domains.domainBoundaries.values():
                self.Rdata += str(i[0]-minMax) + "\n"
        else:
            self.Rdata += str("\t".join([str(i) for i in start.values()])) + "\n"

        self.Rdata += "==chainStop==\n"
        # If there is domains
        if self.domains != None:
            for i in self.domains.domainBoundaries.values():
                self.Rdata += str(i[1]-minMax) + "\n"
        else:
            self.Rdata += str("\t".join([str(i) for i in stop.values()])) + "\n"
        self.Rdata += "==colorStart==\n"
        self.Rdata += "blue\n"
        self.Rdata += "==colorStop==\n"
        self.Rdata += "red\n"
        self.Rdata += "==sequence==\n"
            # In lower case, because in Uppercase, R interprete F and T as
            # FALSE and TRUE
        self.Rdata += "\t".join(sequence).lower() + "\n"
#             print "".join(sequence)

        self.Rdata += "==position==\n"
        self.Rdata += position[0]+"\n"
        self.Rdata += position[1]+"\n"
        self.Rdata += position[2]+"\n"
        self.Rdata += position[3]+"\n"
#         self.Rdata += 

        self.Rdata += "==addInfo==\n"
        self.Rdata += "\t".join(self.additionnalData) + "\n"

#         return self.Rdata
        
    
    def clear_Rdata(self):
        """
        """
        pass
        self.Rdata = ""
    
    def map_Rdata(self):
        """
        """
        pass

        allChain = {}       # Key : chain Obj. ; Value : resID of residues
        totalMatLen = 0     # Number of residues represented
        rDataMat = []       # 2D matrix containing 
        rDataMatLen = []    # 2D matrix containing minimum length between residues
        chainStart = {}     # Key : chain Obj. ; Value : first residue position in 2D Matrix
        chainStop = {}      # Key : chain Obj. ; Value : last residue position in 2D Matrix
        resIndInRmap = {}   # Key : residue Obj. ; Value : position in 2D matrix
        position = []       # position in 3D of CA

        # Get all chains represented in the given map
        # In case of different map represent different things
        myMap = self.map

        for res1 in myMap:
            allChain[res1.get_parent()] = []
            for res2 in myMap[res1]:
                allChain[res2.get_parent()] = []

        # Get all residue indices for each chain
        # To access them later
        # In case of tricky numbering
        for i in allChain:
            # Use of the ID because the object sorting lead to mistakes
            allChain[i] = sorted([j.get_id()[1] for j in i])
        
#         # ResID and position
#         print allChain
#         for i in allChain:
#             print i
#             print [j.get_id()[1] for j in i]
#             print [j.get_id()[1] for j in i]
        
        # Size of the final matrix
        for i in allChain:
            totalMatLen += len(allChain[i])

        chainToRind = 0
        # Find where each chain start in
        # And res. position in res. list
        for i in sorted(allChain.keys()):
            pass
            chainStart[i] = chainToRind +1
            for j in allChain[i]:
                resIndInRmap[i[j]] = chainToRind
                chainToRind += 1
            chainStop[i] = chainToRind

        # Final data matrix to be filled
        rDataMat = [[0 for k in xrange(totalMatLen)] for j in xrange(totalMatLen)]
        # Process Map to get R map
        for res1 in myMap:
            resPos1 = resIndInRmap[res1]
            for res2 in myMap[res1]:
                resPos2 = resIndInRmap[res2]

                for at1 in myMap[res1][res2]:
                    for at2 in myMap[res1][res2][at1]:

                        rDataMat[resPos1][resPos2] += 1

        # Final data matrix to be filled
        rDataMatLen = [[0 for k in xrange(totalMatLen)] for j in xrange(totalMatLen)]
        # Process Map to get R map
        for res1 in myMap:
            resPos1 = resIndInRmap[res1]
            for res2 in myMap[res1]:
                resPos2 = resIndInRmap[res2]

                for at1 in myMap[res1][res2]:
                    for at2 in myMap[res1][res2][at1]:

                        if myMap[res1][res2][at1][at2] < rDataMatLen[resPos1][resPos2]:
                            rDataMatLen[resPos1][resPos2] = myMap[res1][res2][at1][at2]
                        if rDataMatLen[resPos1][resPos2] == 0:
                            rDataMatLen[resPos1][resPos2] = myMap[res1][res2][at1][at2]

        # Sequence
        sequence = [""] * totalMatLen
        try:
            for i in resIndInRmap:
                pos = resIndInRmap[i]
                sequence[pos] = three_to_one(i.get_resname())
        except KeyError:
            type, value, traceback = sys.exc_info()
            print("Bad residue name : " + str(value))
            print("Exit prematurely.")
            quit()
             
        # Position
        try:
            ind = []
            x = []
            y = []
            z = []
            for i in sorted(resIndInRmap, key=lambda i: i.get_id()):
                xyz = i["CA"].get_coord()
                x.append(str(xyz[0]))
                y.append(str(xyz[1]))
                z.append(str(xyz[2]))
                ind.append(str(i.get_id()[1]))

            ind = "\t".join(ind)
            x2 = "\t".join(x)
            y2 = "\t".join(y)
            z2 = "\t".join(z)
            position = [x2, y2, z2, ind]

        except KeyError:
            type, value, traceback = sys.exc_info()
            print("Bad residue name : " + str(value))
            print("Exit prematurely.")
            quit()
             
#         mapText = "==Map"+ str(mapID) +"==\n"
        mapText = ""
        for jnd, j in enumerate(rDataMat):
            for knd, k in enumerate(rDataMat[jnd]):
                mapText += str(rDataMat[jnd][knd]) + "\t"
            mapText = mapText[:-1] + "\n"

        mapTextLen = ""
        for jnd, j in enumerate(rDataMatLen):
            for knd, k in enumerate(rDataMatLen[jnd]):
                mapTextLen += str(rDataMatLen[jnd][knd]) + "\t"
            mapTextLen = mapTextLen[:-1] + "\n"

        return mapText, mapTextLen, chainStart, chainStop, sequence, position
    
#     def add_map(self, mapToAdd):
#         """
#         """
#         self.maps[self.mapID] = mapToAdd
#         self.mapID = self.mapID * 2
    
    def get_Rdata(self):
        """
        """
        pass
        return self.Rdata
    
    def set_Rscript(self, dataPath = "./", dataFiles = "data"):
        """
        """
            # Variables
        self.Rscript = "rm(list=ls())\n"
#         self.Rscript += "Sys.setenv(LANG = \"en\")\n"
        self.Rscript += "setwd(\""+ os.path.dirname(os.getcwd() + "/" + dataPath) +"\")\n"
        self.Rscript += 'data = "'+ dataFiles +'"\n'
#         self.Rscript += "imgName = \""+ dataFolder +"\"\n\n"
        # library
        self.Rscript += "suppressMessages(library(\"gplots\"))\n"
        self.Rscript += "suppressMessages(library(\"gdata\"))\n\n"
        # functions entry point
        self.Rscript += "####    Functions    ####\n"
#         self.Rscript += 'addBoundaries <- function(segPos, color){\n'
#         self.Rscript += '  myLineVec = c()\n'
#         self.Rscript += '  for (pos in 1:length(segPos)){\n'
#         self.Rscript += '    myLineVec = c(myLineVec, abline(v=as.integer(segPos[pos]+0.5,\n'
#         self.Rscript += '                                                 lwd=0.1, col="red")))\n'
#         self.Rscript += '    myLineVec = c(myLineVec, abline(h=as.integer(segPos[pos]+0.5,\n'
#         self.Rscript += '                                                 lwd=0.1, col="blue")))\n'
#         self.Rscript += '  }\n'
#         self.Rscript += '  return (myLineVec)\n'
#         self.Rscript += '}\n'
        self.Rscript += "####    End functions    ####\n\n"
    
        # Data file processing
        self.Rscript += '## Parse data from data file ##\n'
        self.Rscript += 'dataFile = readLines(file(data, "r"))\n'
        self.Rscript += 'curField = ""\n'
        self.Rscript += 'listField = c()\n'
        self.Rscript += 'for (line in dataFile){\n'
        self.Rscript += '  if (substr(line, 0, 1) != "#"){\n'
        self.Rscript += '    if (substr(line, 0, 2) == "=="){\n'
        self.Rscript += '      curField = gsub("==", "", line)\n'
        self.Rscript += '      listField = c(listField, curField)\n'
        self.Rscript += '      assign(curField, "")\n'
        self.Rscript += '      eval(parse(text = paste("write(\\\"\\\", file = \\\"",\n'
        self.Rscript += '                        paste("tmp", curField, sep=".", collapse=""),"\\")",\n'
        self.Rscript += '                        sep = "", collapse= "")))\n'
        self.Rscript += '    }\n'
        self.Rscript += '    else if (nchar(line) > 0){\n'
        self.Rscript += '      eval(parse(text = paste("write(\\\"", line, "\\\", file = \\\"",\n'
        self.Rscript += '                              paste("tmp", curField, sep=".", collapse=""),"\\\", append = TRUE)",\n'
        self.Rscript += '                              sep = "", collapse= "")))\n'
        self.Rscript += '} } }\n\n'
         
        self.Rscript += 'for (myField in listField){\n'
        self.Rscript += '  cur_tmp = paste("tmp", myField, sep=".", collapse="")\n'
        self.Rscript += '  eval(parse(text = paste(myField, " = read.table(\\"",\n'
        self.Rscript += '                           cur_tmp, "\\", h=F)",\n'
        self.Rscript += '                           sep = "", collapse= "")))\n'
        self.Rscript += '  eval(parse(text = paste("file.remove(\\"", cur_tmp, "\\")",\n'
        self.Rscript += '                          sep = "", collapse= "")))\n'
        self.Rscript += '}\n#### #### ####\n\n'
    
        # matrixData For each map given
#         for mapID in self.maps:
        self.Rscript += "# Heatmap \n"
        self.Rscript += "posList = c()\n"
        self.Rscript += "colList = c()\n"
        self.Rscript += "for(i in chainStart){\n"
        self.Rscript += "  posList = c(posList, i)\n"
        self.Rscript += "  colList = c(colList, colorStart[1])\n"
        self.Rscript += "}\n"
        self.Rscript += "for(i in chainStop){\n"
        self.Rscript += "  posList = c(posList, i)\n"
        self.Rscript += "  colList = c(colList, colorStop[1])\n"
        self.Rscript += "}\n\n"

        self.Rscript += "if (sum(map[map != -1]) > 0){\n"
#         self.Rscript += "png(\"" + imgFolder + "/" +baseName+ ".png\", width="
#         self.Rscript +=      str(len(sequence)*7)+", height="+str(len(sequence)*7)+")\n"
#             self.Rscript += "sequence[is.na(sequence)] = \"i\" #  ERROR...\n"
        self.Rscript += "sequence = as.matrix(sequence)\n"
#             self.Rscript += "position = as.matrix(position)\n"
        self.Rscript += "myHeatmap = as.matrix(map)\n"
        self.Rscript += "myHeatmap[myHeatmap <= 0.0] = NA\n"
        self.Rscript += "heatmap.2(myHeatmap, Rowv = F, Colv = F, dendrogram=\"none\", trace = \"none\"\n"
        self.Rscript += "          ,col = colorpanel(100, \"blue\", \"green\", \"red\")#, keysize=1.0\n"
        self.Rscript += "          ,revC = T\n"
        self.Rscript += "          ,labRow=toupper(sequence)\n"
        self.Rscript += "          ,labCol=toupper(sequence)#, revC = T\n"
        self.Rscript += "          ,cexRow=0.07, cexCol=0.07\n"
        self.Rscript += "          ,srtCol=45,srtRow=45\n"
#             self.Rscript += "          ,labRow=toupper(paste(position, sequence, sep=\"-\"))\n"
#             self.Rscript += "          ,labCol=toupper(paste(position, sequence, sep=\"-\")), revC = T\n"
        self.Rscript += "          ####    Complete heatmap    ####\n"
#         self.Rscript += "          ,add.expr=addBoundaries(posList, colList)\n"
        self.Rscript += "          ,add.expr=c(abline(v=c(chainStop[,1])+0.5,\n"
        self.Rscript += "                             col=\"red\", lwd=0.5), \n"
        self.Rscript += "                      abline(h=c(chainStop[,1])+0.5,\n" 
        self.Rscript += "                             col=\"red\", lwd=0.5),\n"
        self.Rscript += "                      abline(v=c(chainStart[,1])+0.5,\n"
        self.Rscript += "                             col=\"blue\", lwd=0.5), \n"
        self.Rscript += "                      abline(h=c(chainStart[,1])+0.5,\n" 
        self.Rscript += "                             col=\"blue\", lwd=0.5))\n"
        self.Rscript += "          ####    End heatmap    ####\n"
        self.Rscript += "    )\n"
#         self.Rscript += "dev.off()"
        self.Rscript += "}\n\n"
#             self.Rscript += "Sys.sleep(5)\n\n"
#         break
    
    def get_Rscript(self):
        """
        """
        pass
        return self.Rscript 

    def add_domains(self, domainObject = None):
        """
        """
        pass
    
        if domainObject == None:
            return
        else:
            self.domains = domainObject
        
    def add_info(self, filePath):
        """
        """
        pass
    
        self.additionnalData = open(filePath, "r").read()
        self.additionnalData = self.additionnalData.replace(" ", "_")
    

class circos_rep(object):
    """
    """
    def __init__(self, dataPath): # , firstMap, firstMap_color, pdbPath
        """
        """
        pass
        print "++++++++++++++"
        print "Hello Circos !"
        
        self.Cscript = ""
        self.Cchain = ""
        self.Clink = ""
        self.Cseq = ""
        self.Cdomains = ""
        self.additionnalData = {}
        
        self.chrData = ""
#         self.bandData = ""
        self.Cconf = ""

        
        self.mapID = 1
        self.maps = {}
        self.colors = {}
#         self.pdbPath = pdbPath
        ##############################
#         self.dataPath = dataPath
#         if (self.dataPath[-1] != "/"):
#             self.dataPath + "/"
        

    def __del__(self):
        """
        """
        pass
        print "Bye Circos..."
        print "-------------"

    def __len__(self):
        """
        """
        return len(self.maps.keys())

    def set_conf(self, dataFileName = "./Cdata.data",
                 linkFileName = "./Cdata.link",
                 seqFileName = "./Cdata.seq",
                 imgPathName = "./circos.png"):
        """
        Prepare configuration file for circos.
        """
        pass
        self.Cconf = ""
        self.Cconf += "karyotype = "+dataFileName+"\n"
        self.Cconf += "chromosomes_units           = 1\n"
        self.Cconf += "chromosomes_display_default = yes\n"

        self.Cconf += "\n\n"
        self.Cconf += "<ideogram>\n"
        self.Cconf += "    # Default ideogram settings\n"
        self.Cconf += "    # thickness and color of ideograms\n"
        self.Cconf += "    thickness        = 25p\n"
        self.Cconf += "    stroke_thickness = 1\n"
        self.Cconf += "    stroke_color     = vdgrey\n"
        self.Cconf += "    # the default chromosome color is set here and any value\n"
        self.Cconf += "    # defined in the karyotype file overrides it\n"
        self.Cconf += "    fill             = yes\n"
        self.Cconf += "    fill_color       = black\n"
        self.Cconf += "    # fractional radius position of chromosome ideogram within image\n"
        self.Cconf += "    radius         = 0.90r\n"
        self.Cconf += "    show_label     = no\n"
        self.Cconf += "    label_font     = default\n"
        self.Cconf += "    label_radius   = dims(ideogram,radius) - 0.085r\n"
        self.Cconf += "    label_size     = 6\n"
        self.Cconf += "    label_parallel = yes\n"
        self.Cconf += "    label_case     = upper\n"
        self.Cconf += "    # Bands = AA\n"
        self.Cconf += "    show_bands = yes\n"
        self.Cconf += "    fill_bands = yes\n\n"
        self.Cconf += "    <spacing>\n"
        self.Cconf += "        #default = 0u\n"
        self.Cconf += "        #default = 10u\n"
        self.Cconf += "        #default = 0r\n"
        self.Cconf += "        default = 0.1r\n\n"
        self.Cconf += "    </spacing>\n"
        self.Cconf += "</ideogram>\n\n"
        self.Cconf += "<image>\n"
        self.Cconf += "    dir*    = "+ os.path.dirname(imgPathName) + "\n"
        self.Cconf += "    file*   = "+ os.path.basename(imgPathName) + "\n"
#         self.Cconf += "    radius* = "+str(imgRadius)+"\n"
        self.Cconf += "    <<include etc/image.conf>>\n"
        self.Cconf += "</image>\n\n"
        self.Cconf += "<links>\n"
        self.Cconf += "    <link>\n"
        self.Cconf += "    file    =    "+linkFileName+"\n"
        self.Cconf += "    radius    =    0.95r\n"
        self.Cconf += "    bezier_radius    =    0.35r\n"
        self.Cconf += "    </link>\n"
        self.Cconf += "</links>\n\n"
        
        


        self.Cconf += "<plots>\n"
        self.Cconf += "<plot>\n"
        self.Cconf += "    # Amino acid name\n"
        self.Cconf += "    type    = text\n"
        self.Cconf += "    file    = "+seqFileName+"\n"
        self.Cconf += "    color    = black\n"
        self.Cconf += "    label_snuggle = yes\n"
        self.Cconf += "\n"
        self.Cconf += "    show_links = no\n"
        self.Cconf += "\n"
        self.Cconf += "    r0 = 0.96r\n"
        self.Cconf += "    r1 = 1.00r\n"
        self.Cconf += "    <rules>\n"
        self.Cconf += "    <rule>\n"
        self.Cconf += "    condition  =    1\n"
        self.Cconf += "    value    =    eval(replace(var(value), \"_\", \" \"))\n"
        self.Cconf += "    </rule>\n"
        self.Cconf += "    </rules>\n"
        self.Cconf += "</plot>\n"
        self.Cconf += "</plots>\n\n"
#         self.Cconf += "<<include etc/colors.conf>>\n"
#         self.Cconf += "<<include etc/colors_fonts_patterns.conf>>\n"
#         self.Cconf += "<<include etc/housekeeping.conf>>\n"
    
    def get_conf(self):
        """
        Return configuration file content.
        """
        pass
    
        addConf  = "<<include etc/colors.conf>>\n"
        addConf += "<<include etc/colors_fonts_patterns.conf>>\n"
        addConf += "<<include etc/housekeeping.conf>>\n"
        return self.Cconf + addConf
        
    def get_data(self):
        """
        """
        pass
        return self.Cchain
        
    def set_data(self):
        """
        Prepare the data file content.
        """
        pass
    
        chromLst = {}
        bandLst = {}
        seqList = []
        self.Cchain = ""
        self.Clink = ""
        self.Cseq = ""
        
    
        for mapID in self.maps:
            pass
            start, stop, sequence = self.map_CchainData(mapID)
            self.Clink += self.map_linkData(mapID)
#             print start
#             print stop
            
            for chain in start:
                if chain.get_id() not in chromLst:
                    chromLst[chain.get_id()] = [str(start[chain]), str(stop[chain]), self.colors[mapID], chain]
                    seqList.append([three_to_one(i) for i in sequence])
                else:
                    continue
                    print "Warning, several chain have the same name. It can cause errors if they does not represent the same protein."

        # Chromosomes info == Represented chains
        for indChr, chr in enumerate(sorted(chromLst.keys())):
            myChain = chromLst[chr][3]
            minMax = [i.get_id()[1] for i in myChain]
            
            self.Cchain += "chr - " + chr + " " + str(indChr)
            self.Cchain += " " + str(min(minMax)-1) + " " + str(max(minMax))
            self.Cchain += " " + chromLst[chr][2] + "\n"

        # Band info == residue sequence
        self.Cchain += "\n# band CHRNAME BANDNAME BANDLABEL START END COLOR\n"
        for indChr, chr in enumerate(sorted(chromLst.keys())):
            myChain = chromLst[chr][3]
            ### Additional data
            if self.additionnalData.has_key(chr):
#                 print len(myChain), len(self.additionnalData[chr])
                if len(myChain) != len(self.additionnalData[chr]):
                    self.additionnalData[chr] = " " * len(myChain)
            else:
                self.additionnalData[chr] = " " * len(myChain)
#             print "###" + self.additionnalData[chr] + "###"
            #####
            for ind, i in enumerate(myChain):
                self.Cchain += "\nband" + " " + chr + " " + i.get_resname()
                self.Cchain += "_" + str(i.get_id()[1])
                ##
                if str(self.additionnalData[chr][ind]) != " ":
                    self.Cchain += "_" + str(self.additionnalData[chr][ind])
                    ##
                self.Cchain += " " + i.get_resname()
                self.Cchain += "_" + str(i.get_id()[1])
                ##
                if str(self.additionnalData[chr][ind]) != " ":
                    self.Cchain += "_" + str(self.additionnalData[chr][ind])
                    ##
                self.Cchain += " " + str(i.get_id()[1] -1)
#                 self.Cchain += " " + str(i.get_id()[1]) + " " + randomChromColor()
#                 self.Cchain += " " + str(i.get_id()[1]) + " " + chromLst[chr][2]
                self.Cchain += " " + str(i.get_id()[1]) + " 200,25,126"

                self.Cseq += chr + " " + str(i.get_id()[1] -1)
                self.Cseq += " " + str(i.get_id()[1])
                self.Cseq += " " + three_to_one(i.get_resname())
                self.Cseq += "_" + str(i.get_id()[1])
                ##
                if str(self.additionnalData[chr][ind]) != " ":
                    self.Cseq += "_" + str(self.additionnalData[chr][ind])
                    ##
                self.Cseq += " label_size=8p"
                self.Cseq += "\n"

#         return self.Rdata
    
    def map_CchainData(self, mapID):
        """
        """
        pass

        allChain = {}       # Key : chain Obj. ; Value : resID of residues
        totalMatLen = 0     # Number of residues represented
        chainStart = {}     # Key : chain Obj. ; Value : first residue position in 2D Matrix
        chainStop = {}      # Key : chain Obj. ; Value : last residue position in 2D Matrix
        resIndInRmap = {}   # Key : residue Obj. ; Value : position in 2D matrix

        # Get all chains represented in the given map
        # In case of different map represent different things
        myMap = self.maps[mapID]

        for res1 in myMap:
            allChain[res1.get_parent()] = []
            for res2 in myMap[res1]:
                allChain[res2.get_parent()] = []

        # Get all residue indices for each chain
        # To access them later
        # In case of tricky numbering
        for i in allChain:
            # Use of the ID because the object sorting lead to mistakes
            allChain[i] = sorted([j.get_id()[1] for j in i])
        # Size of the final matrix
        for i in allChain:
            totalMatLen += len(allChain[i])

        chainToRind = 0
        # Find where each chain start in
        # And res. position in res. list
        for i in sorted(allChain.keys()):
            pass
            chainStart[i] = chainToRind
            for j in allChain[i]:
                resIndInRmap[i[j]] = chainToRind
                chainToRind += 1
            chainStop[i] = chainToRind

        # Final data matrix to be filled
        rDataMat = [[0] * totalMatLen for j in xrange(totalMatLen)]
        # Process Map to get R map
        for res1 in myMap:
            resPos1 = resIndInRmap[res1]
            for res2 in myMap[res1]:
                resPos2 = resIndInRmap[res2]

                for at1 in myMap[res1][res2]:
                    for at2 in myMap[res1][res2][at1]:

                        rDataMat[resPos1][resPos2] += 1

        # Sequence
        sequence = [""] * totalMatLen
        for i in resIndInRmap:
            pos = resIndInRmap[i]
            sequence[pos] = i.get_resname()
             
#         mapText = "==Map"+ str(mapID) +"==\n"
        mapText = ""
        for jnd, j in enumerate(rDataMat):
            for knd, k in enumerate(rDataMat[jnd]):
                mapText += str(rDataMat[jnd][knd]) + "\t"
            mapText = mapText[:-1] + "\n"

        return chainStart, chainStop, sequence#, mapText
    
    def map_linkData(self, mapID):
        """
        """
        pass

        allChain = {}       # Key : chain Obj. ; Value : resID of residues
        totalMatLen = 0     # Number of residues represented
        rDataMat = []       # 2D matrix containing 
        chainStart = {}     # Key : chain Obj. ; Value : first residue position in 2D Matrix
        chainStop = {}      # Key : chain Obj. ; Value : last residue position in 2D Matrix
        resIndInRmap = {}   # Key : residue Obj. ; Value : position in 2D matrix

        # Get all chains represented in the given map
        # In case of different map represent different things
        myMap = self.maps[mapID]
        myColor = self.colors[mapID]
#         print myColor
        link = ""
        
        for res1 in myMap:
            for res2 in myMap[res1]:
                count = 0

                for at1 in myMap[res1][res2]:
                    for at2 in myMap[res1][res2][at1]:
                        count += 1


                link += str(res1.get_parent().get_id()) + " "
                link += str(res1.get_id()[1]-1) + " "
                link += str(res1.get_id()[1]-1 + 1) + " "
                link += str(res2.get_parent().get_id()) + " "
                link += str(res2.get_id()[1]-1) + " "
                link += str(res2.get_id()[1]-1 + 1) + " "
                link += "color=" + myColor
                link += "\n"
                    
        return link
    
    def get_link(self):
        """
        """
        pass
    
        return self.Clink
    
    def get_seq(self):
        """
        """
        pass
    
        return self.Cseq
    
    def add_map(self, mapToAdd, color="blue"):
        """
        """
        
        self.maps[self.mapID] = mapToAdd
        self.colors[self.mapID] = color
        self.mapID = self.mapID * 2

    def add_domains(self, domainObject = None):
        """
        """
        pass
    
        if domainObject == None:
            return
#         print domainObject.domainBoundaries
        
        tmpStr = ""
#         print self.Cchain
        for i in self.Cchain.split("\n"):
            # Empty lines or 
            # lines that does not start by "band"
            if not (len(i) > 0 and  i[0] == "b"):
                tmpStr += i
            else:
                j = i.split(' ')
                # If a residue is in a domain
                if int(j[5]) in domainObject.domainBoundariesInverted:
                    dom = domainObject.domainBoundariesInverted[int(j[5])]
                    j[6] = str(str((dom*200+100) % 255) + "," + 
                               str((dom*100+30) % 255) + "," + 
                               str((dom*50+100) % 255))
                else:
                    j[6] = '0,0,0'
                tmpStr += " ".join(j)
            tmpStr += "\n"

        self.Cchain = tmpStr
        
    def add_info(self, filePath, chain):
        """
        """
        pass
    
        self.additionnalData[chain] = open(filePath, "r").read()
        self.additionnalData[chain] = self.additionnalData[chain].replace(" ", "_")