# !/usr/bin/env python
######################################################################################################################################
# Automatic Zmatrix script generator for BOSS and MCPRO 
#
# Authors: Israel Cabeza de Vaca Lopez and Matthew Robinson
#
# Script based on the README notes written by Dr. Julian Tirado-Rives 
#
# Usage: python pdb2zmat.py -p pdb_file -r residue_name -c cutoff_size
#
#           use --help for further details or instructions
#
# Outputs:
#       It generates a folder called finalZmatrices with the final zmatrix files (all, cap, cap + conrot for flexible backbones) 
# 
# Requirements:
#       BOSS
#       MCPRO
#       Reduce executable (http://kinemage.biochem.duke.edu/software/reduce.php)
#       Propka-3.1 executable (https://github.com/jensengroup/propka-3.1)
#       Babel
#       Python 3.6
#       
######################################################################################################################################

import os
import sys
import subprocess
from biopandas.pdb import PandasPdb
import pandas as pd
import argparse

def main():

    # create parser object
    parser = argparse.ArgumentParser(
        prog='pdb2zmatConrotGenerator_mcr.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
    Automatic Zmatrix script generator for BOSS and MCPRO 

    @author: Israel Cabeza de Vaca Lopez, israel.cabezadevaca@yale.edu   
    @author: Matthew Robinson, matthew.robinson@yale.edu
    @author: Yue Qian, yue.qian@yale.edu
    @author: William L. Jorgensen Lab 

    Example usage: python pdb2zmat.py -p 4WR8_ABC_3TX.pdb -r 3TX -c 18.0

    Or, if you already have an optimized z-matrix:

    Usage: python pdb2zmat.py -p 4WR8_ABC_3TX.pdb -z 3TX_h.z -r 3TX -c 18.0
    
    REQUIREMENTS:
    BOSS (need to set BOSSdir in bashrc and cshrc)
    MCPRO (need to set MCPROdir in bashrc and cshrc)
    reduce executable (from Richardson lab)
    propka-3.1 executable (from Jensen lab)
    Preferably Anaconda python 3.6 with following modules:
    pandas 
    biopandas
    Babel
    """
    )

    #defining arguments for parser object
    parser.add_argument(
        "-p", "--pdb", help="name of the PDB file for the complex - i.e., [complex_name].pdb")
    parser.add_argument(
        "-z", "--zmat", help="name of the zmat of the ligand with .z file descriptor. \
        only need this if want to import the optimized ligand zmat yourself.")
    parser.add_argument(
        "-r", "--resname", help="Residue name of the ligand from PDB FILE", type=str)
    parser.add_argument(
        "-c", "--cutoff", help="Size of the cutoff for chop to cut", type=str)
    parser.add_argument(
        "-v", "--variable_cutoff", help="Size of the variable cutoff for chop to cut", type=str)
    parser.add_argument(
        "-t", "--tautomerize", help="Automatically assign Histidine Tautomerization states", action="store_true")

    #parse the arguments from standard input
    args = parser.parse_args()

    fixedComplexPDB, ligandZmat = prepare_complex(args.zmat, args.pdb, args.resname, args.cutoff, args.variable_cutoff, args.tautomerize)

    # note that prepare_zmats takes in the FIXED PDB
    prepare_zmats(ligandZmat, fixedComplexPDB, args.resname)

def prepare_complex(zmat_arg, pdb_arg, resname_arg, cutoff_arg, variable_cutoff_arg, his_arg):

    # preliminary definitions 
    # please only change these if you know what you are doing!

    complexPDB = '' #'1wbv.pdb'
    ligandZmatOrig = '' #'benzene.z'
    ligandResidueName = '' #'LI3'

    ligandLstToRemoveFromPDB = [] #recommend you start this empty, so I need to get rid of these to start
    residueToBeTheCenterOfTheChopping = '' #'LI3'   # Normally the ligand
    setCapOriginAtom = '' #'LI3'  # LIGAND name ex. 'LI3'   
    setCutOffSize = '' #'18.0'

    titleOfYourSystem = '' # optional

    fixBackBoneSelection = [] #['4 67 70 74 110 144 152 153'] # If you have a lot of residues split the selection in different lines

    if zmat_arg:
        ligandZmatOrig = str(zmat_arg)

    if pdb_arg:
        complexPDB = str(pdb_arg)

    if resname_arg:
        ligandResidueName = str(resname_arg)
        residueToBeTheCenterOfTheChopping = str(resname_arg) 
        setCapOriginAtom = 'c01'#str(args.resname)

    if cutoff_arg:
        setCutOffSize = str(cutoff_arg)

    if variable_cutoff_arg:
        setVariableCutOffSize = str(variable_cutoff_arg)
    else:
        setVariableCutOffSize = '10'

    # *********************** CODE STARTS *********************************************

    checkParameters(complexPDB)

    checkMultipleLigands(complexPDB,ligandResidueName)

    changeLigandResidueNumber(complexPDB,ligandResidueName)

    if zmat_arg:
        _, resnumber, resnumberLigandOriginal = generateLigandPDB(complexPDB,ligandResidueName)
        ligandZmat = fixDummyAtomNamesDirect(ligandZmatOrig)

    else:
        ligandZmat,resnumber,resnumberLigandOriginal = prepareLigandZmatrix(complexPDB,ligandResidueName,mcproPath,BOSSscriptsPath)


    #remove solvent
    no_solvent_PDB = removeSolvent(complexPDB)

    #change Chain ID of ligand
    fixedComplexPDB = fixPDBprotein(complexPDB,ligandResidueName)


    print('MERGE')
    mergeLigandWithPDBProtein(mcproPath,fixedComplexPDB,ligandZmat,resnumberLigandOriginal,ligandResidueName)

    #get histdine lists before chopping
    #can also comment this section out if desired


    HipLstOfResidues = []  # resnumber, Chain ex. ['77a','56b'] #optional
    #HieLstOfResidues = []  # resnumber, Chain this the default!!! original code wrong
    HidLstOfResidues = []
    if his_arg:

        HipLstOfResidues, HidLstOfResidues = makeHisLists(pdb_arg,fixedComplexPDB,HipLstOfResidues,HidLstOfResidues)
        print(HipLstOfResidues)
        print(HidLstOfResidues)

    print('CHOP')
    prepareReducedChopped(fixedComplexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,setVariableCutOffSize, HipLstOfResidues,HidLstOfResidues)

    return fixedComplexPDB, ligandZmat

def changeLigandResidueNumber(complexPDB, ligandName):
    
    print('CHANGING LIGAND RESIDUE NUMBER') #needs to be 100 to work with Clu

    # Read in the file
    with open(complexPDB, 'r') as pdb_file:
        pdb_data = pdb_file.readlines()

    with open('original_'+complexPDB, 'w') as file:
        file.writelines(pdb_data)

    # Replace the res number with 100 for clu to work
    # Replace the res number with 100 for clu to work
    idx = 0
    for line in pdb_data:
        if 'ATOM' in line[:7] or 'HETATM' in line[:7]:
            if ligandName in line:

                newLine = line[:21]+' '+line[22:] #get rid of chain ID
                resnumber = newLine[12:29].split()[2]
                ligandResnumber = resnumber
                if int(resnumber)<=999:
                    ligandResnumber = '100'
                    newLine = newLine[:23] + '100' + newLine[26:]
                    #pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
                elif int(resnumber)>999:
                    ligandResnumber = ' 100'
                    newLine = newLine[:22] + ' 100' + newLine[26:]
                    # pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
            else:
                newLine = line
        else:
            newLine = line

        pdb_data[idx] = newLine
        idx = idx + 1

    # Write the file out again
    with open(complexPDB, 'w') as file:
        file.writelines(pdb_data)


def runPropka(originalPDB,HipLstOfResidues):
    #run propka on the original protein

    try:
        os.system('propka31 ' + originalPDB)
    except:
        print('propka failed on the command:', ('propka31 ' + originalPDB))
        sys.exit()

    propka_output_name = originalPDB[:-4]+'.pka'

    #read propka output
    try:
        with open(propka_output_name) as propka_output:
            pka_data = propka_output.readlines()
    except:
        print('PROPKA FAILED ON THIS PROTEIN')
        return

    #get only the pKa data
    for i in range(len(pka_data)):
        if (pka_data[i][0:7] == 'SUMMARY'):
            pka_data = pka_data[i+1:]
            break

    #read the histidine residues
    propka_hip_list = [] 
    for line in pka_data:
        line = line.rstrip()
        if line.split() == []:
            continue
        res_name = line.split()[0]
        if res_name == 'HIS':
            res_number = line.split()[1] + line.split()[2]
            pka = line.split()[3]
            if (float(pka)>=7.0):
                propka_hip_list.append(res_number)

    #compare the two lists
    for res in propka_hip_list:
        if res not in HipLstOfResidues:
            print("prokpka thinks " + res + " should be HIP, reduce does not")

    for res in HipLstOfResidues:
        if res not in propka_hip_list:
            print("reduce thinks " + res + " should be HIP, propka does not")

def makeHisLists(originalPDB,complexPDB,HipLstOfResidues,HidLstOfResidues):

    pdb_name = complexPDB[:-4]+'_cplx.pdb'

    pdb_out = complexPDB[:-4]+'_cplx_h.pdb'

    # run reduce on the protein
    #subprocess.call("reduce –build " + pdb_name + " > " + pdb_name[:-4]+ "_h.pdb", shell=True)
    try:
        os.system("reduce -Build %s > %s" % (pdb_name, pdb_out))
    except:
        print("reduce failed on the command:", ("reduce -Build %s > %s" % (pdb_name, pdb_out)))
        
    
    #read in the protein with biopandas
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_out)
    atom_df = ppdb.df['ATOM']

    #construct variables for boolean selection
    HID = atom_df['atom_name']=='HD1'
    HIE = atom_df['atom_name']=='HE2'
    HIS = atom_df['residue_name']=='HIS'
    
    #make dataframes
    hid_df = atom_df.loc[HIS & HID]
    hie_df = atom_df.loc[HIS & HIE]
    
    #make booleans to check it has that type of his
    has_hid = hid_df.shape[0] > 0
    has_hie = hie_df.shape[0] > 0
    
    #construct the lists
    hid_list = []
    hie_list = []
    hip_list = []
    
    #first for his
    if has_hid:
        hid_nums = [str(x) for x in list(hid_df['residue_number'])]
        hid_chains = list(hid_df['chain_id'])
        hid_list = [hid_nums[i] + hid_chains[i] for i in range(len(hid_nums))]
    #now for hie
    if has_hie:
        hie_nums = [str(x) for x in list(hie_df['residue_number'])]
        hie_chains = list(hie_df['chain_id'])
        hie_list = [hie_nums[i] + hie_chains[i] for i in range(len(hie_nums))]
    
    #construct the hip list
    if has_hid and has_hie:
        for res in hie_list:
            if res in hid_list:
                hip_list.append(res)

    print('hiplist')
    print(hip_list)
    print('hielist')
    print(hie_list)
    print('hidlist')
    print(hid_list)

    #add to the global lists
    for res in hip_list:
        if (res not in HipLstOfResidues) and (res not in HidLstOfResidues):
            HipLstOfResidues.append(res)

    for res in hid_list:
        if (res not in HipLstOfResidues) and (res not in HidLstOfResidues):
            HidLstOfResidues.append(res)

    #check if propka disagrees with reduce
    runPropka(originalPDB,HipLstOfResidues)

    return HipLstOfResidues, HidLstOfResidues     

def checkParameters(complexPDB):

    if not mcproPath:
        print('Define MCPROdir enviroment variable, please')
        sys.exit()
    else:
        print('Using MCPROdir .....    ',mcproPath)

    if not BOSSPath:
        print('Define BOSSdir enviroment variable, please')
        sys.exit()
    else:
        print('Using BOSSdir .....    ',BOSSPath)

    if not os.path.isfile(complexPDB):
        print('PDB not found ......    ',complexPDB)
        sys.exit()

def checkMultipleLigands(complexPDB,ligandResidueName):

    ppdb = PandasPdb()
    ppdb.read_pdb(complexPDB) 
    atom_df = ppdb.df['ATOM']
    hetatm_df = ppdb.df['HETATM']

    # get only information about the ligand
    lig_df = hetatm_df.loc[hetatm_df['residue_name'] == ligandResidueName]

    # check if lig_df contains multiple chain ids(means we have multiple ligands)
    if lig_df.chain_id.nunique() > 1:
        print("THERE IS MORE THAN COPY OF THE DESIRED LIGAND. PLEASE DELETE EXTRAS")
        sys.exit()


def generateLigandPDB(complex,ligandName):

    print('Generating Ligand PDB')

    namePDBLigand = ligandName.lower()+'.pdb'
    pdbLigandFile = open(namePDBLigand,'w')

    ligandResnumber = ''

    ligandResnumberOriginal = ''

    for line in open(complex):
        if 'ATOM' in line[:7] or 'HETATM' in line[:7]:
            if ligandName in line[15:23]:
                
                newLine = line[:21]+' '+line[22:]
                resnumber = newLine[12:29].split()[2]
                ligandResnumber = resnumber
                ligandResnumberOriginal = resnumber
                if int(resnumber)<=999:
                    ligandResnumber = '100'
                    newLine = newLine[:23] + '100' + newLine[26:]
                    #pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
                    pdbLigandFile.write(newLine)
                elif int(resnumber)>999:
                    ligandResnumber = ' 100'
                    newLine = newLine[:22] + '0100' + newLine[26:]
                    # pdbLigandFile.write(newLine.replace(resnumber,ligandResnumber))
                    pdbLigandFile.write(newLine)
                else:
                    pdbLigandFile.write(newLine)

    pdbLigandFile.close()

    #complexPDBfixName = fixComplexResidueNumber(complex,ligandName)

    return namePDBLigand,ligandResnumber,ligandResnumberOriginal#,complexPDBfixName
    

def protonateLigandWithBabel(namePDBLigand):

    print('Protonating ligand with Babel')

    namePDBLigand_protonated = namePDBLigand[:-4]+'_h.pdb'

    # chimeraScriptFile = open('protonate.cmd','w')
    # chimeraScriptFile.write(chimeraScript.replace('aaaa',namePDBLigand).replace('bbbb',namePDBLigand_protonated))
    # chimeraScriptFile.close()

    #cmd = '/Applications/Chimera.app/Contents/Resources/bin/chimera --nogui protonate.cmd'
    # cmd = 'chimera --nogui protonate.cmd'
    cmd = ('babel %s -O %s -p' % (namePDBLigand,namePDBLigand_protonated))

    try:
        os.system(cmd)
    except:
        print("babel failed to perform the command:", cmd)
        sys.exit()

    return namePDBLigand_protonated 

def fixDummyAtomNames(namePDBLigand_protonated):
    
    print('Fixing DUMMY atoms names')

    filetmp = open('tmp.txt','w')

    ligandZmat = namePDBLigand_protonated[:-4]+'.z'

    for line in open(ligandZmat):
        newLine = line
                
        if len(line.split())==12 and len(line)==77:
            newLine = line[:71]+'    1\n'
        filetmp.write(newLine.replace('   1 DUM','   1 DU1').replace('   2 DUM','   2 DU2').replace('   3 O00  800  800','   3 DU3   -1    0'))

    filetmp.close()

    cmd = 'cp tmp.txt '+ligandZmat
    
    print(cmd)

    os.system(cmd)

    return ligandZmat

def fixDummyAtomNamesDirect(ligandZmatOrig):

    print('Fixing DUMMY atoms names')

    filetmp = open('tmp.txt', 'w')

    ligandZmat = ligandZmatOrig[:-2] + '_fixed.z'

    for line in open(ligandZmatOrig):
        newLine = line

        if len(line.split()) == 12 and len(line) == 77:
            newLine = line[:71] + '    1\n'
        filetmp.write(
            newLine.replace('   1 DUM', '   1 DU1').replace('   2 DUM', '   2 DU2').replace('   3 O00  800  800',
                                                                                            '   3 DU3   -1    0'))

    filetmp.close()

    cmd = 'cp tmp.txt ' + ligandZmat

    print(cmd)

    os.system(cmd)

    return ligandZmat

def optimizeStructure(ligandZmat,mcproPath):

    print('Optimizing Structure')

    cmd = mcproPath+'/scripts/xOPT '+ligandZmat[:-2]

    print(cmd)

    os.system(cmd)

def prepareLigandZmatrix(complex,ligandName,mcproPath,BOSSscriptsPath):

    print('Preparing Ligand Z matrix')

    namePDBLigand, ligandResnumber, ligandResnumberOriginal = generateLigandPDB(complex,ligandName)

    namePDBLigand_protonated = protonateLigandWithBabel(namePDBLigand)

    cmd = BOSSscriptsPath+'/xPDBMCP '+namePDBLigand_protonated[:-4]

    os.system(cmd)

    ligandZmat = fixDummyAtomNames(namePDBLigand_protonated)
    
    optimizeStructure(ligandZmat, mcproPath)

    return ligandZmat,ligandResnumber,ligandResnumberOriginal#,complexPDBfixName

def mergeLigandWithPDBProtein(mcproPath,complexPDB,ligandZmat,resnumber,resname):
    #clu -t:s=5001 2be2.pdb -r r22_h.z -n 2be2_cplx.pdb

    print('Merging ligand with PDB Protein')
    
    clu = mcproPath+'/miscexec/clu '

    outputPDB = complexPDB[:-4]+'_cplx.pdb'
    
    if os.path.isfile(outputPDB): os.remove(outputPDB)

    # THIS RESNUMBER MUST BE THE SAME AS BEOFRE 
    # cmd = clu + ' -t:s='+resnumber+' '+complexPDB+' -r '+ligandZmat+' -n '+outputPDB
    cmd = cmd = clu + ' -t:s=:'+resname+' '+complexPDB+' -r '+ligandZmat+' -n '+outputPDB

    print(cmd)

    os.system(cmd)

def generateScript(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,setVariableCutOffSize, HipLstOfResidues,HidLstOfResidues):

    lastPart = ''' write pdb aaaa.chop.pdb
write pepz all aaaa.chop.all.in
write pepz variable aaaa.chop.var.in
write translation aaaa.chop.tt
EOF

'''

    complexPDBName = complexPDB[:-4]

    tmpScript = '$MCPROdir/miscexec/chop -t -i '+complexPDBName+'_cplx.pdb << EOF\n'

    if len(ligandLstToRemoveFromPDB)!=0:
        for ele in ligandLstToRemoveFromPDB:
            tmpScript = tmpScript + 'delete ligand :'+ele+'\n'

    if len(HipLstOfResidues)!=0:
        lst = ''
        for ele in HipLstOfResidues:
            lst = lst + ':'+ele.strip()+' '
        tmpScript = tmpScript + 'set hip '+lst+'\n\n'


    if len(HidLstOfResidues)!=0:
        lst = ''
        for ele in HidLstOfResidues:
            lst = lst + ':'+ele.strip()+' '
        tmpScript = tmpScript + 'set hid '+lst+'\n\n'

    tmpScript = tmpScript + 'add center :'+residueToBeTheCenterOfTheChopping+'\n'
    tmpScript = tmpScript + 'set cap origin :'+setCapOriginAtom+'\n'

    tmpScript = tmpScript + 'set cut origin ligand \n'
    tmpScript = tmpScript + 'set cut size '+setCutOffSize+'\n'
    #tmpScript += 'delete cut :gol\n'
    if len(ligandLstToRemoveFromPDB)!=0:
        for ele in ligandLstToRemoveFromPDB:
            tmpScript = tmpScript + 'delete cut :'+ele.lower()+'\n'

    tmpScript = tmpScript + 'set variable origin ligand\n '
    tmpScript = tmpScript + 'set variable size '+setVariableCutOffSize+'\n'
    tmpScript = tmpScript + 'set minchain 5\n'
    tmpScript = tmpScript + 'fix chains\n'
    tmpScript = tmpScript + 'cap all\n '
    tmpScript = tmpScript + 'set targetq 0\n '
    tmpScript = tmpScript + 'fix charges\n '
    
    # if len(HipLstOfResidues)!=0:
    #     lst = ''
    #     for ele in HipLstOfResidues:
    #         lst = lst + ':'+ele.strip()+' '
    #     tmpScript = tmpScript + 'set hip '+lst+'\n\n'


    # if len(HidLstOfResidues)!=0:
    #     lst = ''
    #     for ele in HidLstOfResidues:
    #         lst = lst + ':'+ele.strip()+' '
    #     tmpScript = tmpScript + 'set hid '+lst+'\n\n'

    tmpScript = tmpScript + lastPart.replace('aaaa',complexPDBName) 

    return tmpScript

def prepareReducedChopped(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,setVariableCutOffSize, HipLstOfResidues,HidLstOfResidues):

    print('Preparing Reduced Chopped System')

    # update the list of residues to be removed
    complexPDBName = complexPDB[:-4]

    chopScript = generateScript(complexPDB,ligandLstToRemoveFromPDB,residueToBeTheCenterOfTheChopping,setCapOriginAtom,setCutOffSize,setVariableCutOffSize, HipLstOfResidues,HidLstOfResidues)

    fileChopScript = open('chopScript.csh','w')
    fileChopScript.write(chopScript)
    fileChopScript.close()

    os.system('csh chopScript.csh')

def removeSolvent(complexPDB):

    print('REMOVING SOLVENT (HOH, SO4, GOL)')

    fileoutName = complexPDB[:-4]+'_no_solvent.pdb'

    ppdb = PandasPdb().read_pdb(complexPDB)

    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] != 'HOH']
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] != 'SO4']
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] != 'GOL']
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'HOH']
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'SO4']
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'GOL']

    ppdb.to_pdb(path=fileoutName, 
            records=None, 
            gz=False, 
            append_newline=True)
    
    return fileoutName

def fixPDBprotein(complexPDB,ligandResidueName):

    print('Fixing PDB')

    # remove chain in the ligand to avoid errors

    fileoutName = complexPDB[:-4]+'_fixed.pdb'

    fileout = open(fileoutName,'w')

    for line in open(complexPDB[:-4]+'_no_solvent.pdb'):
        newLine = line
        if ligandResidueName in line:
            newLine = line[:21]+' '+line[22:]
        fileout.write(newLine)

    fileout.close()
    
    return fileoutName

def prepare_zmats(zmat_arg, pdb_arg, resname_arg):

    # preliminary definitions 
    # please only change these if you know what you are doing!

    fixedComplexPDB = '' # it should look like name_fixed.pdb
    ligandZmat = '' #'benzene.z'

    titleOfYourSystem = '' # optional
    fixBackBoneSelection = [] #['4 67 70 74 110 144 152 153'] # If you have a lot of residues split the selection in different lines

    if zmat_arg:
        ligandZmat = str(zmat_arg)

    if pdb_arg:
        fixedComplexPDB = str(pdb_arg)

    if resname_arg:
        ligandResidueName = str(resname_arg)

    # *********************** CODE STARTS *********************************************


    print('PREPARING Z-MATRICES')

    # first delete the metals from the PDB file
    checkMetals(fixedComplexPDB, pdb_arg, resname_arg)

    prepareFinalZmatrixWithPEPz(fixedComplexPDB,titleOfYourSystem,ligandZmat,fixBackBoneSelection)

    createZmatrices(fixedComplexPDB,ligandResidueName, pdb_arg, resname_arg)

    relaxProteinLigand(fixedComplexPDB,mcproPath)

    generateFinalStructuresWithCAP(fixedComplexPDB) 

    addProteinBonds(fixedComplexPDB)

    # fixExludedAtomsList(fixedComplexPDB)

    #do some file management
    manageFiles(pdb_arg,resname_arg)   

def checkMetals(fixedComplexPDB, pdb_arg, resname_arg):

    print('Checking for Metals')

    ### THIS SHOULD BE HETATM I THINK ###
    ppdb = PandasPdb().read_pdb(fixedComplexPDB)
    hetatm_df = ppdb.df['HETATM']

    # could make this a for loop with the delete lists
    MG_count = hetatm_df[hetatm_df['residue_name']=='MG'].count()
    MN_count = hetatm_df[hetatm_df['residue_name']=='MN'].count()
    CA_count = hetatm_df[hetatm_df['residue_name']=='MN'].count()
    ZN_count = hetatm_df[hetatm_df['residue_name']=='ZN'].count()

    if (MG_count['residue_name'] != 0) or (MN_count['residue_name'] != 0) or (CA_count['residue_name'] != 0) or (ZN_count['residue_name'] != 0):
        print("METAL FOUND IN PDB, Z-MATRIX NOT AVAILABLE UNLESS YOU DELETE METAL AND RESUBMIT. \
            'CHOPPED' PDB IS ALREADY AVAILABLE IF THAT IS ALL USER REQUIRES.")
        manageFiles(pdb_arg, resname_arg)
        sys.exit()


def getNumberOfTheLastResidueOfTheChoppedSystem(choppedPDBName):
    
    resNumber = ''

    for line in open(choppedPDBName):
        if 'ATOM' in line[:6] or 'HETATM' in line[:6]:
            resNumber = line.split()[4]

    return resNumber 

def prepareFinalZmatrixWithPEPz(complexPDB,titleOfYourSystem,ligandZmat,fixBackBoneSelection):

    titleOfYourSystem = "PREPARED WITH PDB2ZMAT" # keep this for now

    complexPDBName = complexPDB[:-4]

    tmpfile = open(complexPDBName+'.all.in','w')

    for line in open(complexPDBName+'.chop.all.in'):
        newLine = line
        if '[ADD YOUR TITLE HERE]' in line: newLine = line.replace('[ADD YOUR TITLE HERE]',titleOfYourSystem)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: newLine = line.replace('[WRITE NAME OF YOUR solute z-matrix FILE]',ligandZmat)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: newLine = line.replace('[NAME OF THE z-matrix TO BE WRITTEN]',complexPDBName+'all.z')

        tmpfile.write(newLine)

    tmpfile.close()
        
    tmpfile = open(complexPDBName+'.var.in','w')

    for line in open(complexPDBName+'.chop.var.in'):
        newLine = line
        if '[ADD YOUR TITLE HERE]' in line: newLine = line.replace('[ADD YOUR TITLE HERE]',titleOfYourSystem)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: newLine = line.replace('[WRITE NAME OF YOUR solute z-matrix FILE]',ligandZmat)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: newLine = line.replace('[NAME OF THE z-matrix TO BE WRITTEN]',complexPDBName+'var.z')

        tmpfile.write(newLine)

    tmpfile.close()

    # create conrot

    lastResidue = getNumberOfTheLastResidueOfTheChoppedSystem(complexPDBName+'.chop.pdb')
    
    tmpfile = open(complexPDBName+'.var_conrot.in','w')

    #get the residues that need to be fixed for conrot
    res_strings = []
    for line in open(complexPDBName+'.var.in'):
        #get the initial strings of the fixed res
        if 'set fixed all' in line:
            res_strings = res_strings + line.split()[4:]

    #get all the numbers
    res_nums = []
    for res_str in res_strings:
        if '-' in res_str:
            flanks = res_str.split('-')
            res_nums = res_nums + list(range(int(flanks[0]),int(flanks[1])+1))
        else:
            res_nums = res_nums + [int(res_str)]

    #make the list of residues to fix the backbone of 
    fix_list = []
    for i in range(1,len(res_nums)):
        num = res_nums[i]
        prev_num = res_nums[i-1]
        diff = num - prev_num
        if (diff <= 5) and (diff != 1):
            fix_list = fix_list + list(range(prev_num+1,num))

    #put list in correct str format        
    fix_list_tmp = [str(x) for x in fix_list]
    fix_list = [' '.join(fix_list_tmp)]        


    for line in open(complexPDBName+'.var.in'):
        newLine = line
        if 'parameter type ALL *' in line: newLine = line+'$ set conrot\n$ set override domain 1-'+lastResidue+'\n' 
        if '$ set fixed backbone' in line:
            newLine = ''
            #for ele in fixBackBoneSelection: 
            for ele in fix_list:
                newLine = newLine + '$ set fixed backbone '+str(ele) + '\n'
        if '$ write zmatrix '+complexPDBName+'var.z' in line: newLine = line.replace(complexPDBName+'var.z',complexPDBName+'varcon.z')

        tmpfile.write(newLine)

    tmpfile.close()

def fixZmatrix(matrixZ,ligandResidueName, pdb_arg, resname_arg):

    print('Fixing Z matrix TERZ .... ',matrixZ)

    tmpfile = open('tmp.txt','w')

    try: 
        be2allFile = [ele for ele in open(matrixZ)]
    except:
        print("FAILED TO MAKE Z-MATRICES. PLEASE CHECK ABOVE ERROR MESSAGES. \n \
        (likely you have hetatoms that failed with BOSS) \n \
        NOTE THAT THE 'CHOPPED' PDB IS ALREADY PREPARED IF THAT IS ALL YOU DESIRE")
        manageFiles(pdb_arg, resname_arg)
        sys.exit()

    
    for iter in range(len(be2allFile)):
        newLine = be2allFile[iter]
        if 'TERZ' in newLine:
            if ligandResidueName in be2allFile[iter+1] or 'CAP' in be2allFile[iter+1]:
                newLine = be2allFile[iter]
            else:
                continue

        tmpfile.write(newLine)

    tmpfile.close()
                

    os.system('cp '+matrixZ+' lll.txt') 
    os.system('cp tmp.txt '+matrixZ)    
    
def createZmatrices(complexPDB,ligandResidueName, pdb_arg, resname_arg):

    MCPROscriptsPath = os.path.join(mcproPath,'scripts')

    complexPDBName = complexPDB[:-4]

    os.system(MCPROscriptsPath+'/xPEPZ '+complexPDBName+'.all')
    os.system(MCPROscriptsPath+'/xPEPZ '+complexPDBName+'.var')
    os.system(MCPROscriptsPath+'/xPEPZ '+complexPDBName+'.var_conrot')

    fixZmatrix(complexPDBName+'all.z',ligandResidueName, pdb_arg, resname_arg)        
    fixZmatrix(complexPDBName+'var.z',ligandResidueName, pdb_arg, resname_arg)        
    fixZmatrix(complexPDBName+'varcon.z',ligandResidueName, pdb_arg, resname_arg)     

def relaxProteinLigand(complexPDB,mcproPath):
    
    complexPDBName = complexPDB[:-4]

    os.system('mkdir CG9;cp '+complexPDBName+'all.z CG9;cd CG9;'+mcproPath+'/scripts/xCGDD9 '+complexPDBName+'all')
            

def getOptimizedZmatrix():

    # By default the file is called optsum 

    ZmatrixOptimized = []

    for line in open('CG9/optsum'):

        if 'Geometry Variations follow' in line: break

        ZmatrixOptimized.append(line)

    return ZmatrixOptimized

def replaceOptimizedCoordinatesInZmatrixFile(varFileName,capFileName,zmatrixOptimized):

    fileout = open(capFileName,'w')

    # get bottom part of the var file

    read = False
    bottomDataInVarFile = []
    
    for line in open(varFileName):
        if 'Geometry Variations follow' in line: read = True
        if read: bottomDataInVarFile.append(line)

    for ele in zmatrixOptimized: fileout.write(ele)
    for ele in bottomDataInVarFile: fileout.write(ele)

    
def generateFinalStructuresWithCAP(complexPDB):

    print('Generating final structures')

    complexPDBName = complexPDB[:-4]

    zmatrixOptimized = getOptimizedZmatrix()


    replaceOptimizedCoordinatesInZmatrixFile(complexPDBName+'var.z',complexPDBName+'cap.z',zmatrixOptimized)
    replaceOptimizedCoordinatesInZmatrixFile(complexPDBName+'varcon.z',complexPDBName+'capcon.z',zmatrixOptimized)
    
    try:
        os.mkdir('finalZmatrices')
    except:
        os.system('rm -r finalZmatrices')
        os.mkdir('finalZmatrices')

    os.system('cp CG9/optsum finalZmatrices/'+complexPDBName+'all.z')   
    os.system('cp '+complexPDBName+'cap.z '+ complexPDBName+'capcon.z finalZmatrices')

def addProteinBonds(complexPDB):
    ### based on the BONDADD.f script by Yue and Jonah ###

    complexPDBName = complexPDB[:-4]

    # read the all zmat in as optzmat
    with open('finalZmatrices/'+complexPDBName+'all.z') as optzmat:
        optzmat_data = optzmat.readlines()
    with open('finalZmatrices/'+complexPDBName+'capcon.z') as varzmat:
        varzmat_data = varzmat.readlines()

    # find relevant indicies to parse only important part of z-matrices
    line_index = 0
    for line in optzmat_data:

        if 'Tot. E' in line:
            first_atom_index = line_index+1
        if 'TERZ' in line:
            last_atom_index = line_index
            break # so that we only get first TERZ
        line_index = line_index + 1

    line_index = 0
    for line in optzmat_data:
        if 'Variable Bonds follow' in line:
            variable_bonds_index = line_index
        if 'Additional Bonds follow' in line:
            additional_bonds_index = line_index
            break
        line_index = line_index + 1

    # Get the atom numbers in the "variable bonds follow" list
    variable_atom_number_list = []
    for line in optzmat_data[variable_bonds_index:additional_bonds_index+1]:
        variable_atom_number_list.append(line[0:4])

    # print(variable_atom_number_list)

    # now check which lines these are connected to in the top of the zmat
    atom_pairs_list = []
    for line in optzmat_data[first_atom_index:last_atom_index]:
        # print(line[0:4],line[19:23])
        if line[0:4] in variable_atom_number_list:
            connected_atom = line[19:23]
            new_string = line[0:4] + connected_atom
            atom_pairs_list.append(new_string+'\n')

    # print(atom_pairs_list)

    # now add these pairs of atoms into the additional bonds of varzmat
    line_index = 0
    for line in varzmat_data:

        if 'Additional Bonds follow' in line:
            add_index = line_index

        line_index = line_index + 1

    varzmat_data = varzmat_data[:add_index+1] + atom_pairs_list + varzmat_data[add_index+1:]

    output_filename = 'finalZmatrices/'+complexPDBName+'_final_capcon_zmat.z'
    with open(output_filename,'w') as output_file:
        for line in varzmat_data:
            output_file.write(line)

# def fixExcludedAtomsList(complexPDB):

# 	complexPDBName = complexPDB[:-4]
# 	filename = 'finalZmatrices/'+complexPDBName+'_final_capcon_zmat.z'

# 	with open(filename) as zmat:
#         zmat_data = zmat.readlines()

#     index = 0
#     found = False
#     for line in zmat_data:
#     	if 'Excluded Atoms List follows' in line:
#     		found = True
#     		break
#     	index += 1

#     if found:
#     	if '-1' not in zmat_data[index+1][0:3]:




def manageFiles(original_pdb, ligand_resname):
    pdb_id = pdb_id = os.path.splitext(os.path.basename(original_pdb))[0][:-6] # [:-6] gets rid of '_fixed', which is used in zmat section
    ligand_id = str(ligand_resname).lower()

    try:
        os.makedirs(pdb_id + '_files')
        os.makedirs(ligand_id + '_files')
    except:
        print('folders already made')
        pass

    #move the pdb files (for some reaosn fixed is in the name here)
    pdb_output_files = [filename for filename in os.listdir('.') if filename.startswith(pdb_id)]
    # pdb_output_files.append(pdb_id + '.pdb')
    # pdb_output_files.append(pdb_id + '.pka')
    # pdb_output_files.append(pdb_id + '_no_solvent.pdb')
    # pdb_output_files.append(pdb_id + '.propka_input')
    #remove the folder name
    if (pdb_id + '_files') in pdb_output_files:
        pdb_output_files.remove(pdb_id + '_files')
    #mv these files to the output folder
    for file in pdb_output_files:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + '_files/')   
        except:
            pass

    #now same thing for the ligand
    ligand_output_files = [filename for filename in os.listdir('.') if filename.startswith(ligand_id)]
    #remove the folder name
    if (ligand_id + '_files') in ligand_output_files:
        ligand_output_files.remove(ligand_id + '_files')
    #mv these files to the output folder
    for file in ligand_output_files:
        try:
            os.system('mv ' + file + ' ./' + ligand_id + '_files/')   
        except:
            pass

    #take care of the other files
    try:
        os.makedirs(pdb_id + 'other_output')
    except:
        pass
    other_files_l = ['sum','log','out','lll.txt','plt.pdb','chopScript.csh','protonate.cmd','tmp.txt']
    for file in other_files_l:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + 'other_output/')
        except:
            pass

    # put it all in one folder
    os.system("mkdir "+pdb_id+'_folder')

    list_of_files = os.listdir(os.getcwd())

    excluded_files = ['.DS_Store', '.ipynb_checkpoints', '__init__.py', 'test_pdb2zmat.py', 'pdb_files', 'further_tests']
    folder_names = [filename for filename in os.listdir('.') if filename.endswith('_folder')]
    for folder in folder_names:
        excluded_files.append(folder)

    for file in list_of_files:

        if file.startswith('original_'):
            os.system("mv "+file+' ./pdb_files/'+pdb_id+'.pdb')

        elif file not in excluded_files:
            os.system("mv "+file+' ./'+pdb_id+'_folder')

if __name__ == "__main__":

    # get paths of needed dirs
    mcproPath = os.environ.get('MCPROdir')
    BOSSPath = os.environ.get('BOSSdir')
    BOSSscriptsPath = os.path.join(BOSSPath,'scripts')

    # run the main function
    main()

    