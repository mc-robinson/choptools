# !/usr/bin/env python
###############################################################################
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
#       It generates a folder called finalZmatrices with the final
#       z-matrix files (all, cap, cap + conrot for flexible backbones) 
# 
# Requirements:
#       BOSS
#       MCPRO
#       Reduce executable (http://kinemage.biochem.duke.edu/software/reduce.php)
#       Propka-3.1 executable (https://github.com/jensengroup/propka-3.1)
#       Babel
#       Python 2.7/3.6
#       
###############################################################################

import os
import sys
import subprocess
from biopandas.pdb import PandasPdb
import argparse

def main():
    '''
        The main function, this is the specific function that is the entry 
        point for the conda script. 
    '''
    # create parser object
    parser = argparse.ArgumentParser(
        prog='pdb2zmat_rewrite.py',
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
    Preferably Anaconda python 2.7/3.6 with following modules:
    biopandas
    Babel
    """
    )

    #defining arguments for parser object
    parser.add_argument(
        "-p", "--pdb", 
        help="name of the PDB file for the complex - i.e, [complex_name].pdb")
    parser.add_argument(
        "-z", "--zmat", 
        help="name of the zmat of the ligand with .z file descriptor. Only \
        need this if you want to import the optimized ligand zmat yourself.")
    parser.add_argument(
        "-r", "--resname",
        help="Residue name of the ligand from PDB FILE",
        type=str)
    parser.add_argument(
        "-c", "--cutoff",
        help="Size of the cutoff for chop to cut",
        type=str)
    parser.add_argument(
        "-v", "--variable_cutoff", 
        help="Size of the variable cutoff for chop to cut", 
        type=str)
    parser.add_argument(
        "-s", "--remove_solvent",
        help="Automatically attempt to remove solvent",
        action="store_true")
    parser.add_argument(
        "-t", "--tautomerize",
        help="Automatically assign Histidine Tautomerization states",
        action="store_true")
    parser.add_argument(
        "--hip",
        help="List of Histidine residues to make HIP. e.g. --hip 77A 56B ",
        nargs='*')
    parser.add_argument(
        "--hid",
        help="List of Histidine residues to make HID. e.g. --hid 77A 56B ",
        nargs='*')

    #parse the arguments from standard input
    args = parser.parse_args()

    # get paths of needed dirs
    MCPRO_path = os.environ.get('MCPROdir')
    BOSS_path = os.environ.get('BOSSdir')
    BOSS_scripts_path = os.path.join(BOSS_path,'scripts')

    # prepare the pdb of the complex using chop
    chopped_pdb_filename, ligand_zmat_filename = prepare_complex(args.zmat, 
                                                                args.pdb,
                                                                args.resname,
                                                                args.cutoff,
                                                                args.variable_cutoff,
                                                                args.tautomerize,
                                                                args.hip,
                                                                args.hid,
                                                                args.remove_solvent)

    # make the zmats from the fixed pdb
    # note that prepare_zmats takes in the FIXED PDB
    prepare_zmats(args.pdb, ligand_zmat_filename, chopped_pdb_filename, args.resname)

def prepare_complex(zmat_arg, complex_pdb_arg, ligand_resname_arg, 
                    cutoff_arg, variable_cutoff_arg, his_arg, hip_list, hid_list, remove_solvent_arg):
    """ Prepares the complex pdb for transfomation into a z-matrix
    
    The preparation of the complex involves several steps.
        - First, Checking for errors:
            * Check if BOSSdir and MCPROdir defined
            * Check that complex pdb only contains one ligand
        - Change the ligand res number for Clu to work properly
        - Preparing the LIGAND zmat
        - Removing the solvent from the complex pdb
        - Removing the chain ID of the ligand in the complex pdb
        - Merging the ligand with the complex pdb
        - [optional] Automatically handle histidine tautomers
        - Use chop to cut down the complex
    
    Args:
        zmat_arg: Filename of user supplied zmat file, if it was supplied
        pdb_arg: Filename of the user supplied pdb of the complex
        ligand_resname_arg: The residue name of ligand in the complex pdb
        cutoff_arg: The user specified cutoff for chopping
        his_arg: If the user would like to handle histidines automatically
        hip_list: List of Histidine residues to make HIP e.g. ['77a','56b']
        hid_list: List of Histidine residues to make HID e.g. ['77a','56b']
        remove_solvent_arg: Boolean. True if the program should automatically attempt to remove solven.


    Returns:
        chopped_complex_pdb: Filename of the prepared, chopped complex
        ligand_zmat: Filename of the prepared ligand zmat
    """

    ### preliminary definitions, that can be overidden 
    ### please only change these if you know what you are doing!

    complex_pdb_filename = '' # e.g. '1wbv.pdb'
    original_ligand_zmat = '' # e.g. 'benzene.z'
    ligand_resname = '' # e.g. 'LI3'

    ligand_list_to_remove_from_pdb = [] # recommend you start this empty
    residue_to_be_center_of_chopping = '' # e.g. 'LI3' # Normally the ligand
    cap_origin_atom = '' # e.g. 'LI3'  
    cutoff_size = '' # e.g. '18.0'

    title_of_system = '' # optional

    fixed_backbone_selection = [] # e.g. ['4 67 70 74 110 144 152 153'] 
    # If you have a lot of residues split the selection in different lines

    ### end of preliminary definitions

    # turn argparse arugments into strings 
    if zmat_arg:
        original_ligand_zmat_filename = str(zmat_arg)

    if complex_pdb_arg:
        complex_pdb_filename = str(complex_pdb_arg)

    if ligand_resname_arg:
        ligand_resname = str(ligand_resname_arg)
        residue_to_be_center_of_chopping = str(ligand_resname_arg) 
        cap_origin_atom = 'c01' # str(args.resname)

    if cutoff_arg:
        cutoff_size = str(cutoff_arg)
    else:
        cutoff_size = '20'

    if variable_cutoff_arg:
        variable_cutoff_size = str(variable_cutoff_arg)
    else:
        variable_cutoff_size = '10'

    # ********** MAIN CODE STARTS *********** # 

    ######## Why not just work with the data

    # get the paths for BOSS and MCPRO
    BOSS_path, MCPRO_path = check_paths(complex_pdb_filename)

    check_for_multiple_ligands(complex_pdb_filename, ligand_resname)

    # why am I not returning data
    complex_pdb_data = change_ligand_resnumber(complex_pdb_filename, 
                                                ligand_resname)

    # prepare the ligand zmat
    if zmat_arg:
        _, ligand_resnumber, ligand_original_resnumber = \
        generate_ligand_pdb(complex_pdb_data, ligand_resname)
        # use seperate function to fix the dummy atom names directly
        ligand_zmat_filename = fix_dummy_atom_names_directly(
                                                original_ligand_zmat_filename)
    else:
        ligand_zmat_filename, ligand_resnumber, ligand_original_resnumber = \
        prepare_ligand_zmat(complex_pdb_data, ligand_resname, MCPRO_path, 
                            BOSS_path)

   
    if remove_solvent_arg:
        no_solvent_pdb_data, no_solvent_pdb_filename = remove_solvent(complex_pdb_data, 
                                                                    complex_pdb_filename)
    else:
        no_solvent_pdb_data = complex_pdb_data
    
    # change Chain ID of ligand
    pre_clu_pdb_filename = fix_ligand_chain_id(no_solvent_pdb_data,
                                                complex_pdb_filename,
                                                ligand_resname)

    post_clu_pdb_filename = merge_ligand_with_protein(MCPRO_path, 
                                                    pre_clu_pdb_filename, 
                                                    ligand_zmat_filename, 
                                                    ligand_resname) # how is there no original resnumber here
    # get histdine lists before chopping
    # can also comment this section out if desired
    hip_residues_list = []  # resnumber, Chain ex. ['77a','56b'] #optional
    # hie_residues_list = []# resnumber, This the default!!! original code wrong
    hid_residues_list = []

    if his_arg:
        hip_residues_list, hid_residues_list, post_reduce_pdf_filename = \
        make_histidine_lists(post_clu_pdb_filename,
                                hip_residues_list,
                                hid_residues_list)

        print("List of HIP Residues from Reduce:", hip_residues_list)
        print("List of HID Residues from Reduce:", hid_residues_list)

        pre_chop_pdb_filename = post_reduce_pdf_filename
        tautomerized = True
    else:
        pre_chop_pdb_filename = post_clu_pdb_filename
        tautomerized = False

    if hip_list:
        for res in hip_list:
            if res not in hip_residues_list:
                hip_residues_list.append(res)
            if res in hid_residues_list:
                hid_residues_list.remove(res)
    if hid_list:
        for res in hid_list:
            if res not in hid_residues_list:
                hid_residues_list.append(res)
            if res in hip_residues_list:
                hip_residues_list.remove(res)

    print("List of HIP Residues:", hip_residues_list)
    print("List of HID Residues:", hid_residues_list)



    print('CHOP')
    chopped_pdb_filename = prepare_chopped_pdb(pre_chop_pdb_filename, 
                                                ligand_list_to_remove_from_pdb,
                                                residue_to_be_center_of_chopping ,
                                                cap_origin_atom, 
                                                cutoff_size,
                                                variable_cutoff_size,
                                                hip_residues_list,
                                                hid_residues_list,
                                                tautomerized)

    return chopped_pdb_filename, ligand_zmat_filename



def check_paths(complex_pdb_filename):
    """ Checks that files and BOSSdir and MCPROdir properly defined
    
    Args:
        complex_pdb_filename: String. Filename of the complex pdb

    Returns:
        BOSS_path: Full path to MCPROdir as defined in .bashrc
        MCPRO_path: Full path to MCPROdir as defined in .bashrc

    Raises:
        Errors for not having the path in the .bashrc file
        Promply exits after printing message.
        NOTE THAT THIS CURRENT STYLE IS FROWNED UPON, BUT WORKS
    """

    # attempt to get paths of needed dirs
    MCPRO_path = os.environ.get('MCPROdir')
    BOSS_path = os.environ.get('BOSSdir')

    if not MCPRO_path:
        print('Define MCPROdir enviroment variable, please')
        sys.exit(1)
    else:
        print('Using MCPROdir ...    ',MCPRO_path)

    if not BOSS_path:
        print('Define BOSSdir enviroment variable, please')
        sys.exit(1)
    else:
        print('Using BOSSdir ...    ',BOSS_path)

    if not os.path.isfile(complex_pdb_filename):
        print('PDB not found ...   ',complex_pdb_filename)
        sys.exit(1)

    return BOSS_path, MCPRO_path

def check_for_multiple_ligands(complex_pdb_filename, ligand_resname):
    """ Checks that supplied pdb only contains one copy of desired ligand
    
    Args:
        complex_pdb_filename: String. Filename of the complex pdb
        ligand_resname: String. Residue name of desired ligand

    Raises:
        Errors for having more than one copy of the ligand
        Promply exits after printing message.
        NOTE THAT THIS CURRENT STYLE IS FROWNED UPON, BUT WORKS
    """

    # create a BioPandas object
    ppdb = PandasPdb()
    ppdb.read_pdb(complex_pdb_filename) 
    atom_df = ppdb.df['ATOM']
    hetatm_df = ppdb.df['HETATM']

    # get only information about the ligand
    lig_df = hetatm_df.loc[hetatm_df['residue_name'] == ligand_resname]

    # check if lig_df contains multiple chain ids (contains multiple ligands)
    if lig_df.chain_id.nunique() > 1:
        print("THERE IS MORE THAN ONE COPY OF THE DESIRED LIGAND." \
            "PLEASE DELETE EXTRAS")
        sys.exit(1)

def change_ligand_resnumber(complex_pdb_filename, ligand_resname):
    """ Changes the ligand residue number in the original pdb
    
    For some reason, CLU seems to need the residue number to be 100
    Therefore, this code goes through and does that.
    
    Args:
        complex_pdb_filename: String. Filename of the complex pdb
        ligand_resname: String. Residue name of desired ligand

    Returns:
        complex_pdb_data: List. The updated data from the complex pdb file
    """
    
    print('CHANGING LIGAND RESIDUE NUMBER') #needs to be 100 to work with Clu

    # Read in the file
    with open(complex_pdb_filename, 'r') as complex_pdb_file:
        complex_pdb_data = complex_pdb_file.readlines()

    # Replace the res number with 100 for clu to work
    idx = 0
    for line in complex_pdb_data:
        if 'ATOM' in line[:7] or 'HETATM' in line[:7]:
            if ligand_resname in line:

                new_line = line[:21]+' '+line[22:] # get rid of chain ID
                resnumber = new_line[12:29].split()[2]
                ligand_resnumber = resnumber

                if int(resnumber)<=999:
                    ligand_resnumber = '100'
                    new_line = new_line[:23] + '100' + new_line[26:]

                elif int(resnumber)>999:
                    ligand_resnumber = ' 100'
                    # note the extra space we include here
                    new_line = new_line[:22] + ' 100' + new_line[26:]
            else:
                new_line = line
        else:
            new_line = line

        complex_pdb_data[idx] = new_line
        idx = idx + 1

    # return the updated data
    return complex_pdb_data

def prepare_ligand_zmat(complex_pdb_data, ligand_resname, MCPRO_path, 
                        BOSS_path):
    """ Prepares the ligand zmat, starting from the complex pdb

    The process involves extracting the ligand pdb info from the complex pdb,
    then protonating it and converting it to a z-matrix. Lastly, it's optimized.
    
    Args:
        complex_pdb_data: List. Data from the complex pdb file
        ligand_resname: String. Residue name of desired ligand
        BOSS_path: Full path to MCPROdir as defined in .bashrc
        MCPRO_path: Full path to MCPROdir as defined in .bashrc

    Returns:
        ligand_zmat_filename: String. Filename of the ligand zmat
        ligand_resnumber: String. The resnumber of the lig after changes (100)
        ligand_original_resnumber: String. Resnumber of lig in complex_pdb_data
    """

    print('Preparing Ligand Z matrix')

    ligand_pdb_filename, ligand_resnumber, ligand_original_resnumber = \
    generate_ligand_pdb(complex_pdb_data, ligand_resname)

    protonated_ligand_pdb_filename = \
    protonate_ligand_with_Babel(ligand_pdb_filename)

    # Convert pdb file to BOSS zmat using BOSS
    convert_pdb_to_zmat(protonated_ligand_pdb_filename, BOSS_path)

    # fix the names of the dummy atoms in the newly created zmat
    ligand_zmat_filename = fix_dummy_atom_names(protonated_ligand_pdb_filename)
    
    optimize_zmat(ligand_zmat_filename, MCPRO_path)

    return ligand_zmat_filename, ligand_resnumber, ligand_original_resnumber

def generate_ligand_pdb(complex_pdb_data, ligand_resname):
    """ Generates the pdb file for the ligand by extracting ligand information
    from the pdb file of the complex.
    
    Args:
        complex_pdb_data: List. Data from the complex pdb file
        ligand_resname: String. Residue name of desired ligand

    Returns:
        ligand_pdb_filename: String. The ligand pdb filename
        ligand_resnumber: String. The resnumber of the lig after changes (100)
        ligand_original_resnumber: String. Resnumber of lig in complex_pdb_data
    """

    print('Generating Ligand PDB')

    ligand_pdb_filename = ligand_resname.lower()+'.pdb'
    ligand_pdb_file = open(ligand_pdb_filename,'w')

    # the rest of this is rather sloppy an could be tightened up

    ligand_resnumber = ''
    ligand_original_resnumber = ''

    for line in complex_pdb_data:
        if 'ATOM' in line[:7] or 'HETATM' in line[:7]:
            if ligand_resname in line[15:23]:
                
                new_line = line[:21]+' '+line[22:]
                resnumber = new_line[12:29].split()[2]
                ligand_resnumber = resnumber
                ligand_original_resnumber = resnumber

                if int(resnumber)<=999:
                    ligand_resnumber = '100'
                    new_line = new_line[:23] + '100' + new_line[26:]
                    ligand_pdb_file.write(new_line)

                elif int(resnumber)>999:
                    ligandResnumber = ' 100'
                    new_line = new_line[:22] + '0100' + new_line[26:]
                    ligand_pdb_file.write(new_line)

                else:
                    ligand_pdb_file.write(new_line)

    ligand_pdb_file.close()

    return ligand_pdb_filename, ligand_resnumber, ligand_original_resnumber

def protonate_ligand_with_Babel(ligand_pdb_filename):
    """ Uses Babel to protonate the ligand pdb file
    
    Args:
        ligand_pdb_filename: String. Name of ligand pdb file

    Returns:
        protonated_ligand_pdb_filename: String. ligand_pdb_filename + '_h'

    Raises:
        Error if Babel fails. Promptly exits
    """

    print('Protonating ligand with Babel')

    protonated_ligand_pdb_filename = ligand_pdb_filename[:-4]+'_h.pdb'

    # chimeraScriptFile = open('protonate.cmd','w')
    # chimeraScriptFile.write(chimeraScript.replace('aaaa',namePDBLigand) \ 
    # .replace('bbbb',namePDBLigand_protonated))
    # chimeraScriptFile.close()

    #cmd = '/Applications/Chimera.app/Contents/Resources/bin/chimera \
    # --nogui protonate.cmd'
    # cmd = 'chimera --nogui protonate.cmd'

    protonate_cmd = ('babel %s -O %s -p' % (ligand_pdb_filename, 
                                    protonated_ligand_pdb_filename))

    try:
        os.system(protonate_cmd)
    # bad error handling style below, but works for now
    except Exception as e:
        print(e)
        print("Babel failed to perform the command:", protonate_cmd)
        sys.exit(1)

    return protonated_ligand_pdb_filename

def convert_pdb_to_zmat(protonated_ligand_pdb_filename, BOSS_path):
    """ Converts the ligand pdb to a BOSS z-matrix
    
    Args:
        protonated_ligand_pdb_filename: String. name of ligand pdb file
    
    Raises:
        Error if the BOSS pdb to zmat conversion command fails
    """

    print('Converting ligand pdb to BOSS zmat')

    convert_to_zmat_cmd = \
    BOSS_path + '/scripts/xPDBMCP ' + protonated_ligand_pdb_filename[:-4]

    try:
        os.system(convert_to_zmat_cmd)
    except Exception as e:
        print(e)
        print("BOSS failed to perform the command:", convert_to_zmat_cmd)
        sys.exit(1)  

def fix_dummy_atom_names(protonated_ligand_pdb_filename):
    """ Fixes the dummy atom names in the ligand zmat produced by BOSS

    The issue here is that ?????????
    
    Args:
        protonated_ligand_pdb_filename: String. name of ligand pdb file + '_h'

    Returns:
        ligand_zmat_filename: String. Name of ligand zmat file. Note this is
                                protonated_ligand_pdb_filename with .z, not .pdb
    """
    
    print('Fixing DUMMY atoms names')

    ligand_zmat_filename = protonated_ligand_pdb_filename[:-4]+'.z'

    # copy over zmat data to a tmp file and change dummy names in process
    tmp_file = open('tmp.txt','w')
    for line in open(ligand_zmat_filename):
        new_line = line       
        if len(line.split())==12 and len(line)==77:
            new_line = line[:71]+'    1\n'
        tmp_file.write(new_line.replace('   1 DUM','   1 DU1') \
                        .replace('   2 DUM','   2 DU2') \
                        .replace('   3 O00  800  800','   3 DU3   -1    0'))
    tmp_file.close()

    # copy the file over to desired filename
    copy_cmd = 'cp tmp.txt ' + ligand_zmat_filename 
    os.system(copy_cmd)

    return ligand_zmat_filename

def fix_dummy_atom_names_directly(original_ligand_zmat_filename):

    ### I am pretty sure this function is unnecessary)!!!!! ###

    print('Fixing DUMMY atoms names directly')

    ligand_zmat_filename = original_ligand_zmat_filename[:-2] + '_fixed.z'

    # copy over zmat data to a tmp file and change dummy names in process
    tmp_file = open('tmp.txt','w')
    for line in open(ligand_zmat_filename):
        new_line = line       
        if len(line.split())==12 and len(line)==77:
            new_line = line[:71]+'    1\n'
        tmp_file.write(new_line.replace('   1 DUM','   1 DU1') \
                        .replace('   2 DUM','   2 DU2') \
                        .replace('   3 O00  800  800','   3 DU3   -1    0'))
    tmp_file.close()

    # copy the file over to desired filename
    copy_cmd = 'cp tmp.txt ' + ligand_zmat_filename 
    os.system(copy_cmd)

    return ligand_zmat_filename

def optimize_zmat(ligand_zmat_filename, MCPRO_path):
    """ Optimizes the zmat created from BOSS, once Dummy atoms have been fixed
    
    Args:
        protonated_ligand_pdb_filename: String. name of ligand pdb file + '_h'

    Returns:
        ligand_zmat_filename: String. Name of ligand zmat file. Note this is
                                protonated_ligand_pdb_filename with .z, not .pdb
    
    Raises:
        Error if the MCPRO optimization command fails
    """

    print('Optimizing Structure')

    optimize_cmd = MCPRO_path+'/scripts/xOPT '+ligand_zmat_filename[:-2]

    try:
        os.system(optimize_cmd)
    # bad error handling style below, but works for now
    except Exception as e:
        print(e)
        print("MCPRO failed to perform the command:", optimize_cmd)
        sys.exit(1)

def remove_solvent(complex_pdb_data, complex_pdb_filename):
    """ Removes solvents (HOH, SO4, GOL) from the pdb file

    Note that this function does a bit of work creating and deleting files 
    in order to work with BioPandas. The goal is to only work with the data
    until strictly necessary to switch to files. This function thus cleans the
    files as it goes.
    
    Args:
        complex_pdb_data: List. Data from the complex pdb file
        complex_pdb_filename: String. Filename of the complex pdb

    Returns:
        no_solvent_pdb_data: List. Data from the pdb with no solvent
    
    """

    print('REMOVING SOLVENTS (HOH, SO4, GOL)')

    # create a temporary file
    tmp_filename = 'no_solvent_tmp.pdb'
    tmp_file = open(tmp_filename,'w')
    for line in complex_pdb_data:
        tmp_file.write(line + '\n')
    tmp_file.close

    fileout_name = complex_pdb_filename[:-4]+'_no_solvent.pdb'

    ppdb = PandasPdb().read_pdb(tmp_filename)

    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] != 'HOH']
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] != 'SO4']
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] != 'GOL']
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'HOH']
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'SO4']
    ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name'] != 'GOL']

    #print(ppdb.df['HETATM']['residue_name'])

    ppdb.to_pdb(path=fileout_name, 
            records=None, 
            gz=False, 
            append_newline=True)

    # read the new file
    with open(fileout_name, 'r') as no_solvent_pdb_file:
        no_solvent_pdb_data = no_solvent_pdb_file.readlines()

    # clean up the created pdb files
    os.system('rm ' + tmp_filename)
    # os.system('rm ' + fileout_name)
    
    return no_solvent_pdb_data, fileout_name

def fix_ligand_chain_id(no_solvent_pdb_data, complex_pdb_filename, 
                        ligand_resname):
    """ Removes the ligand chain ID, which is necessary for working with clu
    and chop.
    
    Args:
        no_solvent_pdb_data: List. Data from the complex pdb file with no solvent
        complex_pdb_filename: String. Filename of the complex pdb
        ligand_resname: String. Residue name of desired ligand

    Returns:
        pre_club_pdb_filename: String. Filename of the complex pdb prepared
                                        for merging with Clu.
    
    """

    print('Fixing ligand chain id of PDB')

    # remove chain id in the ligand to avoid errors
    pre_clu_pdb_filename = complex_pdb_filename[:-4]+'_pre_clu.pdb'
    fileout = open(pre_clu_pdb_filename,'w')

    for line in no_solvent_pdb_data:
        new_line = line
        if ligand_resname in line:
            new_line = line[:21]+' '+line[22:]
        fileout.write(new_line)

    fileout.close()
    
    return pre_clu_pdb_filename

def merge_ligand_with_protein(MCPRO_path, pre_clu_pdb_filename, 
                                ligand_zmat_filename, ligand_resname):
    """ Uses clu to replace the original ligand in the complex with the 
    protonated, optimized ligand.

    Args:
        MCPRO_path: MCPRO_path: String. Full path to MCPROdir defined in .bashrc
        pre_clu_pdb_filename: String. Filename of the complex pdb pre clu
        ligand_zmat_filename: String. Filename of the ligand zmat
        ligand_resname: String. Residue name of the ligand

    Returns:
        post_clu_pdb_filename: String. Filename of the complex PDB after merging with clu.
    
    """

    #clu -t:s=5001 2be2.pdb -r r22_h.z -n 2be2_cplx.pdb

    print('Merging ligand with PDB Protein')
    
    clu_path = MCPRO_path +'/miscexec/clu '

    post_clu_pdb_filename = pre_clu_pdb_filename.replace('_pre_clu.pdb','_post_clu.pdb')
    if os.path.isfile(post_clu_pdb_filename): 
        os.remove(post_clu_pdb_filename)

    # THIS RESNUMBER MUST BE THE SAME AS BEOFRE 
    # cmd = clu + ' -t:s='+resnumber+' '+complexPDB+' -r '+ligandZmat+' -n '+outputPDB
    clu_cmd = \
    clu_path + ' -t:s=:'+ligand_resname+' '+pre_clu_pdb_filename+ \
    ' -r '+ligand_zmat_filename+' -n '+post_clu_pdb_filename
    print(clu_cmd)

    try:
        os.system(clu_cmd)
    except Exception as e:
        print(e)
        print("MCPRO failed to perform the command:", clu_cmd)
        sys.exit(1) 

    return post_clu_pdb_filename

def make_histidine_lists(post_clu_pdb_filename, hip_residues_list, hid_residues_list):
    """
    Uses reduce to decide the correct histidine tautomerization states of the complex.

    Args:
        post_clu_pdb_filename: String. Filename of the complex PDB after merging with clu.
        hip_residues_list: List. List of residue numbers for HIP residues, ex. ['77a','56b']
        hid_residues_list: List. List of residue numbers for HIE residues, ex. ['77a','56b']
    Returns:
        reduce_output_filename: String. Filename of the complex PDB after merging.
    """

    pdb_input_filename = post_clu_pdb_filename
    reduce_output_filename = post_clu_pdb_filename[:-4]+'_h.pdb'

    # run reduce on the protein
    #subprocess.call("reduce –build " + pdb_name + " > " + pdb_name[:-4]+ "_h.pdb", shell=True)
    try:
        os.system("reduce -Build %s > %s" % (pdb_input_filename, reduce_output_filename))
    except Exception as e:
        print(e)
        print("reduce failed on the command:", ("reduce -Build %s > %s" % (pdb_input_filename, reduce_output_filename)))
        sys.exit()
        
    
    #read in the protein with biopandas
    ppdb = PandasPdb()
    ppdb.read_pdb(reduce_output_filename)
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

    print('HIP list')
    print(hip_list)
    print('HIE list')
    print(hie_list)
    print('HID list')
    print(hid_list)

    #add to the global lists
    for res in hip_list:
        if (res not in hip_residues_list) and (res not in hid_residues_list):
            hip_residues_list.append(res)

    for res in hid_list:
        if (res not in hip_residues_list) and (res not in hid_residues_list):
            hid_residues_list.append(res)

    #check if propka disagrees with reduce
    run_propka(post_clu_pdb_filename, hip_residues_list)

    return hip_residues_list, hid_residues_list, reduce_output_filename     

def run_propka(post_clu_pdb_filename, hip_residues_list):
    """
    Uses propka to show any residues for which propka and reduce disagree on the
    histidine tautomerization state.

    Args:
        post_clu_pdb_filename: String. Filename of the complex PDB after merging with clu.
        hip_residues_list: List. List of residue numbers for HIP residues, ex. ['77a','56b']
    """

    original_pdb_filename = post_clu_pdb_filename.replace('_post_clu.pdb','.pdb')

    try:
        os.system('propka31 ' + original_pdb_filename)
    except Exception as e:
        print(e)
        print('propka failed on the command:', ('propka31 ' + original_pdb_filename))
        sys.exit()

    propka_output_name = original_pdb_filename[:-4]+'.pka'

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
        if res not in hip_residues_list:
            print("prokpka thinks " + res + " should be HIP, reduce does not")

    for res in hip_residues_list:
        if res not in propka_hip_list:
            print("reduce thinks " + res + " should be HIP, propka does not")

def prepare_chopped_pdb(pre_chop_pdb_filename, 
                        ligand_list_to_remove_from_pdb,
                        residue_to_be_center_of_chopping ,
                        cap_origin_atom, 
                        cutoff_size,
                        variable_cutoff_size,
                        hip_residues_list,
                        hid_residues_list,
                        tautomerized):
    """
    Creates a chop script with all of the directions for chopping the protein.

    Args:
        pre_chop_pdb_filename: String. Filename of the complex PDB after clu and/or reduce
        ligand_list_to_remove_from_pdb: List. Ligands to be removed from the protein.
        residue_to_be_center_of_chopping: String. Normally the ligand, e.g. 'LI3' 
        cap_origin_atom: String. Normally 'c01'
        cutoff_size: String. Cutoff for the whole protein. e.g. '20'
        variable_cutoff_size: String. Cutoff for the cutoff regions. e.g. '10'
        hip_residues_list: List. List of residue numbers for HIP residues, e.g. ['77a','56b']
        hid_residues_list: List. List of residue numbers for HIE residues, e.g. ['77a','56b']
        tautomerized: Boolean. True if histidines were added with Reduce.

    Returns:
        chopped_pdb_filename: String. Filename of the protein post chopping. 
    """

    print('Preparing Chopped System')

    # update the list of residues to be removed
    # complexPDBName = complexPDB[:-4]

    chop_script = generate_chop_script(pre_chop_pdb_filename, 
                                        ligand_list_to_remove_from_pdb,
                                        residue_to_be_center_of_chopping ,
                                        cap_origin_atom, 
                                        cutoff_size,
                                        variable_cutoff_size,
                                        hip_residues_list,
                                        hid_residues_list,
                                        tautomerized)

    chop_script_file = open('chop_script.csh','w')
    chop_script_file.write(chop_script)
    chop_script_file.close()

    try:
        os.system('csh chop_script.csh')
    except Exception as e:
        print(e)
        print('chop failed on the command: csh chop_script.csh')
        sts.exit()

    if tautomerized:
        chopped_pdb_filename = pre_chop_pdb_filename.replace('_post_clu_h.pdb','.chop.pdb')
    else:
        chopped_pdb_filename = pre_chop_pdb_filename.replace('_post_clu.pdb','.chop.pdb')

    return chopped_pdb_filename


def generate_chop_script(pre_chop_pdb_filename, 
                        ligand_list_to_remove_from_pdb,
                        residue_to_be_center_of_chopping ,
                        cap_origin_atom, 
                        cutoff_size,
                        variable_cutoff_size,
                        hip_residues_list,
                        hid_residues_list,
                        tautomerized):
    
    """
    Creates a chop script with all of the directions for chopping the protein.

    Args:
        pre_chop_pdb_filename: String. Filename of the complex PDB after clu and/or reduce
        ligand_list_to_remove_from_pdb: List. Ligands to be removed from the protein.
        residue_to_be_center_of_chopping: String. Normally the ligand, e.g. 'LI3' 
        cap_origin_atom: String. Normally 'c01'
        cutoff_size: String. Cutoff for the whole protein. e.g. '20'
        variable_cutoff_size: String. Cutoff for the cutoff regions. e.g. '10'
        hip_residues_list: List. List of residue numbers for HIP residues, e.g. ['77a','56b']
        hid_residues_list: List. List of residue numbers for HIE residues, e.g. ['77a','56b']
        tautomerized: Boolean. True if histidines were added with Reduce.

    Returns:
        chop_script: String. A long string that can be written to a script file. 
    """
    

    script_ending = \
''' 
write pdb aaaa.chop.pdb
write pepz all aaaa.chop.all.in
write pepz variable aaaa.chop.var.in
write translation aaaa.chop.tt
EOF

'''

    # complexPDBName = complexPDB[:-4]

    tmp_script = '$MCPROdir/miscexec/chop -t -i '+pre_chop_pdb_filename+' << EOF\n'

    if len(ligand_list_to_remove_from_pdb)!=0:
        for ele in ligand_list_to_remove_from_pdb:
            tmp_script = tmp_script + 'delete ligand :'+ele+'\n'

    if len(hip_residues_list)!=0:
        lst = ''
        for ele in hip_residues_list:
            lst = lst + ':'+ele.strip()+' '
        tmp_script = tmp_script + 'set hip '+lst+'\n\n'


    if len(hid_residues_list)!=0:
        lst = ''
        for ele in hid_residues_list:
            lst = lst + ':'+ele.strip()+' '
        tmp_script = tmp_script + 'set hid '+lst+'\n\n'

    tmp_script = tmp_script + 'add center :'+residue_to_be_center_of_chopping+'\n'
    tmp_script = tmp_script + 'set cap origin :'+cap_origin_atom+'\n'

    tmp_script = tmp_script + 'set cut origin ligand \n'
    tmp_script = tmp_script + 'set cut size '+cutoff_size+'\n'
    #tmpScript += 'delete cut :gol\n'
    if len(ligand_list_to_remove_from_pdb)!=0:
        for ele in ligand_list_to_remove_from_pdb:
            tmp_script = tmp_script + 'delete cut :'+ele.lower()+'\n'

    tmp_script = tmp_script + 'set variable origin ligand\n '
    tmp_script = tmp_script + 'set variable size '+variable_cutoff_size+'\n'
    tmp_script = tmp_script + 'set minchain 5\n'
    tmp_script = tmp_script + 'fix chains\n'
    tmp_script = tmp_script + 'cap all\n '
    tmp_script = tmp_script + 'set targetq 0\n '
    tmp_script = tmp_script + 'fix charges\n '
    
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
    if tautomerized:
        tmp_script = tmp_script + script_ending.replace('aaaa',pre_chop_pdb_filename.replace('_post_clu_h.pdb','')) 
    else:
        tmp_script = tmp_script + script_ending.replace('aaaa',pre_chop_pdb_filename.replace('_post_clu.pdb','')) 
    chop_script = tmp_script

    return chop_script

def prepare_zmats(complex_pdb_arg, ligand_zmat_filename, chopped_pdb_filename, ligand_resname_arg):
    """
    Prepares the final zmatrices from the chopped PDB and the ligand zmat

    Args:
        complex_pdb_arg: Filename of the original pdb complex
        ligand_zmat_filename: String. Filename of the ligand zmat
        chopped_pdb_filename: String. Filename of the protein post chopping. 
        ligand_resname_arg: The residue name of ligand in the complex pdb

    """
    if ligand_resname_arg:
        ligand_resname = str(ligand_resname_arg)

    if complex_pdb_arg:
        original_pdb_filename = str(complex_pdb_arg)

    ### ********************* preliminary definitions ********************************
    # please only change these if you know what you are doing!

    system_title = '' # optional
    fixed_backbone_residues = [] # e.g. ['4 67 70 74 110 144 152 153']
                                 # If you have a lot of residues split the selection in different lines

    # *********************** CODE STARTS *********************************************


    print('PREPARING Z-MATRICES')

    # first delete the metals from the PDB file
    check_metals(chopped_pdb_filename, original_pdb_filename, ligand_resname)

    prepare_final_zmats_for_PEPZ(original_pdb_filename, chopped_pdb_filename, system_title, 
                                ligand_zmat_filename, fixed_backbone_residues)

    create_zmats(chopped_pdb_filename, ligand_resname, original_pdb_filename)

    relax_protein_ligand_system(original_pdb_filename)

    generate_final_structures_with_cap(original_pdb_filename) 

    add_protein_bonds(original_pdb_filename)

    # fixExludedAtomsList(fixedComplexPDB)

    #do some file management
    manage_files(original_pdb_filename,ligand_resname)   

def check_metals(chopped_pdb_filename, original_pdb_filename, ligand_resname):
    """
    Checks the chopped pdb file for metals, which MCPRO will not be able to deal with.

    Args:
        chopped_pdb_filename: String. Filename of the chopped pdb complex
        original_pdb_filename: String. Filename of the original pdb complex
        ligand_resname: String. The residue name of ligand in the complex pdb

    """

    print('Checking for Metals')

    ### THIS SHOULD BE HETATM I THINK ###
    ppdb = PandasPdb().read_pdb(chopped_pdb_filename)
    hetatm_df = ppdb.df['HETATM']

    # could make this a for loop with the delete lists
    MG_count = hetatm_df[hetatm_df['residue_name']=='MG'].count()
    MN_count = hetatm_df[hetatm_df['residue_name']=='MN'].count()
    CA_count = hetatm_df[hetatm_df['residue_name']=='CA'].count()
    ZN_count = hetatm_df[hetatm_df['residue_name']=='ZN'].count()

    if (MG_count['residue_name'] != 0) or (MN_count['residue_name'] != 0) or \
        (CA_count['residue_name'] != 0) or (ZN_count['residue_name'] != 0):

        print("METAL FOUND IN PDB, Z-MATRIX NOT AVAILABLE UNLESS YOU DELETE METAL AND RESUBMIT. \
            'CHOPPED' PDB IS ALREADY AVAILABLE IF THAT IS ALL USER REQUIRES.")
        manage_files(original_pdb_filename, ligand_resname)
        sys.exit()

def prepare_final_zmats_for_PEPZ(original_pdb_filename, chopped_pdb_filename, system_title, 
                                ligand_zmat_filename, fixed_backbone_residues):
    """
    Prepares all the .in files for the .var, .all. and .conrot z-matrices to be made

    Args:
        original_pdb_filename: String. Filename of the original pdb complex
        chopped_pdb_filename: String. Filename of the chopped pdb complex
        system_title: String. Title to be put at the top of all z-matrices.
        ligand_zmat_filename: String. Filename of the ligand zmat.
        fixed_backbone_residues: List. Residues to keep fixed in the z-matrices.
    """

    # original_pdb_filename = chopped_pdb_filename.replace('.chop.pdb','.pdb')

    system_title = "PREPARED WITH PDB2ZMAT" # keep this for now

    tmpfile = open(original_pdb_filename[:-4]+'.all.in','w')

    for line in open(original_pdb_filename[:-4]+'.chop.all.in'):
        new_line = line
        if '[ADD YOUR TITLE HERE]' in line: 
            new_line = line.replace('[ADD YOUR TITLE HERE]', system_title)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: 
            new_line = line.replace('[WRITE NAME OF YOUR solute z-matrix FILE]',ligand_zmat_filename)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: 
            continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: 
            new_line = line.replace('[NAME OF THE z-matrix TO BE WRITTEN]',original_pdb_filename[:-4]+'.all.z')

        tmpfile.write(new_line)

    tmpfile.close()
        
    tmpfile = open(original_pdb_filename[:-4]+'.var.in','w')

    for line in open(original_pdb_filename[:-4]+'.chop.var.in'):
        new_line = line
        if '[ADD YOUR TITLE HERE]' in line: 
            new_line = line.replace('[ADD YOUR TITLE HERE]', system_title)
        if '[WRITE NAME OF YOUR solute z-matrix FILE]' in line: 
            new_line = line.replace('[WRITE NAME OF YOUR solute z-matrix FILE]',ligand_zmat_filename)
        if '[NAME OF pdb file TO BE WRITTEN]' in line: 
            continue
        if '[NAME OF THE z-matrix TO BE WRITTEN]' in line: 
            new_line = line.replace('[NAME OF THE z-matrix TO BE WRITTEN]',original_pdb_filename[:-4]+'.var.z')

        tmpfile.write(new_line)

    tmpfile.close()

    # create conrot
    last_residue = get_number_of_last_residue(chopped_pdb_filename)
    
    tmpfile = open(original_pdb_filename[:-4]+'.var_conrot.in','w')

    # get the residues that need to be fixed for conrot
    res_strings = []
    for line in open(original_pdb_filename[:-4]+'.var.in'):
        #get the initial strings of the fixed res
        if 'set fixed all' in line:
            res_strings = res_strings + line.split()[4:]

    # get all the numbers
    res_nums = []
    for res_str in res_strings:
        if '-' in res_str:
            flanks = res_str.split('-')
            res_nums = res_nums + list(range(int(flanks[0]),int(flanks[1])+1))
        else:
            res_nums = res_nums + [int(res_str)]

    # make the list of residues to fix the backbone of 
    fix_list = []
    for i in range(1,len(res_nums)):
        num = res_nums[i]
        prev_num = res_nums[i-1]
        diff = num - prev_num
        if (diff <= 5) and (diff != 1):
            fix_list = fix_list + list(range(prev_num+1,num))

    # put list in correct str format        
    fix_list_tmp = [str(x) for x in fix_list]
    fix_list = [' '.join(fix_list_tmp)]        


    for line in open(original_pdb_filename[:-4]+'.var.in'):
        new_line = line
        if 'parameter type ALL *' in line: 
            new_line = line+'$ set conrot\n$ set override domain 1-'+last_residue+'\n' 
        if '$ set fixed backbone' in line:
            new_line = ''
            #for ele in fixBackBoneSelection: 
            for ele in fix_list:
                new_line = new_line + '$ set fixed backbone '+str(ele) + '\n'
        if '$ write zmatrix '+original_pdb_filename[:-4]+'.var.z' in line: 
            new_line = line.replace(original_pdb_filename[:-4]+'.var.z',original_pdb_filename[:-4]+'.varcon.z')

        tmpfile.write(new_line)

    tmpfile.close()

def get_number_of_last_residue(pdb_filename):
    """
    Gets the residue number of the last residue in the pdb file.

    Args:
        pdb_filename: String. Filename of the pdb.
    """
    res_number = ''

    for line in open(pdb_filename):
        if 'ATOM' in line[:6] or 'HETATM' in line[:6]:
            res_number = line.split()[4]

    return res_number 

def create_zmats(chopped_pdb_filename, ligand_resname, original_pdb_filename):
    """
    Creates the z-matrices using xPEPZ from MCPRO.

    Args:
        chopped_pdb_filename: String. Filename of the chopped pdb complex
        ligand_resname: String. Residue name of the ligand
        original_pdb_filename: String. Filename of the original PDB.
    """

    MCPRO_scripts_path = os.path.join(os.environ.get('MCPROdir'),'scripts')

    os.system(MCPRO_scripts_path+'/xPEPZ '+original_pdb_filename[:-4]+'.all')
    os.system(MCPRO_scripts_path+'/xPEPZ '+original_pdb_filename[:-4]+'.var')
    os.system(MCPRO_scripts_path+'/xPEPZ '+original_pdb_filename[:-4]+'.var_conrot')

    fix_zmat(original_pdb_filename[:-4]+'.all.z', ligand_resname, original_pdb_filename)        
    fix_zmat(original_pdb_filename[:-4]+'.var.z', ligand_resname, original_pdb_filename)        
    fix_zmat(original_pdb_filename[:-4]+'.varcon.z', ligand_resname, original_pdb_filename)   

def fix_zmat(zmat_filename, ligand_resname, original_pdb_filename):
    """
    Need to remember why this is necessary. Perhaps too many TERZ are written.

    Args:
        zmat_filename: String. Filename of zmat to be fixed
        ligand_resname: String. Residue name of the ligand
        original_pdb_filename: String. Filename of the original PDB.
    """

    print('Fixing Z matrix TERZ .... ', zmat_filename)

    tmpfile = open('tmp.txt','w')

    try: 
        be2allFile = [ele for ele in open(zmat_filename)]
    except:
        print("FAILED TO MAKE Z-MATRICES. PLEASE CHECK ABOVE ERROR MESSAGES. \n \
        (likely you have hetatoms that failed with BOSS) \n \
        NOTE THAT THE 'CHOPPED' PDB IS ALREADY PREPARED IF THAT IS ALL YOU DESIRE")
        sys.exit()
        manage_files(original_pdb_filename, ligand_resname)
        sys.exit()

    
    for iter in range(len(be2allFile)):
        new_line = be2allFile[iter]
        if 'TERZ' in new_line:
            if ligand_resname in be2allFile[iter+1] or 'CAP' in be2allFile[iter+1]:
                new_line = be2allFile[iter]
            else:
                continue

        tmpfile.write(new_line)

    tmpfile.close()
                

    os.system('cp '+zmat_filename+' lll.txt') 
    os.system('cp tmp.txt '+zmat_filename) 

def relax_protein_ligand_system(original_pdb_filename):
    """
    Run an optimization of the system using the xCGDD9 script from MCPRO.

    Args:
        original_pdb_filename: String. Filename of the original pdb complex
    """
    
    # complexPDBName = complexPDB[:-4]
    MCPRO_scripts_path = os.path.join(os.environ.get('MCPROdir'),'scripts')

    os.system('mkdir CG9;cp '+original_pdb_filename[:-4]+'.all.z CG9;cd CG9;'+MCPRO_scripts_path+'/xCGDD9 '+original_pdb_filename[:-4]+'.all')

def generate_final_structures_with_cap(original_pdb_filename):
    """
    Generates the final zmats.

    Args:
        original_pdb_filename: String. Filename of the original pdb complex
    """

    print('Generating final structures')

    optimized_zmat = get_optimized_zmat()


    replace_optimized_coordinates_in_zmat_file(original_pdb_filename[:-4]+'.var.z',
                                                original_pdb_filename[:-4]+'.cap.z',
                                                optimized_zmat)
    replace_optimized_coordinates_in_zmat_file(original_pdb_filename[:-4]+'.varcon.z',
                                                original_pdb_filename[:-4]+'.capcon.z',
                                                optimized_zmat)
    
    try:
        os.mkdir('finalZmatrices')
    except:
        os.system('rm -r finalZmatrices')
        os.mkdir('finalZmatrices')

    os.system('cp CG9/optsum finalZmatrices/'+original_pdb_filename[:-4]+'.all.z')   
    os.system('cp '+original_pdb_filename[:-4]+'.cap.z '+ original_pdb_filename[:-4]+'.capcon.z finalZmatrices')

def get_optimized_zmat():
    """
    Gets all of the info following the 'Geometry Variations follow' part of optsum

    """

    # By default the output file of xCGDD9 is called optsum 

    optimized_zmat = []

    for line in open('CG9/optsum'):

        if 'Geometry Variations follow' in line: break

        optimized_zmat.append(line)

    return optimized_zmat

def replace_optimized_coordinates_in_zmat_file(var_filename, cap_filename, optimized_zmat):
    """
    Args:
        var_filename: String. Filename of var zmat. e.g. chopped_pdb_filename[:-4]+'var.z
        cap_filename: String. Filename of the capped zmat. e.g. chopped_pdb_filename[:-4]+'capcon.z'
        optimized_zmat: List. List of lines in the optimized zmat
    """

    output_file = open(cap_filename,'w')

    # get bottom part of the var file

    read = False
    data_at_bottom_of_var_file = []
    
    for line in open(var_filename):
        if 'Geometry Variations follow' in line: read = True
        if read: data_at_bottom_of_var_file.append(line)

    for ele in optimized_zmat: output_file.write(ele)
    for ele in data_at_bottom_of_var_file: output_file.write(ele)

def add_protein_bonds(original_pdb_filename):
    """
    ### based on the BONDADD.f script by Yue and Jonah ### Basically,
    we are adding the additional bonds to the zmat. 

    Args:
        original_pdb_filename: String. Filename of the chopped pdb complex
    """



    # read the all zmat in as optzmat
    with open('finalZmatrices/'+original_pdb_filename[:-4]+'.all.z') as optzmat:
        optzmat_data = optzmat.readlines()
    with open('finalZmatrices/'+original_pdb_filename[:-4]+'.capcon.z') as varzmat:
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

    output_filename = 'finalZmatrices/'+original_pdb_filename[:-4]+'_final_capcon_zmat.z'
    with open(output_filename,'w') as output_file:
        for line in varzmat_data:
            output_file.write(line)

# def fixExcludedAtomsList(complexPDB):

#   complexPDBName = complexPDB[:-4]
#   filename = 'finalZmatrices/'+complexPDBName+'_final_capcon_zmat.z'

#   with open(filename) as zmat:
#         zmat_data = zmat.readlines()

#     index = 0
#     found = False
#     for line in zmat_data:
#       if 'Excluded Atoms List follows' in line:
#           found = True
#           break
#       index += 1

#     if found:
#       if '-1' not in zmat_data[index+1][0:3]:

def manage_files(original_pdb_filename, ligand_resname):
    """
    Sorts all of the produced files into organized directories.

    Args:
        original_pdb_filename: String. Filename of the original pdb complex
        ligand_resname: String. Residue name of the ligand.

    """


    pdb_id = pdb_id = os.path.splitext(os.path.basename(original_pdb_filename))[0] # [:-6] gets rid of '_fixed', which is used in zmat section
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
        os.makedirs(pdb_id + '_other_output')
    except:
        pass
    other_files_l = ['sum','log','out','lll.txt','plt.pdb','chop_script.csh','tmp.txt']
    for file in other_files_l:
        try:
            os.system('mv ' + file + ' ./' + pdb_id + '_other_output/')
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

    # run the main function
    main()


