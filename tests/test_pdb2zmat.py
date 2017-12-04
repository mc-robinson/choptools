# !/usr/bin/env python
# File name: smiles_to_ligpargen.py
# Author: Matt Robinson
# Date created: 06/07/2017
# Date last modified: 06/07/2017
# Python Version: 3.6
"""
usage: python test_pdb2zmat.py

"""
import os

test_list = [['rt_2.pdb','LIG'],
            ['1zu0.pdb','CBS'],
            ['2am1.pdb','1LG'],
            ['5hvs_one_lig.pdb','65V'],
            ['1a9u.pdb','SB2'],
            ['1a28_chain_A.pdb','STR'],
            ['2dbl.pdb','S5H'],
            ['7tim.pdb','PGH'],
            ['1a6W.pdb','NIP']]

#test_list = [['MIF-180_cplx_orig.cm5.pdb','LIG']]
#test_list = [['RT_1.pdb','L11']]

for l in test_list:
    pdb = l[0]
    lig = l[1]

    #os.system('mv pdb_files/' + pdb + ' ./')
    os.system('mv pdb_files/' + pdb + ' ./')

    try:
        os.system('python ../choptools/pdb2zmat.py -p ' + pdb + ' -r ' + lig + ' -c 18.0 -t')
        print('SUCCEEDED ON ' + pdb)
    except:
        print('FAILED ON' + pdb)

    # os.system("mkdir "+pdb[0:-4]+'_folder')

    # list_of_files = os.listdir(os.getcwd())

    # # need to defined outside of loop, excluded_files = ['.DS_Store', '.ipynb_checkpoints', '__init__.py', 'test_pdb2zmat.py', 'pdb_files', 'further_tests']
    # excluded_files.append(pdb[0:-4]+'_folder')

    # for file in list_of_files:

    #     if file.startswith('original_'):
    #         os.system("mv "+file+' ./pdb_files/'+pdb)

    #     elif file not in excluded_files:
    #         os.system("mv "+file+' ./'+pdb[0:-4]+'_folder')

    # try:
    #   os.system('rm -rf ./finalZmatricies')
    #   os.system('rm -rf ./CG9')
    #   #os.system('rm -rf)
    # except:
    #   pass
