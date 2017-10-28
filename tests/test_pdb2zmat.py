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

test_list = [['1zu0.pdb','CBS'],
                ['2am1.pdb','1LG'],
                    ['5hvs.pdb','65V'],
                        ['RT_1.pdb','L11'],
                            ['rt_2.pdb','LIG']]

for l in test_list:
    pdb = l[0]
    lig = l[1]

    os.system('mv pdb_files/' + pdb + ' ./')

    try:
        os.system('python ../choptools/chop_protein.py -p ' + pdb + ' -r ' + lig + ' -c 18.0')
        print('SUCCEEDED ON ' + pdb)
    except:
        print('FAILED ON' + pdb)

    # try:
    #   os.system('rm -rf ./finalZmatricies')
    #   os.system('rm -rf ./CG9')
    #   #os.system('rm -rf)
    # except:
    #   pass






