# -*- coding: utf-8 -*-
"""
@Time:     Created on 2020/11/26 14:49
@author:   Tianbiao Yang & Mingyue Zheng
@Email:    Tianbiao_Yang@163.com
@Filename: fetch_pdbs.py
@Software: Spyder & Python
"""

import argparse
import os
import numpy as np
import wget
from mpi4py import MPI

# =============================================================================
# Aims:     Load the PDBs files
# README:   The input files was BSMset, VERTEXset, BARELIERset
#           The output was the download pdb files and the correct and error pdbs
# =============================================================================

class GetPDBList():
    
    def __init__(self, dpath):
        self.dpath = dpath
        
    def PDBsList(self,pdbs_dict,bsmset_path):
        """Get the PDBs list."""
        with open(self.dpath + bsmset_path) as rpklf:
            for i_data in rpklf:
                i = i_data.rstrip('\n').split('\t')
                if i[0] == 'PDB_A':
                    pass
                else:
                    pdbs_dict[i[0].split('_')[0]] = ' '
                    pdbs_dict[i[1].split('_')[0]] = ' '
        return pdbs_dict
    
    def FetchPDBs(self):
        bsmset_path = 'BSMset/BSMset.txt'
        vertex_path = 'VERTEXset/Revised_Vertex_Updata02.txt'
        barelier_path = 'BARELIERset/Barelier_V1_Updata02.txt'
        set_paths = [bsmset_path,vertex_path,barelier_path] 
        pdbs_dict = dict()
        for p in set_paths:
            pdbs_dict = self.PDBsList(pdbs_dict,p)
            
        return pdbs_dict
    
    def SaveLogs(self,out_dict):
        with open('download_logs.txt','a+') as wpklf:
            for k,v in out_dict.items():
                wpklf.write(k + '\t' + v  +'\n')
        return 'Download was finished!'
    
    def LoadPDBs(self, pdbs_dict):
        """Fetch pdb file from Protein Data Bank repo and store it in `path`."""
        # Download PDB file from database.
        correct_pdbs,error_pdbs = dict(),dict()
        n = 0
        for pdb_id in pdbs_dict:
            n = n + 1
            print(n, pdb_id)
            pdbs_path = self.dpath + 'Structure/'
            url = "https://files.rcsb.org/download/"  # URL used to fetch PDB files
            file_path = pdbs_path + pdb_id.lower() + ".pdb"
            try:
                if not os.path.exists(file_path):
                    file_path = wget.download(url + pdb_id + ".pdb", out=pdbs_path, bar=None)
                correct_pdbs[pdb_id] = 'correct'
            except:
                error_pdbs[pdb_id] = 'error'

        return self.SaveLogs(correct_pdbs), self.SaveLogs(error_pdbs)
    
if __name__ == '__main__':
    dpath = '/home/tbyang/Desktop/GraphBSM/data/'
    DownloadPDBs = GetPDBList(dpath)
    pdbs_dict = list(DownloadPDBs.FetchPDBs().keys())
    correct_pdbs, error_pdbs = DownloadPDBs.LoadPDBs(pdbs_dict)
    