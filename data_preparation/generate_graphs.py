# -*- coding: utf-8 -*-
"""
@Time:     Created on 2020/10/12 14:49
@author:   Tianbiao Yang & Mingyue Zheng
@Email:    Tianbiao_Yang@163.com
@Filename: generate_graphs.py
@Software: Spyder & Python
"""

import os
import numpy as np
from collections import defaultdict
from scipy.spatial.distance import cosine, euclidean
from scipy.stats import percentileofscore as perc

# =============================================================================
# Aims:     Parse PDB Pocket files
# README:   This script parses pdb for residue types and (x,y,z) coordinates 
#           for construction of protein pocket graphs.
# =============================================================================

class GetPocketList():
    
    def __init__(self, dpath):
        self.dpath = dpath
        
    def PocketsList(self,pdbs_dict,bsmset_path):
        """Get the PDBs list."""
        with open(self.dpath + bsmset_path) as rpklf:
            for i_data in rpklf:
                i = i_data.rstrip('\n').split('\t')
                if i[0] == 'PDB_A':
                    pass
                else:
                    pdbs_dict[i[0]] = ' '
                    pdbs_dict[i[1]] = ' '
        return pdbs_dict
    
    def GetSitePath(self):
        bsmset_path = 'Datasets/BSMset.txt'
        vertex_path = 'Datasets/Vertex.txt'
        barelier_path = 'Datasets/Barelier.txt'
        set_paths = [bsmset_path,vertex_path,barelier_path] 
        pdbs_dict = dict()
        for p in set_paths:
            pdbs_dict = self.PocketsList(pdbs_dict,p)
        head = self.dpath + 'Structure/'
        site_path_list = [head + i.split('_')[0] + '/' + i + '_site.pdb' for i in list(pdbs_dict.keys())]
        return site_path_list
    
    def ParsePocket(self,site_path):
        """Parse PDB Pocket file information.

        Parameters
        ----------
        site_path: str
        return :   Index, Type, Position, Sidechain
        """
        # Parse residue, atom type and atomic coordinates
        residues = ["ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX", "GLY",
                    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                    "TYR", "UNK", "VAL",]
        seq_data, protein_data, sidechain_data = [], [], []
        complex_data = defaultdict(list)
        res_, res_i, res_c,chain_id = None, None, None, None
        sidechain_flag, sidechain_counter = False, 0

        with open(site_path,'r') as rpklf:
            lines = rpklf.readlines()
            for row in lines:
                if row[:4] == "ATOM":
                    # Check if for chain
                    chain_id = row[21].upper()
                    if res_i is None:
                        # residue number in the sequence
                        res_i = row[22:26]
                    if row[22:26] == res_i:
                        if row[12:17] in [" CA  ", " CA A"]:
                            # carbon alpha of the atom
                            res_ = row[17:20]
                            seq_data.append([row[21],residues.index(res_)])
                            # coordinates of the atom
                            res_c = [row[30:38].strip(),row[38:46].strip(),row[47:54].strip(),]
                            sidechain_flag = True
                            sidechain_counter += 1
                        elif row[12:17] in [" CB  ", " CB A"]:
                            # carbon alpha of the atom
                            res_b_ = row[17:20]
                            # coordinates of the atom
                            res_b_c = [row[30:38].strip(),row[38:46].strip(),row[47:54].strip(),]
                        else:
                            if sidechain_flag:
                                # a carbon alpha is asociated to a side chain. other atoms.
                                if sidechain_counter > 2:
                                    sidechain_data.append([row[30:38].strip(),
                                                           row[38:46].strip(),row[47:54].strip(),])
                                else:
                                    sidechain_counter += 1
                    else: 
                        try:
                            ress = residues.index(res_)
                        except:
                            ress = residues.index("UNK")
                        if len(sidechain_data) > 0:
                            sidechain_data = np.array(sidechain_data).astype("float")
                            sidechain_c = np.mean(sidechain_data, axis=0).tolist()
                            sidechain_data = []
                        else:
                            sidechain_c = res_c
                        sidechain_flag = False
                        sidechain_counter = 0
                        if res_c is not None:
                            res_data = [res_i, ress] + res_c + sidechain_c + res_b_c
                            protein_data.append(res_data)
                            complex_data[chain_id].append(res_data)                        
                        res_i = row[22:26]                

                if row[:3] == "TER":  # The last residues selected by the way
                    try:
                        ress = residues.index(res_)
                    except:
                        ress = residues.index("UNK")
                    if len(sidechain_data) > 0:
                        sidechain_data = np.array(sidechain_data).astype("float")
                        sidechain_c = np.mean(sidechain_data, axis=0).tolist()
                        sidechain_data = []
                    else:
                        sidechain_c = res_c
                    sidechain_flag = False
                    sidechain_counter = 0
                    if res_c is not None:
                        res_data = [res_i, ress] + res_c + sidechain_c + res_b_c
                        protein_data.append(res_data)
                        complex_data[chain_id].append(res_data)
                    res_i = row[22:26]

        data = {}
        for ii in complex_data.keys():
            chain_data = np.array(complex_data[ii])
            chain_c = chain_data[:, 2:5].astype("float")
            chain_sc_c = chain_data[:, 5:-3].astype("float")
            chain_centroid = np.mean(chain_c, axis=0)
            residue_depth = np.array([euclidean(chain_centroid, c) for c in chain_c])
            residue_depth_percentile = [1 - perc(residue_depth, d) / 100.0 for d in residue_depth]
            chain_c = chain_c - chain_centroid
            chain_sc_c = chain_sc_c - chain_centroid
            chain_sc_c = chain_sc_c - chain_c
            chain_c = -(chain_c)
            residue_orientation = list(np.nan_to_num([1 - cosine(chain_c[i],chain_sc_c[i]) 
                                        for i in range(len(chain_c))]))
            data[ii] =  [complex_data[ii][n] + [residue_depth_percentile[n]] + [residue_orientation[n]] 
                                        for n in range(0,len(complex_data[ii]))]    

        return data
    
    def SaveGraphs(self,site_path,site_data):
        site_name = site_path.split('/')[-1].split('_site.pdb')[0]
        with open(self.dpath + 'Graphs/' + site_name + '.txt','w') as wpklf:
            for k in site_data.keys():
                for res_value in site_data[k]:
                    res_value = [ str(round(float(i),6)) for i in res_value]
                    wpklf.write('\t'.join(res_value) + '\n')
        return site_name + ' is finished!'
    
if __name__ == "__main__":
    
    dpath = '/home/tbyang/Desktop/GraphBSM202001/data/'
    GetPockets = GetPocketList(dpath)
    site_path_list = GetPockets.GetSitePath()
    n,m = 0,0
    for site_path in site_path_list:
        try:
            site_data = GetPockets.ParsePocket(site_path)
            output = GetPockets.SaveGraphs(site_path,site_data)
            n = n + 1
            print(n,output)
        except:
            n = n + 1
            m = m + 1
            print(n,m,site_path,'Error!!!')