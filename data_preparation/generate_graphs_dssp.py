# -*- coding: utf-8 -*-
"""
@Time:     Created on 2021/01/11 14:49
@author:   Tianbiao Yang & Mingyue Zheng
@Email:    Tianbiao_Yang@163.com
@Filename: generate_graphs_dssp.py
@Software: Spyder & Python
@Aims:     Change the pockets graph via the one hot and dssp programer
"""

import os
import oddt
import warnings
import numpy as np
from collections import defaultdict
from scipy.spatial.distance import cosine, euclidean
from scipy.stats import percentileofscore as perc
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import numpy as np
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

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
    
    
class GetResiduesFeature():
    
    def __init__(self, dpath):
        self.dpath = dpath
        
    def one_of_k_encoding_unk(self,x, allowable_set):
        """Maps inputs not in the allowable set to the last element."""
        if x not in allowable_set:
            x = allowable_set[-1]
        return list(map(lambda s: x == s, allowable_set))
    
    def one_of_k_encoding(self,x, allowable_set):
        if x not in allowable_set:
            raise Exception("input {0} not in allowable set{1}:".format(
                x, allowable_set))
        return list(map(lambda s: x == s, allowable_set))

    def get_central_coo(self, r,):
        atoms=r.atoms
        all_coor=[]
        for a in atoms:
            c=a.coords
            all_coor.append(c)
        all_coor=np.array(all_coor)
        central_coo=np.mean(all_coor,0)  
        return central_coo

    def Oddt_fearure(self,pdb_path):

        res_list = ['ASN', 'PHE', 'PRO', 'ARG', 'THR', 'ALA', 'VAL', 'ASP', 'MET', 'SER', 
                    'HIS', 'GLY', 'ILE', 'CYS', 'GLN', 'TRP', 'GLU', 'LEU', 'TYR', 'LYS']
        oddt_f = list()
        protein = next(oddt.toolkit.readfile('pdb', pdb_path))
        protein.protein = True
        for r in protein.residues:
            res_chain = r.chain
            res_index = r.number
            res_name = r.name
            res_center = [str(v) for v in self.get_central_coo(r,)]
            res_onehot = self.one_of_k_encoding_unk(res_name,res_list)
            res_onehot = [np.array(key, dtype=np.float) for key in res_onehot]
            oddt_f.append([res_chain,res_index,res_name,res_onehot,] + res_center)

        return oddt_f
    
    def DSSP_feature(self,pdb_path):
        possible_ss = ['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']
        p = PDBParser()
        structure = p.get_structure('temp', pdb_path)
        model = structure[0]
        dssp = DSSP(model, pdb_path)
        output = {'ss': [],
                  'relative_ASA': [],
                  'phi': [],
                  'psi': []}
        for item in dssp:
            output['ss'].append(item[2])
            output['relative_ASA'].append(item[3])
            output['phi'].append(item[4])
            output['psi'].append(item[5])
        output['ss'] = [self.one_of_k_encoding(single, possible_ss) for single in output['ss']]
        return {key: np.array(output[key], dtype=np.float) for key in output}
    
    def FIXerPDB(self,pdb_path):
        fixer = PDBFixer(filename=pdb_path)
        # Without fixer.missingResidues = {}, fixer.addMissingAtoms() throw an exception
        # and if I call fixer.findMissingResidues() several terminal residues are added
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        error_pdb_path = dpath + 'Structure/temp/' + pdb_path.split('/')[-1]
        PDBFile.writeFile(fixer.topology, fixer.positions, open(error_pdb_path, 'w'))
        return error_pdb_path
    
    def Residues_feature(self,pdb_path):
        dssp_f = self.DSSP_feature(pdb_path)
        ottd_f = self.Oddt_fearure(pdb_path)
        if len(dssp_f['relative_ASA']) == len(ottd_f):
            output = np.concatenate([ottd_f, dssp_f['ss'],
                                     dssp_f['relative_ASA'].reshape(-1, 1),
                                     dssp_f['phi'].reshape(-1, 1),
                                     dssp_f['psi'].reshape(-1, 1)], axis=1)
        else:
            error_pdb_path = self.FIXerPDB(pdb_path)
            dssp_f = self.DSSP_feature(error_pdb_path)
            ottd_f = self.Oddt_fearure(error_pdb_path)
            output = np.concatenate([ottd_f, dssp_f['ss'],
                                     dssp_f['relative_ASA'].reshape(-1, 1),
                                     dssp_f['phi'].reshape(-1, 1),
                                     dssp_f['psi'].reshape(-1, 1)], axis=1)
        return output

    def ParsePocket(self,site_path_list):
        
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        pocket_feature = dict()
        n = 0 
        error = list()
        for site_path in  site_path_list:
            n = n  +1
            output = self.Residues_feature(site_path)
            pocket_feature[site_path] = output
            print(n,site_path,'was finished !!!')
            
        return pocket_feature

    def SaveGraphs(self,site_path_list,pocket_feature):
        
        for site_path in site_path_list:
            site_name = site_path.split('/')[-1].split('_site.pdb')[0]
            with open(dpath + 'Graphs/' + site_name + '.txt','w') as wpklf:
                for res_value in pocket_feature[site_path]:
                    res_value = [str(v) for v in res_value[3]] + [str(v) for v in res_value[7:]] +  list(res_value[4:7])
                    wpklf.write('\t'.join(res_value) + '\n')
                    
        return print('All were finished !!!')
    
    
    
if __name__ == "__main__":
    dpath = '/home/tbyang/Desktop/GraphBSM3/data/'
    GetPockets = GetPocketList(dpath)
    site_path_list = GetPockets.GetSitePath()
    GetRes = GetResiduesFeature(dpath)
    pocket_feature = GetRes.ParsePocket(site_path_list)
    GetRes.SaveGraphs(site_path_list,pocket_feature)
    