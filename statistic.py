# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 09:42:38 2018

@author: xjw1001001
"""
import numpy as np
from IGCexpansion.CodonGeneconFunc import *

#树文件位置
EDNECP_newicktree ='/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/Zhang2002_data ECPEDN/from gene bank/primate_EDN_ECP.newick'
Yeast_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/YeastTree.newick'
ERa_ERb_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/SR_Thornton/ER/species.newick'
ARa_ERa_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/SR_Thornton/ARa_ERa/ERa_ARa_species.newick'
ARMRGRPR_newicktree = '/Users/xjw1001001/Documents/GitHub/IGCexpansion2/reconstruction_data/SR_Thornton/AR_MR_GR_PR/species_common/species_common.newick'

#四种simulation
paralogs_pair = [['MR', 'GR'],['ERa','ERb'],['EDN','ECP'],['YDR418W','YEL054C']]

#真序列位置
true_series_path = '/Users/xjw1001001/Documents/GitHub/IGCsimulation2/' #+_'.join(pair)+'/'

#重建序列位置
constructed_series_path = '/Users/xjw1001001/Documents/GitHub/IGCsimulation2/test/Ancestral_reconstruction/series/' #+_'.join(pair)+'/'

#序列选择文件
#for tau in tau_list:
#    for IGCgeo in IGCgeo_list:
#        for sim in range(30):
#             + 'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) + '/ancestral_reconstruction_' + _'.join(pair) + '_MG94.fasta'


                           
                           
                           
                           
                           