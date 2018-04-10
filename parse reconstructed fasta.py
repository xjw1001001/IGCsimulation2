# -*- coding: utf-8 -*-
"""List of extant and reconstructed sequences
Overall accuracy of the 24 ancestral sequences:
Created on Thu Aug 10 08:23:33 2017

@author: xjw1001001
"""'''
#only when PAML in desktop is available,the yeast version only
import os
def read_rst(filename,filename1,filename2):
    f = open(filename,"r")  
    f1 = open(filename1,"w+")  
    f2 = open(filename2,"w+")  
    lines = f.readlines()#读取全部内容

    tree1_flag1 = 0
    tree2_flag1 = 0
    tree1_flag2 = 0
    tree2_flag2 = 0
    tree1_flag3 = 0
    flag = 0
    lines_index = 0
    for line in lines:
        lines_index = lines_index + 1
        if flag == 0:
            if line == 'List of extant and reconstructed sequences\n':
                tree1_flag1 = lines_index-1#flag1+4
                flag = flag + 1
                continue
        if flag == 1:
            if line == 'Overall accuracy of the 23 ancestral sequences:\n':
                tree1_flag2 = lines_index-1#flag2-3
                continue
            if line == 'List of extant and reconstructed sequences\n' :
                flag = flag + 1
                continue
        if flag == 2:
            if line == 'List of extant and reconstructed sequences\n':
                tree2_flag1 = lines_index-1#flag1+4
                flag = flag + 1
                continue
        if flag == 3: 
            if line == 'Overall accuracy of the 23 ancestral sequences:\n' :
                tree2_flag2 = lines_index-1#flag2-3
                break 
    lines_index = 0                    
    for line in lines:  
        if len(line)>= 5:
            if lines_index >= tree1_flag1+4 and lines_index <= tree1_flag2-3:
                if len(line)>=10:
                    f1.write(line)
            if lines_index >= tree2_flag1+4 and lines_index <= tree2_flag2-3:
                if len(line)>=10:
                    f2.write(line)
        lines_index = lines_index + 1
    f.close()
    f1.close()
    f2.close()


pair = ['EDN','ECP']
tau_list = [0.0,0.1,  0.7, 1.0, 0.4079238,3.0,6.0, 10.0, 20.0]
IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]

path = '/Users/xjw1001001/Documents/GitHub/IGCsimulation2/PAMLresult/'+'_'.join(pair)+'/'

for tau in tau_list:
    for IGCgeo in IGC_geo_list:
        for sim in range(30):
            primalline=[]
            fastaline=[]
            rstpath = path   +'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) +'/'
            read_rst(rstpath + 'rst',rstpath + 'rst_tree1',rstpath + 'rst_tree2')
            with open(rstpath + 'rst_tree1','r') as f:
                for line in f.readlines():
                    primalline.append(line)
                    sline = '>' + line
                    sline=sline.replace('node #26','N0'+pair[0])
                    sline=sline.replace(' ','')
                    sline=sline.replace('\n','')
                    for i in range(11):
                        sline=sline.replace('node#' + str(26+1+i),'N'+str(1+i)+pair[0])
                        sline=sline.replace('node#' + str(37+1+i),'N'+str(1+i)+pair[1])
                    sline=sline.replace(pair[0],pair[0] + '\n')
                    sline=sline.replace(pair[1],pair[1] + '\n')
                    fastaline.append(sline)    
            f1 = open('/Users/xjw1001001/Documents/GitHub/IGCsimulation2/test/Ancestral_reconstruction/series/'+'_'.join(pair)+'/'
                      +'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) + '/ancestral_reconstruction_' + '_'.join(pair) + '_PAML.fasta','w+')
            for line in fastaline:
                f1.write(line)
                f1.write('\n')
            f1.close()


''''''
# -*- coding: utf-8 -*-
"""List of extant and reconstructed sequences
Overall accuracy of the 24 ancestral sequences:
Created on Thu Aug 10 08:23:33 2017

@author: xjw1001001
"""
#only when PAML in desktop is available,the yeast version only
import os
def read_rst(filename,filename1,filename2):
    f = open(filename,"r")  
    f1 = open(filename1,"w+")  
    f2 = open(filename2,"w+")  
    lines = f.readlines()#读取全部内容

    tree1_flag1 = 0
    tree2_flag1 = 0
    tree1_flag2 = 0
    tree2_flag2 = 0
    tree1_flag3 = 0
    flag = 0
    lines_index = 0
    for line in lines:
        lines_index = lines_index + 1
        if flag == 0:
            if line == 'List of extant and reconstructed sequences\n':
                tree1_flag1 = lines_index-1#flag1+4
                flag = flag + 1
                continue
        if flag == 1:
            if line == 'Overall accuracy of the 11 ancestral sequences:\n':
                tree1_flag2 = lines_index-1#flag2-3
                continue
            if line == 'List of extant and reconstructed sequences\n' :
                flag = flag + 1
                continue
        if flag == 2:
            if line == 'List of extant and reconstructed sequences\n':
                tree2_flag1 = lines_index-1#flag1+4
                flag = flag + 1
                continue
        if flag == 3: 
            if line == 'Overall accuracy of the 11 ancestral sequences:\n' :
                tree2_flag2 = lines_index-1#flag2-3
                break 
    lines_index = 0                    
    for line in lines:  
        if len(line)>= 5:
            if lines_index >= tree1_flag1+4 and lines_index <= tree1_flag2-3:
                if len(line)>=10:
                    f1.write(line)
            if lines_index >= tree2_flag1+4 and lines_index <= tree2_flag2-3:
                if len(line)>=10:
                    f2.write(line)
        lines_index = lines_index + 1
    f.close()
    f1.close()
    f2.close()


pair = ['YDR418W','YEL054C']
tau_list = [0.0,0.1, 0.3, 0.5, 0.7, 1.0, 1.409408, 10.0, 20.0]
IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]

path = '/Users/xjw1001001/Documents/GitHub/IGCsimulation2/PAMLresult/'+'_'.join(pair)+'/'

for tau in tau_list:
    for IGCgeo in IGC_geo_list:
        for sim in range(30):
            primalline=[]
            fastaline=[]
            rstpath = path   +'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) +'/'
            read_rst(rstpath + 'rst',rstpath + 'rst_tree1',rstpath + 'rst_tree2')
            with open(rstpath + 'rst_tree1','r') as f:
                for line in f.readlines():
                    primalline.append(line)
                    sline = '>' + line
                    sline=sline.replace('node #14','N0'+pair[0])
                    sline=sline.replace(' ','')
                    sline=sline.replace('\n','')
                    for i in range(5):
                        sline=sline.replace('node#' + str(14+1+i),'N'+str(1+i)+pair[0])
                        sline=sline.replace('node#' + str(19+1+i),'N'+str(1+i)+pair[1])
                    sline=sline.replace(pair[0],pair[0] + '\n')
                    sline=sline.replace(pair[1],pair[1] + '\n')
                    fastaline.append(sline)    
            f1 = open('/Users/xjw1001001/Documents/GitHub/IGCsimulation2/test/Ancestral_reconstruction/series/'+'_'.join(pair)+'/'
                      +'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) + '/ancestral_reconstruction_' + '_'.join(pair) + '_PAML.fasta','w+')
            for line in fastaline:
                f1.write(line)
                f1.write('\n')
            f1.close()'''
import os
def read_rst(filename,filename1,filename2):
    f = open(filename,"r")  
    f1 = open(filename1,"w+")  
    f2 = open(filename2,"w+")  
    lines = f.readlines()#读取全部内容

    tree1_flag1 = 0
    tree2_flag1 = 0
    tree1_flag2 = 0
    tree2_flag2 = 0
    tree1_flag3 = 0
    flag = 0
    lines_index = 0
    for line in lines:
        lines_index = lines_index + 1
        if flag == 0:
            if line == 'List of extant and reconstructed sequences\n':
                tree1_flag1 = lines_index-1#flag1+4
                flag = flag + 1
                continue
        if flag == 1:
            if line == 'Overall accuracy of the 19 ancestral sequences:\n':
                tree1_flag2 = lines_index-1#flag2-3
                continue
            if line == 'List of extant and reconstructed sequences\n' :
                flag = flag + 1
                continue
        if flag == 2:
            if line == 'List of extant and reconstructed sequences\n':
                tree2_flag1 = lines_index-1#flag1+4
                flag = flag + 1
                continue
        if flag == 3: 
            if line == 'Overall accuracy of the 19 ancestral sequences:\n' :
                tree2_flag2 = lines_index-1#flag2-3
                break 
    lines_index = 0                    
    for line in lines:  
        if len(line)>= 5:
            if lines_index >= tree1_flag1+4 and lines_index <= tree1_flag2-3:
                if len(line)>=10:
                    f1.write(line)
            if lines_index >= tree2_flag1+4 and lines_index <= tree2_flag2-3:
                if len(line)>=10:
                    f2.write(line)
        lines_index = lines_index + 1
    f.close()
    f1.close()
    f2.close()


pair = ['MR','GR']
tau_list = [0.0,0.5,  1.0, 0.1630137]
IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]

path = '/Users/xjw1001001/Documents/GitHub/IGCsimulation2/PAMLresult/'+'_'.join(pair)+'/'

for tau in tau_list:
    for IGCgeo in IGC_geo_list:
        for sim in range(30):
            primalline=[]
            fastaline=[]
            rstpath = path   +'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) +'/'
            read_rst(rstpath + 'rst',rstpath + 'rst_tree1',rstpath + 'rst_tree2')
            with open(rstpath + 'rst_tree1','r') as f:
                for line in f.readlines():
                    primalline.append(line)
                    sline = '>' + line
                    sline=sline.replace('node #22','N0'+pair[0])
                    sline=sline.replace(' ','')
                    sline=sline.replace('\n','')
                    for i in range(9):
                        sline=sline.replace('node#' + str(22+1+i),'N'+str(1+i)+pair[0])
                        sline=sline.replace('node#' + str(31+1+i),'N'+str(1+i)+pair[1])
                    sline=sline.replace(pair[0],pair[0] + '\n')
                    sline=sline.replace(pair[1],pair[1] + '\n')
                    fastaline.append(sline)    
            f1 = open('/Users/xjw1001001/Documents/GitHub/IGCsimulation2/test/Ancestral_reconstruction/series/'+'_'.join(pair)+'/'
                      +'tau' + str(tau) + '/IGCgeo_' + str(IGCgeo) + '/sim_' + str(sim) + '/ancestral_reconstruction_' + '_'.join(pair) + '_PAML.fasta','w+')
            for line in fastaline:
                f1.write(line)
                f1.write('\n')
            f1.close()