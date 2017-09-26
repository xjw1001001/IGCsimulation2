import os

if __name__ == '__main__':
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'
    paralog = [paralog1, paralog2]
    tau_list = [0.0, 1.0, 1.409408, 10.0, 20.0]
    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    sim_num_list = range(30)

    sh_line = 'sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/'


    for tau in tau_list:
        IGC_geo_sh_file = './' + '_'.join(paralog) + '_tau_'+str(tau)+ '.sh'
        with open(IGC_geo_sh_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for IGC_geo in IGC_geo_list: 
                f.write(sh_line + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh \n')
                with open('./ShFiles/' + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    for num_sim in sim_num_list:
                        g.write('python RunSimulation.py' +' --P1 '+ str(paralog1) +' --P2 '+ str(paralog2) + ' --Tau '+ str(tau) +  ' --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')

    paralog1 = 'EDN'
    paralog2 = 'ECP'
    paralog = [paralog1, paralog2]
    tau_list = [0.0, 1.0, 0.4079238, 3.0,6.0,10.0, 20.0]
    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    sim_num_list = range(30)

    sh_line = 'sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/'


    for tau in tau_list:
        IGC_geo_sh_file = './' + '_'.join(paralog) + '_tau_'+str(tau)+ '.sh'
        with open(IGC_geo_sh_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for IGC_geo in IGC_geo_list: 
                f.write(sh_line + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh \n')
                with open('./ShFiles/' + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    for num_sim in sim_num_list:
                        g.write('python RunSimulation.py' +' --P1 '+ str(paralog1) +' --P2 '+ str(paralog2) + ' --Tau '+ str(tau) +  ' --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')
    paralog1 = 'ERa'
    paralog2 = 'ERb'
    paralog = [paralog1, paralog2]
    tau_list = [0.0, 1.0, 0.27788, 3.0,6.0, 10.0, 20.0]
    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    sim_num_list = range(30)

    sh_line = 'sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/'


    for tau in tau_list:
        IGC_geo_sh_file = './' + '_'.join(paralog) + '_tau_'+str(tau)+ '.sh'
        with open(IGC_geo_sh_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for IGC_geo in IGC_geo_list: 
                f.write(sh_line + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh \n')
                with open('./ShFiles/' + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    for num_sim in sim_num_list:
                        g.write('python RunSimulation.py' +' --P1 '+ str(paralog1) +' --P2 '+ str(paralog2) + ' --Tau '+ str(tau) +  ' --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')

    paralog1 = 'MR'
    paralog2 = 'GR'
    paralog = [paralog1, paralog2]
    tau_list = [0.0, 1.0, 0.1630137, 3.0,6.0, 10.0, 20.0]
    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    sim_num_list = range(30)

    sh_line = 'sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/'


    for tau in tau_list:
        IGC_geo_sh_file = './' + '_'.join(paralog) + '_tau_'+str(tau)+ '.sh'
        with open(IGC_geo_sh_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for IGC_geo in IGC_geo_list: 
                f.write(sh_line + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh \n')
                with open('./ShFiles/' + '_'.join(paralog) + '_tau_'+str(tau)+'_IGCgeo_' + str(IGC_geo) + '.sh', 'w+') as g:
                    g.write('#!/bin/bash' + '\n')
                    for num_sim in sim_num_list:
                        g.write('python RunSimulation.py' +' --P1 '+ str(paralog1) +' --P2 '+ str(paralog2) + ' --Tau '+ str(tau) +  ' --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')

##    sh_line = 'sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles_multiply3/'
##
##    for IGC_geo in IGC_geo_list:
##        IGC_geo_sh_file = './' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_multiply3.sh'
##        with open(IGC_geo_sh_file, 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for num_sim in sim_num_list:
##                f.write(sh_line + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh \n')
##                with open('./ShFiles_multiply3/' + '_'.join(paralog) + '_IGCgeo_' + str(IGC_geo) + '_sim_' + str(num_sim) + '_multiply3.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    g.write('python RunSimulation_multiply3.py --Geo ' + str(IGC_geo) + ' --sim_num ' + str(num_sim) + '\n')
