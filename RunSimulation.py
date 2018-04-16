#-*- coding:utf-8 -*-

from CodonGeneconv import ReCodonGeneconv
import os
import argparse
'''
def main(args):
    paralog1 = args.p1
    paralog2 = args.p2
    paralog = [paralog1, paralog2]
    newicktree = args.treefile

    tau     = args.tau
    IGC_geo = args.IGC_geo
    sim_num = args.sim_num

    alignment_file = './' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_leaf.fasta'
    save_name = './SimulationSave/' + '_'.join(paralog) +'/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/' + '_'.join(paralog) + '_MG94_geo_'  + str(IGC_geo) + '_Sim_' + str(sim_num) + '_save.txt'
    summary_path = './SimulationSummary/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
    summary_name1 = './SimulationSummary/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' +'IGCsummary.txt'
    summary_name2 = './SimulationSummary/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' +'tau=0summary.txt'
    # generate folder if not exist
    if not os.path.isdir('./SimulationSave/' + '_'.join(paralog) + '/'):
        os.makedirs('./SimulationSave/' + '_'.join(paralog) + '/')
    
    if not os.path.isdir('./SimulationSave/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/'):
        os.makedirs('./SimulationSave/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/')
    
    if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog) + '/'):
        os.makedirs('./SimulationSummary/' + '_'.join(paralog) + '/')
    
    if not os.path.isdir(summary_path):
        os.makedirs(summary_path)
        
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False, save_name = save_name,IGC_geo = IGC_geo, sim_num = sim_num,realtau = tau)
    
    
    #test.get_mle(False, True, 0, 'BFGS')
    test.site_reconstruction()
    test.get_individual_summary(summary_path = './SimulationSummary/' + '_'.join(paralog) + '/', file_name = summary_name1)
    test.save_likelihood(summary_path = summary_path)
    test.Expected_tau_for_sitewise_and_branchwise(summary_path = summary_path)
    #MG94_tau_series = MG94_tau.reconstruction_series
    #MG94_tau_likelihooddict = MG94_tau.likelihood_dict
    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = False, save_name = save_name,IGC_geo = IGC_geo, sim_num = sim_num,realtau = tau)
    test2.site_reconstruction()
    test2.save_likelihood(summary_path = summary_path)
    test2.get_individual_summary(summary_path = './SimulationSummary/' + '_'.join(paralog) + '/', file_name = summary_name2)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--P1', dest = 'p1', help = 'paralog1')
    parser.add_argument('--P2', dest = 'p2', help = 'paralog2')
    parser.add_argument('--treefile', dest = 'treefile', help = 'treefile')
    parser.add_argument('--Tau', dest = 'tau', help = 'IGC tau')
    parser.add_argument('--Geo', dest = 'IGC_geo', help = 'IGC tract length parameter')
    parser.add_argument('--sim_num', dest = 'sim_num', help = 'Simulation number')
    
    main(parser.parse_args())    

'''
from CodonGeneconv import ReCodonGeneconv
import os
import argparse


paralog = ['MR', 'GR']
newicktree = './Thornton_MRGRARPR.newick'
tau     = 0.0
IGC_geo = 3.0
sim_num = 0

alignment_file = './' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '_leaf.fasta'
save_name = './SimulationSave/' + '_'.join(paralog) +'/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/' + '_'.join(paralog) + '_MG94_geo_'  + str(IGC_geo) + '_Sim_' + str(sim_num) + '_save.txt'
summary_path = './SimulationSummary/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'
summary_name1 = './SimulationSummary/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' +'IGCsummary.txt'
summary_name2 = './SimulationSummary/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' +'tau=0summary.txt'
# generate folder if not exist
if not os.path.isdir('./SimulationSave/' + '_'.join(paralog) + '/'):
    os.makedirs('./SimulationSave/' + '_'.join(paralog) + '/')

if not os.path.isdir('./SimulationSave/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/'):
    os.makedirs('./SimulationSave/' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/')

if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog) + '/'):
    os.makedirs('./SimulationSummary/' + '_'.join(paralog) + '/')

if not os.path.isdir(summary_path):
    os.makedirs(summary_path)
    
test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False, save_name = save_name,IGC_geo = IGC_geo, sim_num = sim_num,realtau = tau)


#test.get_mle(False, True, 0, 'BFGS')
test.site_reconstruction()
test.get_individual_summary(summary_path = './SimulationSummary/' + '_'.join(paralog) + '/', file_name = summary_name1)
test.save_likelihood(summary_path = summary_path)
test.Expected_tau_for_sitewise_and_branchwise(summary_path = summary_path)
#MG94_tau_series = MG94_tau.reconstruction_series
#MG94_tau_likelihooddict = MG94_tau.likelihood_dict
test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = {5:0.0}, clock = False, save_name = save_name,IGC_geo = IGC_geo, sim_num = sim_num,realtau = tau)
test2.site_reconstruction()
test2.save_likelihood(summary_path = summary_path)
test2.get_individual_summary(summary_path = './SimulationSummary/' + '_'.join(paralog) + '/', file_name = summary_name2)
