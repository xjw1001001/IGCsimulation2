from IGCCodonSimTree import *
import os
#-*- coding:utf-8 -*-

# This file is used to generate simulation data
# Xiang Ji
# xji3@ncsu.edu

# again, I am too lazy and I know it...

if __name__ == '__main__':
    '''
    # constants in the simulation
    # mostly from  inference with IGC expansion see MG94_YDR418W_YEL054C_nonclock_summary.txt
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'
    paralog = [paralog1, paralog2]
    newicktree = './YeastTree.newick'
    num_exon = 163
    tau_list = [0.0, 1.0, 1.409408, 10.0, 20.0]
    IGC_threshold = -0.1

    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    #IGC_geo_list = [3.0]
    for IGC_geo in IGC_geo_list:
        for tau in tau_list:
            IGC_geo_codon = IGC_geo / 3.0
            IGC_init = tau / IGC_geo_codon
            x_IGC = [IGC_init, 1.0 / IGC_geo_codon, IGC_threshold]  # These values vary for the simulation study
            
            #sim_num = 1
            for sim_num in range(30):
                log_folder = './' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/log/'
                div_folder = './' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/div/'
    
                # Now create folder if they don't exist
                
                if not os.path.isdir('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/'):
                    os.makedirs('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/')
    
                if not os.path.isdir('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'):
                    os.makedirs('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/')
    
                if not os.path.isdir(log_folder):
                    os.makedirs(log_folder)
    
                if not os.path.isdir(div_folder):
                    os.makedirs(div_folder)
    
                
                seq_file = './' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
                log_file = './' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.log'
    
                save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
    
                #x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
                pi_a = 2.983653651297249465e-01
                pi_c = 2.095729139806868369e-01
                pi_g = 1.871536409219771990e-01
                pi_t = 3.049080799676109899e-01
                kappa = 8.404333642032199236
                #omega = 7.618529705823594289e-02
                omega = 1.0
                x_exon = [pi_a + pi_g, pi_a / (pi_a + pi_g), pi_c / (pi_c + pi_t), kappa, omega]  # parameter from MG94_YDR418W_YEL054C_nonclock_summary.txt
    
                save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
                x_rates = np.exp([-3.925914311698581738, -1.533949335440421891, -1.564219363086546410, -1.762095498829083562, -3.660826292734839171, -3.626559649929371076,
                           -3.438639957077882059, -2.460890888502840657, -3.690966112631306473, -2.870638256407853639, -2.844813362532052192, -3.822236386160084542])
                seed =    sim_num + 10000     
        
                test = TreeIGCCodonSimulator(num_exon, newicktree, paralog, seq_file, log_file, x_exon, x_IGC, log_folder, div_folder, seed)
                test.unpack_x_rates(x_rates)
                try:
                    test.sim()
                except:
                    print 'failed at sim  ' + str(sim_num) + '  IGC_geo = ' + str(IGC_geo)
                    test.write_log()
     '''   
    paralog1 = 'EDN'
    paralog2 = 'ECP'
    paralog = [paralog1, paralog2]
    newicktree = './primate_EDN_ECP.newick'
    num_exon = 156
    tau_list = [0.0, 1.0, 0.4079238, 10.0, 20.0]
    IGC_threshold = -0.1

    IGC_geo_list = [3.0, 10.0, 50.0, 100.0, 500.0]
    #IGC_geo_list = [3.0]
    for IGC_geo in IGC_geo_list:
        for tau in tau_list:
            IGC_geo_codon = IGC_geo / 3.0
            IGC_init = tau / IGC_geo_codon
            x_IGC = [IGC_init, 1.0 / IGC_geo_codon, IGC_threshold]  # These values vary for the simulation study
            
            #sim_num = 1
            for sim_num in range(30):
                log_folder = './' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/log/'
                div_folder = './' + '_'.join(paralog) + '/tau'+ str(tau) +'/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/div/'
    
                # Now create folder if they don't exist
                
                if not os.path.isdir('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/'):
                    os.makedirs('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/')
    
                if not os.path.isdir('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/'):
                    os.makedirs('./' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/')
    
                if not os.path.isdir(log_folder):
                    os.makedirs(log_folder)
    
                if not os.path.isdir(div_folder):
                    os.makedirs(div_folder)
    
                
                seq_file = './' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.fasta'
                log_file = './' + '_'.join(paralog) + '/tau'+ str(tau) + '/IGCgeo_' + str(IGC_geo) + '/sim_' + str(sim_num) + '/' + '_'.join(paralog) + '_MG94_geo_' + str(IGC_geo) + '_Sim_' + str(sim_num) + '.log'
    
                save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
    
                #x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
                [pi_a,pi_c,pi_g,pi_t,kappa,omega]=[0.29606234407101895,
                 0.23774514722484202,
                 0.21594679926409824,
                 0.2502457094400408,
                 2.1252169762043467,
                 0.9167230032466169]
                x_exon = [pi_a + pi_g, pi_a / (pi_a + pi_g), pi_c / (pi_c + pi_t), kappa, omega]  # parameter from MG94_YDR418W_YEL054C_nonclock_summary.txt

                #cc = {akey[i]:avalue[i] for i in range(len(akey))}   
                # bb = [cc[i] for i in b]
                save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
                x_rates = np.array([0.1716208341560527, 0.1505609697476033, 0.04570771310319821, 
                                    0.017205023600801027, 0.010397410526678958, 0.036796522016806235, 
                                    0.007923586265531414, 0.021274867086765358, 2.593023836040686e-07, 
                                    0.03862769030856453, 2.8104293720398485e-07, 0.34165322297554845,
                                    0.009150267085027971, 0.011006437128204512, 0.008104361160863311,
                                    0.006661357109308685, 0.03767757792720381, 0.007085919094976171,
                                    0.035573130437498574, 0.1444062323563613, 0.09003798531799379,
                                    0.01699842377528958, 0.008462328496304105, 0.014185585255067459])
                seed =    sim_num + 10000     
        
                test = TreeIGCCodonSimulator(num_exon, newicktree, paralog, seq_file, log_file, x_exon, x_IGC, log_folder, div_folder, seed)
                test.unpack_x_rates(x_rates)
                try:
                    test.sim()
                except:
                    print 'failed at sim  ' + str(sim_num) + '  IGC_geo = ' + str(IGC_geo)
                    test.write_log()

