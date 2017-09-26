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
    tau_list = [0.0, 1.0, 0.4079238, 3.0,6.0,10.0, 20.0]
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
                    
    paralog1 = 'ERa'
    paralog2 = 'ERb'
    paralog = [paralog1, paralog2]
    newicktree = './ThorntonERaERb.newick'
    num_exon = 310
    tau_list = [0.0, 1.0, 0.27788, 3.0,6.0, 10.0, 20.0]
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
                [pi_a,pi_c,pi_g,pi_t,kappa,omega]=[0.2170714117031952,
                                                     0.29568916383769184,
                                                     0.28354385627921275,
                                                     0.20369556817990023,
                                                     1.9064363361739103,
                                                     0.27617358730655541]#0.07617358730655541
                x_exon = [pi_a + pi_g, pi_a / (pi_a + pi_g), pi_c / (pi_c + pi_t), kappa, omega]  # parameter from MG94_YDR418W_YEL054C_nonclock_summary.txt

                #cc = {akey[i]:avalue[i] for i in range(len(akey))}   
                # bb = [cc[i] for i in b]
                save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
                x_rates = np.array([0.8078797759350917, 0.7907469852126532, 0.6380568187603329, 0.16134985943814728, 0.25828192353500096, 0.18897733989396742, 0.0028000090411836455, 0.4756967056541251, 0.20288238122081748, 0.20638849479219384, 0.30802527677786434, 0.4390494944806729, 0.07218239329460224, 0.03032528609934204, 5.404106008091292, 0.463523405702731, 0.143696737956612, 0.3137954895063944, 0.06256580618298248, 0.016935030009915572, 0.4962560207474495, 0.4290286230158284, 0.3470923161901951, 0.6634417040326285, 0.05107702867371775, 0.018278866996082287, 0.08555949996275779, 0.12693773714550619, 0.1463793073824961, 0.31875999896962465])
                seed =    sim_num + 10000     
        
                test = TreeIGCCodonSimulator(num_exon, newicktree, paralog, seq_file, log_file, x_exon, x_IGC, log_folder, div_folder, seed)
                test.unpack_x_rates(x_rates)
                try:
                    test.sim()
                except:
                    print 'failed at sim  ' + str(sim_num) + '  IGC_geo = ' + str(IGC_geo)
                    test.write_log()

    paralog1 = 'MR'
    paralog2 = 'GR'
    paralog = [paralog1, paralog2]
    newicktree = './Thornton MRGRARPR.newick'
    num_exon = 342
    tau_list = [0.0, 1.0, 0.1630137, 3.0,6.0, 10.0, 20.0]
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
                [pi_a,pi_c,pi_g,pi_t,kappa,omega]=[0.24513021771794125,
                                                         0.26931311846342154,
                                                         0.27926356228265115,
                                                         0.20629310153598615,
                                                         1.737042856801874,
                                                         0.5911879834911865]#0.05911879834911865 为什么omega总很小
                x_exon = [pi_a + pi_g, pi_a / (pi_a + pi_g), pi_c / (pi_c + pi_t), kappa, omega]  # parameter from MG94_YDR418W_YEL054C_nonclock_summary.txt

                #cc = {akey[i]:avalue[i] for i in range(len(akey))}   
                # bb = [cc[i] for i in b]
                save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
                x_rates = np.array([1.533977222521731,0.6950365800577575, 0.6899231387787704, 0.2769601077973367,
                                 0.36638478535655744, 0.032961870325235794, 0.095529997423931, 0.011915635029788544, 0.21005606283675432,
                                 10.346138780358743, 0.7132296368611716, 0.6674228387142348, 0.528853410762456, 0.4447167241135865,
                                 0.15554041100814398, 0.13742806331482885, 0.09068147404529761, 0.08016059317911103, 0.05539654982103279, 0.031024603497229703])
                seed =    sim_num + 10000     
        
                test = TreeIGCCodonSimulator(num_exon, newicktree, paralog, seq_file, log_file, x_exon, x_IGC, log_folder, div_folder, seed)
                test.unpack_x_rates(x_rates)
                try:
                    test.sim()
                except:
                    print 'failed at sim  ' + str(sim_num) + '  IGC_geo = ' + str(IGC_geo)
                    test.write_log()