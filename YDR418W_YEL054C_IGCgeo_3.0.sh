#!/bin/bash
sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_IGCgeo_3.0_sim_0.sh 
sbatch -o IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_IGCgeo_3.0_sim_1.sh 
