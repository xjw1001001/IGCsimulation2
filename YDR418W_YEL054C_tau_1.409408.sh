#!/bin/bash
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_tau_1.409408_IGCgeo_3.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_tau_1.409408_IGCgeo_10.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_tau_1.409408_IGCgeo_50.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_tau_1.409408_IGCgeo_100.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/YDR418W_YEL054C_tau_1.409408_IGCgeo_500.0.sh 
