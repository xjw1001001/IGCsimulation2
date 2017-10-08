#!/bin/bash
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_0.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_0.4079238.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_0.1.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_0.7.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_1.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_3.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_6.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_10.0.sh 
sbatch -o ./cluster_outs/IGCSim-%j.out --mail-type=FAIL --mail-user=xjw1001001@qq.com ./ShFiles/PAMLEDNECP_20.0.sh 
