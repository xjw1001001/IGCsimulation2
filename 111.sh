#!/bin/bash
sbatch -p bigmem -o ./cluster_outs/IGCSim-%j.out --mail-type=ALL --mail-user=xjw1001001@qq.com ./ShFiles/generate.sh 