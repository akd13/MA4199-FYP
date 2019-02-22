#!/bin/bash


for i in {18..36}
do
	python convertBL_human_subprocess.py Data/CGATGT-s_6_1_seed29_genome."$i".bl 
	rm -f Data/*.track
done 