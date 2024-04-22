# Main function of alignment. 
# This script was written by Justin, and modified by Eve for web web service.
# Currently not used

import pandas as pd
import numpy as np
import os
import matplotlib
import multiprocessing
import Alignment_WebService_func
import time

start_time = time.time()

Fast = True
bin_size = 64
verbose = True
save = True # True will create a plot and save a csv
cutoff_percentile = 0.002 
sample_size = 50000

project = "ANNOTAPI"

directory_path = f'../data/{project}/test_data'
file_list = os.listdir(directory_path)
file_list.sort()
num_files = len(file_list)

os.makedirs(f"../results/{project}/plots_pairwise_comparison", exist_ok=True)
os.makedirs(f"../results/{project}/pairwise_shift", exist_ok=True)


results = Alignment_WebService_func.alignment_2D_flowdata(file1=f"{directory_path}/{file_list[1]}", 
														  file2=f"{directory_path}/{file_list[0]}", 
														  j=1, 
														  k=0, 
														  project=project,
														  Fast=Fast, 
														  bin_size=bin_size, 
														  verbose=verbose, 
														  save=save, sample_size=sample_size)
if verbose:
	print(results)
				

end_time = time.time()
elapsed_time = end_time - start_time
if (verbose):
	print(f"Time to complete alignment analysis of {project}: {round(elapsed_time, 3)} seconds")


