# By Justin Meskas

import pandas as pd
import numpy as np
import os
from scipy.ndimage import shift, correlate
from scipy.optimize import minimize
from scipy.signal import correlate2d
#from scipy.ndimage import fourier_shift
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
import multiprocessing
import time

# normalize a matrix or vector
def normalize(array, min_val=None, max_val=None):

	if min_val is None:
		min_val = np.min(array)
	if max_val is None:
		max_val = np.max(array)

	if min_val == max_val:
		# Handle case where range is zero to avoid division by zero
		normalized_array = np.zeros_like(array)
	else:
		normalized_array = (array - min_val) / (max_val - min_val)

	return normalized_array

# first capping the density of the data to emphasize the low density populations and then calculate the correlations/convolutions of the two matrices to find the optimal shift.
def align_arrays(array1, array2, capping_value, bin_size):
	normalized_array1 = normalize(array1)
	normalized_array2 = normalize(array2)
			
	# capping_value = 1 is no capping
	mask = normalized_array1 >= capping_value
	normalized_array1[mask] = capping_value
	mask = normalized_array2 >= capping_value
	normalized_array2[mask] = capping_value
	
	correlation = correlate2d(normalized_array1, normalized_array2, mode='full')

	# force correlation only to have values that are bin_size/3 away from the centre. This avoids very large moves. TODO make code faster by only going over parts that are kept.
	zero_correlation = np.zeros_like(correlation)
	middle_rows_indices = slice(bin_size-round(bin_size/3), bin_size+round(bin_size/3))  # 43rd to 84th row indices (inclusive)
	middle_cols_indices = slice(bin_size-round(bin_size/3), bin_size+round(bin_size/3))  # 43rd to 84th column indices (inclusive)
	zero_correlation[middle_rows_indices, middle_cols_indices] = correlation[middle_rows_indices, middle_cols_indices]
	correlation = zero_correlation
		
	correlation_df = pd.DataFrame(correlation)
			
	max_index = np.unravel_index(np.argmax(correlation), correlation.shape)
	
	optimal_shift = (max_index[1]+1 - round(np.ceil(correlation.shape[1] / 2)), max_index[0]+1 - round(np.ceil(correlation.shape[0] / 2)))

	# Apply the optimal shift to array2
	aligned_array2 = np.roll(array2, optimal_shift[0], axis=1)
	aligned_array2 = np.roll(aligned_array2, optimal_shift[1], axis=0)

	return aligned_array2, optimal_shift

# Calculate how similar two matrices are.
def calculate_similariry_score(matrix1, matrix2):


	min_matrix = np.minimum(matrix1, matrix2)
	max_matrix = np.maximum(matrix1, matrix2)
		
	similarity_score = 1 - np.sum(max_matrix-min_matrix) / np.sum(max_matrix)

	return similarity_score

# help function for plotting. Choose between fast and slow implementations.
def scatter_setup(x3, y3, Fast):

	if(not Fast):	# Old method
		xy = np.vstack([x3,y3])
		z = gaussian_kde(xy)(xy) #,bw_method=1	

	if(Fast):	# New quicker method
		xy = np.vstack([x3,y3])
		kde = gaussian_kde(xy)

		# Regular grid to evaluate kde upon
		x3_flat = np.r_[min(x3):max(x3):64j]
		y3_flat = np.r_[min(y3):max(y3):64j]
		xs,ys = np.meshgrid(x3_flat,y3_flat)
		grid_coords = np.append(xs.reshape(-1,1),ys.reshape(-1,1),axis=1)

		hist = kde(grid_coords.T)
		hist = hist.reshape(64,64)

		x3_cut = np.digitize(x3, x3_flat, right=True)
		y3_cut = np.digitize(y3, y3_flat, right=True)

		z = np.array([])
		for i in range(0, len(x3)):
			dens_i = hist[y3_cut[i]-1, x3_cut[i]-1]
			z = np.append(z, dens_i)

	idx = z.argsort()
	x, y, z = x3[idx], y3[idx], z[idx]
	return x, y, z

# plot the scatter plot
def plot_scatter(x, y, z, xedges, yedges, title, data, cmap, limits=None):
	plt.scatter(x, y, c=z, s=1, cmap=cmap)
	
	if limits is None:
		plt.xlim(xedges[0], xedges[-1])
		plt.ylim(yedges[0], yedges[-1])
	else: 
		plt.xlim(limits[0], limits[2])
		plt.ylim(limits[3], limits[1])
	plt.title(f'{title}')
	plt.xlabel(f'{data.columns[0]}')
	plt.ylabel(f'{data.columns[1]}')
	
# remove events on the left, top, right and bottom according to percentile	
def remove_outliers(x, y, data, percentile=0.05):

	left_percentile = np.quantile(x, percentile)
	right_percentile = np.quantile(x, 1-percentile)
	bottom_percentile = np.quantile(y, percentile)
	top_percentile = np.quantile(y, 1-percentile)
	
	ind_left = np.where(x <= left_percentile)[0]
	ind_right = np.where(x >= right_percentile)[0]
	ind_bottom = np.where(y <= bottom_percentile)[0]
	ind_top = np.where(y >= top_percentile)[0]
	
	ind_all = np.concatenate((ind_left, ind_right, ind_bottom, ind_top))
	ind_unique = np.array(list(set(ind_all)))

	perc_removed = round(100*len(ind_unique)/data.shape[0],2)	
	
	x = x[~np.isin(np.arange(data.shape[0]), ind_unique)]
	y = y[~np.isin(np.arange(data.shape[0]), ind_unique)]
	x = x.reset_index(drop=True)
	y = y.reset_index(drop=True)

	return x, y, perc_removed

# Calculate the alignment between each pairwise files and save the shift
def alignment_2D_flowdata(file1, file2, j, k, project, Fast=False, bin_size=64, verbose=True, save=False, sample_size=20000, cutoff_percentile=0.002):

	filename1 = os.path.basename(file1)  # Extracts the filename including the extension
	filename1_wout_ext = os.path.splitext(filename1)[0]  # Removes the extension
	filename2 = os.path.basename(file2)  # Extracts the filename including the extension
	filename2_wout_ext = os.path.splitext(filename2)[0]  # Removes the extension

	if(verbose):
		print(f"Starting alignment analysis of {filename2_wout_ext} being aligned to {filename1_wout_ext}")

	plt.figure(figsize=(40, 16))
	matplotlib.rcParams['axes.labelsize'] = 18
	matplotlib.rcParams['xtick.labelsize'] = 18
	matplotlib.rcParams['ytick.labelsize'] = 18
	matplotlib.rcParams['font.size'] = 12
	plt.rcParams['figure.max_open_warning'] = 50

	cmap = plt.cm.jet
	cmaplist = [cmap(i) for i in range(cmap.N)]
	cmaplist = cmaplist[:230]
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

	data1 = pd.read_csv(f"{file1}",index_col=False)
	if (sample_size >= data1.shape[0]):
		sample_size = data1.shape[0]
	np.random.seed(42)
	sampled_rows = np.random.choice(range(data1.shape[0]), size=sample_size, replace=False)
	data1 = data1.loc[sampled_rows]
	data1.reset_index(drop=True, inplace=True)
	x1 = data1.iloc[:,0]
	y1 = data1.iloc[:,1]

	h1, x1edges, y1edges = np.histogram2d(x1, y1, bins=bin_size)	
	if (save==True):
		plt.subplot(2, 5, 1)
		x, y, z = scatter_setup(x1, y1, Fast=Fast)
		plot_scatter(x, y, z, x1edges, y1edges, title=filename1_wout_ext, data=data1, cmap=cmap)
		
	x1_org = x1
	y1_org = y1
	# skipping removing outliers for the master map
	#x1, y1, perc_removed_1 = remove_outliers(x1_org, y1_org, data1, cutoff_percentile)
	perc_removed_1 = 0
	h1, x1edges, y1edges = np.histogram2d(x1, y1, bins=bin_size)
	array1 = h1.T
		
	data2 = pd.read_csv(f"{file2}",index_col=False)
	if (sample_size >= data2.shape[0]):
		sample_size = data2.shape[0]
	np.random.seed(42)
	sampled_rows = np.random.choice(range(data2.shape[0]), size=sample_size, replace=False)
	data2 = data2.loc[sampled_rows]
	data2.reset_index(drop=True, inplace=True)
	x2 = data2.iloc[:,0]
	y2 = data2.iloc[:,1]
	
	h2, x2edges, y2edges = np.histogram2d(x2, y2, bins=bin_size)
	if (save==True):
		plt.subplot(2, 5, 6)
		x, y, z = scatter_setup(x2, y2, Fast=Fast)
		plot_scatter(x, y, z, x2edges, y2edges, title=filename2_wout_ext, data=data2, cmap=cmap)
	
	x2_org = x2
	y2_org = y2
	x2, y2, perc_removed_2 = remove_outliers(x2_org, y2_org, data2, cutoff_percentile) # 0.002 is best for subset_Longevity_Memory__CD3__CD4_CD8
	h2, x2edges, y2edges = np.histogram2d(x2, y2, bins=bin_size)
	array2 = h2.T
	
	if (save==True):
		plt.subplot(2, 5, 2)
		x, y, z = scatter_setup(x1, y1, Fast=Fast)
		plot_scatter(x, y, z, x1edges, y1edges, title=filename1_wout_ext, data=data1, cmap=cmap)

		plt.subplot(2, 5, 7)
		x, y, z = scatter_setup(x2, y2, Fast=Fast)
		plot_scatter(x, y, z, x2edges, y2edges, title=filename2_wout_ext, data=data2, cmap=cmap)
		
	# Align the arrays
	capping_value = 0.2
	aligned_array2, optimal_shift = align_arrays(array1, array2, capping_value, bin_size)
	
	if(verbose):
		print(f"Optimal Shift: {optimal_shift} of {filename2_wout_ext} being aligned to {filename1_wout_ext}")

	if (save==True):
		plt.subplot(2, 5, 3)
		plt.imshow(np.flipud(array1), cmap='jet', interpolation='nearest')

		plt.subplot(2, 5, 8)
		plt.imshow(np.flipud(array2), cmap='jet', interpolation='nearest')

	normalized_array1 = normalize(array1)
	normalized_array2 = normalize(array2)
	
	# capping_value = 1 is no capping
	mask = normalized_array1 >= capping_value
	normalized_array1[mask] = capping_value
	mask = normalized_array2 >= capping_value
	normalized_array2[mask] = capping_value

	if (save==True):
		plt.subplot(2, 5, 4)
		plt.imshow(np.flipud(normalized_array1), cmap='jet', interpolation='nearest')

		plt.subplot(2, 5, 9)
		plt.imshow(np.flipud(normalized_array2), cmap='jet', interpolation='nearest')

	x3 = x2 + optimal_shift[0]*(x2edges[-1]-(x2edges[0]))/bin_size
	y3 = y2 + optimal_shift[1]*(y2edges[-1]-(y2edges[0]))/bin_size

	if (save==True):
		plt.subplot(2, 5, 10)
		x, y, z = scatter_setup(x3, y3, Fast=Fast)
		plot_scatter(x, y, z, x2edges, y2edges, title="Aligned", data=data2, cmap=cmap)

	similarity_score = calculate_similariry_score(array1, array2)
	similarity_score_align = calculate_similariry_score(array1, aligned_array2)
	
	shift_too_large = False
	shift_allowed = 0.33
	if(abs(optimal_shift[0]) > shift_allowed*bin_size):
		shift_too_large = True
	if(abs(optimal_shift[1]) > shift_allowed*bin_size):
		shift_too_large = True
	if(abs(np.sqrt(optimal_shift[0]**2+optimal_shift[1]**2)) > shift_allowed*np.sqrt(2)*bin_size):
		shift_too_large = True
					
	if (save==True):			
		plt.subplot(2, 5, 5)
		plt.text(0.5, 0.8, f"Number of Bins: {bin_size}", ha='center', va='center', fontsize=18, color='black')
		plt.text(0.5, 0.7, f"Optimal Shift: {optimal_shift}", ha='center', va='center', fontsize=18, color='black')
		plt.text(0.5, 0.6, f"Similarity Score: {round(similarity_score,3)}", ha='center', va='center', fontsize=18, color='black')
		plt.text(0.5, 0.5, f"Aligned Similarity Score: {round(similarity_score_align,3)}", ha='center', va='center', fontsize=18, color='black')
		plt.text(0.5, 0.4, f"Events removed from file 1: {perc_removed_1}%", ha='center', va='center', fontsize=18, color='black')
		plt.text(0.5, 0.3, f"Events removed from file 2: {perc_removed_2}%", ha='center', va='center', fontsize=18, color='black')
		plt.text(0.5, 0.2, f"Shift too large: {shift_too_large}", ha='center', va='center', fontsize=18, color='black')
		plt.axis('off')

		output_plots = 'plots_pairwise_comparison'
		plt.tight_layout()
		plt.savefig(f'../results/{project}/{output_plots}/{filename2_wout_ext}__to__{filename1_wout_ext}.png')  # Save the plot
		plt.close()
			

	data_results = np.concatenate((optimal_shift, np.array([bin_size]), np.array([x2edges[0]]), np.array([x2edges[-1]]), np.array([y2edges[0]]), np.array([y2edges[-1]])))
	output_optshift = "pairwise_shift"
	optimal_shift_df = pd.DataFrame(np.column_stack(data_results), columns=["x_shift","y_shift", "Bin_Size", "x_min", "x_max", "y_min", "y_max"])
	if (save==True):
		optimal_shift_df.to_csv(f'../results/{project}/{output_optshift}/from__{filename2_wout_ext}__to__{filename1_wout_ext}.csv', index=False)
			
	return optimal_shift_df 

