# Neural Array Processing and Analysis Systems
Analysis code for processing and analyzing intracortical neural data

## Module Documentation

NAPAS_spikesort : Module for organizing analog snike snippets into discrete units
- NAPAS_SHELL_spikesort : Shell script for initializing PCA-based spike sorting
	- NAPAS_spikesort_returnData : Main wrapper function to organize data and perform PCA-sorting
		- NAPAS_spikesort_runPCA : Uses PCA to run EM algorithm on each channel's snippets
	- NAPAS_spikesort_makeChanMask : Creates mask for useful units

