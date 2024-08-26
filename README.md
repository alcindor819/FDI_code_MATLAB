# Fragment Dispersity Index: A cfDNA fragmentation pattern precise describing chromatin accessibility and its application in early cancer detection

![AppVeyor](https://img.shields.io/badge/MATLAB2020a-red)


## Introduction
This repository contains the pre-training data, source code, and package for the paper "Fragment Dispersity Index: A cfDNA fragmentation pattern precise describing chromatin accessibility and its application in early cancer detection".
 **FDI** is a  cfDNA fragmentation pattern which precisely characterizes chromatin accessibility. Utilizing FDI, we identified key regulatory elements involved in cancer development.
 We provide evidence that FDI has the potential to facilitate early cancer detection, tissue-of-origin prediction, cancer subtyping, and prognosis. 
This code provides modules for **FDI calculation, dispersed regions identification and evaluation**, as well as modules for **K-fold cross-validation, independent validation and tissue-of-origin inference for cancer classification**.



## Overview of FDI model
<p align="center">
  <img src="/FDI/fig/Figure1.png" width="100%"/> 
</p>

## Table of Contents
 - [Environment](#Environment)
 - [Preparation work](#Preparation)
 - [Data preprocessing](#Data_preprocessing)
 - [Identifying dispersed regions](#Identifying_dispersed_regions)
 - [Diagnostic](#Diagnostic)
 - [Citation](#citation)

## 1 Environment

First, Please install MATLAB
Then, get the code.
```
git clone --recursive https://github.com/alcindor819/FDI_code_MATLAB.git
```

## 2 Preparation

```
cd FDI_code_MATLAB/FDI/Basic_info/
```
```
wget -c https://zenodo.org/record/3928546/files/GC.zip
```
```
wget -c https://zenodo.org/record/3928546/files/mappability.zip
```
Then, unzip them.
```
unzip GC.zip
```
```
unzip mappability.zip
```
```
cd FDI_code_MATLAB/FDI/Basic_info//Dispersed_region_identify/bed_folder/
wget -c BH01.bed
```
There is a sample bed file in /FDI/Basic_info//Dispersed_region_identify/bed_folder/ .

## 3 Data_preprocessing
You need to convert all the bed files you want to use into mat files that can be used directly by the code.
First, cd to the Data Preprocessing folder.
```
cd FDI\Data preprocessing
```
Calling the user_processing function.m .
```
user_processing(bed_info_path,bed_folder,numWorkers,res_folder,name_col)
```
You need an xlsx file with all the sample information, including sample names, positive and negative labels of the samples.
You should give the position of the column with the sample name in the xlsx file.
```
bed_info_path = '\FDI\bed\Cristiano_dataset_ID.xlsx';%A sample could be this
name_col = 1;
```
You need a folder dedicated to the bed files.
```
bed_folder = '\FDI\bed\';%A sample could be this
```
You will need a folder dedicated to the generated mat files
```
res_folder = '\FDI\mat\';%A sample could be this
```
You can choose numWorkers according to the performance of your computer, the higher its value the more parallel resources, the faster it runs.
```
numWorkers = 2;%A sample could be this
```
After the code finishes, you'll get the mat format data for all the samples in the res_folder.

## 4 Identifying_dispersed_regions

You Need to call the function user_identify_dispersedregions.m.
```
cd /FDI/Dispersed_region_identify/
```
Calling this function in MATLAB.
```
user_identify_dispersedregions(bed_folder,mat_path,distribution,P,FDR,res_path,X,Y,chr_n,numWorkers)
```
The **parameters** to be passed to the function user_identify_dispersedregions are as follows:
The contains the name of the **parameter: bed_folder,mat_path,distribution,P,FDR,res_path,X,Y,chr_n,numWorkers**, an example of the parameter and the comment on the parameter.
```
bed_folder = '/FDI/Dispersed_region_identify/bed_folder/';%This folder contains the .bed files.

bed_name = 'SRR16574631';%This is the .bed filename prefix for the samples needed to identify the dispersed region.

bed_name = 'SRR16574631';%This is the .bed filename prefix for the samples needed to identify the dispersed region.

distribution = 'Beta';%This is the type of distribution, and you can choose a distribution such as normal depending on the actual distribution of the data.

P =0.01;%Global and local P-value.

FDR = 0.05;%FDR threshold.

res_path = '\FDI\Dispersed_region_identify\call_new\';%Folders for dispersed regions and TSS,CTCF results.

X = 0.5;%A parameter of FDI, default is 0.5, range is (0-1].

Y = 10;%A parameter of FDI, default is 10, range is [10-100].

chr_n = 1:22;%Identify chromosomes from 1 to 22.

numWorkers = 2;%The number of cores used in parallel is recommended to be no higher than 2 for an average local computer.
```
Each completed chromosome returns the prompt: current chromosome completed.
When it finishes running, it generates a folder named bed_name in res_path, where the results are stored.

**Results for identifying dispersed regions.**
<p align="center">
  <img src="/FDI/fig/Example_tssctcf.png" width="100%"/> 
</p>

## 5 diagnostic
### 5.1 
