# SFP project scripts

Exhaustive codebase to run the pleasantness related analyses on the NEMO dataset

## 1. System Requirements

Operating System: Windows

MATLAB Version: R2023B

### Dependencies:
1. [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
2. [GLMSingle](https://github.com/cvnlab/GLMsingle)

The code was tested on a machine with the following recommended specifications:
RAM: 64GB; processor: AMD Ryzen Threadripper PRO 4.00 GHz or higher;
No non-standard hardware required.

## 2. Installation Guide
### Instructions
To install: 
```bash
git clone https://github.com/viveksgr/ARC.git
```
Install MATLAB Toolboxes directly from the source links provided. 

### Typical Install Time
The installation should take approximately a few minutes on a standard desktop computer.

## 3. Demo
Demonstration to run representational similarity analyses are provided in the examples folder. The output should be similar to corresponding analyses in Sagar et. al., 2024 "Separable and integrated pleasantness coding for appetitive and aversive odors across olfactory and ventral prefrontal cortices" (In preparation)
Runtime for the demo scripts varies by the analyses but most RSA analyses should be executable within an hour. 

### Instructions to Run on Demo Data
1. Clone the repository and add all dependencies to filepath
2. Launch MATLAB and add <common_functions> to path
3. For each demo example, make sure file path for variable <mainroot> is set as that of cloned repository on local machine.
4. Currently, demo files are provided to recreate crucial behavioral results from representational similarity analyses.

## 4. Instructions for Use and Reproduction:
To repeat the complete analyses using these scripts, acquire the complete raw [dataset](https://www.nature.com/articles/s41593-023-01414-4#data-availability) and perform suggested preprocessing.
The following scripts should be executed in order after preprocessing:
1. ARC_createsingletrials: PExtract voxelwise single trials responses using the GLM single package
2. ARC_RSA_analyses: Perform representational similarity analysis for valence, salience, appetitive and aversive pleasantness. 
3. ARC_decoding_analyses: Decoding analyses based on support vector machine to predict odor valence, salience, appetitive or aversive pleasantness.

*Detailed documentation for all scripts are in progress*

The full analysis may take approximately over 24-48 hours depending on your system specifications and GLM_single settings.

If you encounter issues, contact sgr.vivek@gmail.com.

## 5. Additional Information:
Terms of use: This content is licensed under MIT License.