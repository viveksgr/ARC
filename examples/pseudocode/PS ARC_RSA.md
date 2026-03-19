Pseudocode for RSA Analysis

Initialize Settings
  - Set paths and parameters for data loading and analysis
  - Define behavioral and neural data paths
  - Configure analysis settings (binning, controls, etc.)

Load Data
  - Load behavioral ratings and neural data for each subject
  - Load anatomical and functional masks

Prepare Data
  - For each subject:
    - Normalize behavioral data
    - Apply masks to neural data to select relevant voxels

Construct Representational Similarity Matrices (RSMs)
  - For valence and salience:
    - Bin behavioral data based on valence and salience scores
    - Create RSMs based on binned data

Check Analysis Mode
  - If valsep is FALSE:
    - Perform RSA for valence and salience
  - If valsep is TRUE:
    - Perform RSA for appetitive and aversive domains

Compute Neural Response RSMs
  - For each anatomical region:
    - Extract neural responses
    - Bin and normalize neural data
    - Compute correlation matrix to form neural RSM

Statistical Analysis
  - Implement permutation testing to assess the significance of RSA weights
  - Apply FDR correction for multiple comparisons
  - Calculate voxel-wise and region-level statistics

Visualize Results
  - Generate bar plots, similarity matrices, and scatter density plots for visual analysis
  - Save figures and RSA maps as NIfTI files for further use

Save and Export Results
  - Store RSA statistics, weights, and permutation results
  - Export processed data and figures to specified paths


