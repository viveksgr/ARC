Pseudocode for Decoding Analysis

Initialize Settings
  - Set paths for data and masks
  - Define configurations like bin sizes and mask applicability
  - Determine the decoding strategy (e.g., Cross-decoding)

Prepare Environment
  - Load necessary libraries and functions
  - Configure output directories and logging

Load Data for each subject
  - Load behavioral and neural data
  - Normalize and median split behavioral ratings if required
  - Apply anatomical and functional masks to neural data

Construct Neural and Behavioral Vectors for each ROI
  - Extract neural activity data
  - Adjust for noise and signal cutoffs as configured
  - Optionally apply PCA or similar dimensionality reduction

Prepare Data for Decoding
  - Discretize behavioral data if necessary
  - Segment data based on conditions (e.g., positive vs negative valence)

Set Up Decoding
  - Prepare neural data matrices for various decoding schemes
  - Configure labels for decoding based on behavioral outcomes
  - Handle cross-decoding setups where training on one condition and testing on another

Perform Decoding
  - Execute decoding using configured models (e.g., SVM, regression models)
  - Optionally perform permutation tests for statistical validation

Calculate and Store Results
  - Compute decoding accuracy, statistical tests, and other relevant metrics
  - Log results for each subject and ROI

Visualize and Save Outputs
  - Generate figures such as bar plots to visualize decoding accuracies
  - Save figures and statistical results to specified paths

Clean Up
  - Optionally clear variables to free up memory
  - Save final configurations and results for documentation

