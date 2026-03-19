## Analyses approach

### Run representational similarity analyses
This example provides instructions to run RSA for odor valence and salience
1. Open ARC_RSA_analyses.m and change variable <mainroot> to the directory of the repository on the local computer
2. Make sure the variable demomode is set to "true"
3. Set the name <modelname> of the results folder 
4. Run ARC_RSA_analyses.m 

It should produce the following RSA results in repository>results>modelname
Results may very slightly differ because of stochasticity and noise.
1. ARC_RSAwt.png: Beta weights of RSA coefficients of valence and salience
2. imagescr.png: Representational similarity matrices for all ROIs and subjects

NOTE: p values may be imprecise because of limited number of iterations in the permutation test.
For actual results, please run the analysis with nshuff = 1000 and demomode set to false. 

### Run Decoding analyses

This example provides instructions to run decoding analyses in Figures 3-5 of the associated manuscript for valence, salience, appetitive-pleasantness, aversive-pleasantness and cross-decoding analyses.
1. Open ARC_decoding_analyses.m and change variable <mainroot> to the directory of the repository on the local computer
2. Adjust <switcher> to change decoding analysis between modes: basic = valence/salience; domainpart = appetitive and aversive coding; crossdec = crossdecoding analyses
3. Set the name <savename> of the results folder 
4. Run ARC_decoding_analyses.m

It should produce the decoding results in repository>examples
Results may very slightly differ because of stochasticity and noise.
1. ARC_decoding.png: Prediction accuracies of 
2. imagescr.png: Representational similarity matrices for all ROIs and subjects

NOTE: p values may be imprecise because of limited number of iterations in the permutation test.
