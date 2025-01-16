## Run basic representational similarity analyses

This example provides instructions to run RSA for odor valence and salience
1. Open ARC_RSA_analyses.m and change variable <mainroot> to the directory of the repository on the local computer
2. Make sure the variable demomode is set to "true"
3. Set the name <modelname> of the results folder 
4. Run ARC_RSA_analyses.m 

It should produce the following RSA results in repository>results>modelname
Results may very slightly differ because of stochasticity and noise.
1. ARC_RSAwt.png: Beta weights of RSA coefficients of valence and salience
2. imagescr.png: Representational similarity matrices for all ROIs and subjects
3. voxprop.png: Fraction of voxels tuned to appetitive or aversive only pleasantness
4. voxwise_similarity.png: Negative or null correlation between RSA coefficients for the pleasantness of appetitive and aversive odors, across all voxels in different regions 
5. ARC_dens.png: . Distribution of voxel-wise RSA coefficients for the pleasantness of appetitive and aversive odors. 