particle-tracking
=================
#### A linking algorithm for particle tracking in n-dimensions, implementing a kinematic model and a memory feature to account for occasional misses.
 
Multi particle tracking is a common problem, often separated in two parts:

1. Particle detection and localisation
2. Linking these sets of positions in time-lapse sequences into trajectories.

The python code provided here implements an algorithm to address the second one; the matlab script accompanying Kelley and Ouellette (Am. J. Phys., 79(3):267, 2011) served as a starting point. 

Among the modifications and features which have been introduced:
+ implemented in Python and supports pandas DataFrame input and output
+ generalised to n-dimensions
+ the kinematic model, in which future positions are inferred from the already linked past positions, accounts for velocities and accelerations based on finite-differences
+ a memory feature was added to account for the occasional loss of tracers, and it was optimised for better performance. 

The algorithm implemented here accounts for the physical process of particles advected by a smooth chaotic flow and for the uncertainties. These arise from the chaotic in time nature of the flow – “physical noise” – as well as from localisation and past linking errors – “experimental noise”.

##### Related projects and publications:
+ E. Afik and V. Steinberg. [Pair dispersion in a chaotic flow reveals the role of the memory of
initial velocity](http://arxiv.org/abs/1502.02818). _ArXiv e-prints arXiv:1502.02818_. submitted.
+ E. Afik. [Robust and highly performant ring detection algorithm for 3d particle tracking using 2d microscope imaging](http://www.nature.com/articles/srep13584). Sci. Rep. 5, 13584; doi: 10.1038/srep13584 (2015)
