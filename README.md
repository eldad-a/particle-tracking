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
+ E. Afik and V. Steinberg. [On the role of initial velocities in pair dispersion in a microfluidic chaotic flow](https://www.nature.com/articles/s41467-017-00389-8). _Nature Communications_ __8__, Article number: 468 (2017) [doi: 10.1038/s41467-017-00389-8](http://dx.doi.org/10.1038/s41467-017-00389-8).
+ E. Afik and V. Steinberg. A Lagrangian approach to elastic turbulence in a curvilinear microfluidic channel. figshare [doi: 10.6084/m9.figshare.5112991](http://dx.doi.org/10.6084/m9.figshare.5112991) (2017).
+ E. Afik. [Robust and highly performant ring detection algorithm for 3d particle tracking using 2d microscope imaging](http://www.nature.com/articles/srep13584). Sci. Rep. 5, 13584; doi: 10.1038/srep13584 (2015)
+ [ridge-directed-ring-detector](https://github.com/eldad-a/ridge-directed-ring-detector) -- a Cython implementation of the ring detection algorithm for 3d particle tracking
+ [natural-cubic-smoothing-splines](https://github.com/eldad-a/natural-cubic-smoothing-splines)
