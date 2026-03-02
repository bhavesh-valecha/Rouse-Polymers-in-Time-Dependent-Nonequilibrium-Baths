# Rouse Polymers in Time Dependent Nonequilibrium Baths

This repository contains code and data for [Rouse Polymers in Time Dependent Nonequilibrium Baths](%add_arXiv_link_here). 

Authors: [Bhavesh Valecha](https://www.uni-augsburg.de/de/fakultaet/mntf/physik/groups/theo2/team/phds/bhavesh-valecha/), [Jens-Uwe Sommer](https://www.ipfdd.de/de/forschung/bereich-theorie-der-polymere/people/prof-dr-jens-uwe-sommer/) and [Abhinav Sharma](https://www.uni-augsburg.de/de/fakultaet/mntf/physik/groups/theo2/team/professors/prof-dr-abhinav-sharma/)

## Abstract
Directed transport is a characteristic feature of numerous biological systems in response to signals such as nutrient and chemical gradients. These signals are often time-dependent owing to the high complexity of interactions in these systems. In this study, we focus on the steady-state behavior of polymeric systems responding to such time dependent signals. We model them as ideal Rouse polymers submerged in a time-dependent nonequilibrium bath, which is described by a spatially and temporally varying self-propulsion wave field. Through a coarse-graining analysis, we show that these polymers display rich emergent response to the temporal stimuli as a function of their length and topology. In particular, short polymers drift against the self-propulsion wave, whereas, longer polymers drift in the direction of the wave signal. Furthermore, monomers connected in a star or ring topology also surf with the wave crests, as opposed to a fully connected structure, which drifts with the wave troughs. We confirm these analytical predictions with robust numerical simulations, showing that response of polymeric systems to temporal stimuli can be controlled by the topology or the length of the polymer.
        
## Description

1. Codes/ - folder containing all codes and header files for the Langevin dynamics simulations.

2.  steady_state_density_polymer_waves.cpp - code to generate the steady state density profile for an active polymer in a traveling activity wave.

3. drift_velocity.cpp - code to generate the drift of the active polymer for different wave speeds of the activity signal.

4. particle_polymer_waves.h - header file that defines the class polymer used in the codes steady_state_density_polymer_waves.cpp and drift_velocity.cpp.

5. random_polymer_waves.h - header file that defines the random number generator used in the codes steady_state_density_polymer_waves.cpp and drift_velocity.cpp.

6. Force_chain.h - header file for calculating the forces between monomers for a linear polymer.

7. Force_ring.h - header file for calculating the forces between monomers for a ring polymer.

8. Force_star.h - header file for calculating the forces between monomers for a star polymer.

9. Force_clique.h - header file for calculating the forces between monomers for a clique polymer.

10. plot_figures.ipynb - notebook used to plot Figs. 2-5 in the main text of the article.

11. Data/ - folder containing all the data files used to create Figs. 2-5 in the main text of the article.


## Citation

If you found this code or paper helpful, please cite us as:
%ADDCITATION

## Contact

For any questions or comments, please contact us at <bhavesh.valecha@uni-a.de>
