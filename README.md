# Mapping bimodal colloidal depletion systems onto the monodisperse gelation phase diagram
This is a simulation and analysis pipeline for creating a phase diagram of bimodal colloidal gels.

Results from this project are [available on arXiV](https://doi.org/10.48550/arXiv.2610.33394) and are currently under review.

Simulations can be run using the [sim-scripts](./sim-scripts) and the DPDMorse extension for HOOMD-blue v4.2.1, available as [hoomd4.2.1-mod](https://github.com/procf/hoomd4.2.1-mod).

To confirm that a system is at equilibrium before gelation, use [analysis-eq](./analysis-scripts/analysis-eq)

## What to expect
For systems of ~10,000 colloidal particles the following scripts are available:

[analysis-mono-colloids](./analysis-scripts/analysis-mono-colloids) and [analysis-bi-colloids](./analysis-scripts/analysis-bi-colloids):
- Coordination number distribution and average coordination number, average and decomposed by population
- Mean squared displacement (MSD), average and decomposed by population
- Void size distribution and average void size
- Voronoi volume distribution
- Network edgelist, the size of the largest connected component, and the physical network diameter ("span" as a proportion of total box size)
- Cluster size distribution for components not in the LCC 

[collate\*.py](./analysis-scripts):
- Collate data for comparison and plotting

[phase-classification\_and\_plot.ipynb](./analysis-scripts/phase-classification_and_plot.ipynb)
- Calculate reduced second virial coefficient
- Phase classification
- Plotting

Note: Full analyses typically take less than 20min

In bimodal depletion gels, colloid-colloid interactions scale with particle size. Therefore, the attraction strength for S-S, S-L, and L-L contacts will be different.
For size ratio 1:2 this scales roughly with the arithmetic mean, such that $D_0^{SL} = 1.5 D_0^{SS}$ and $D_0^{LL} = 2 D_0^{SS}$.

## Software/package requirements
In this project, the following packages are actively used:
1. GNU Fortran (GCC) 11.4.1
2. `python` v3.10.16 
3. `gsd` v3.2.1
4. `numpy` v1.24.40
5. `pandas` v2.2.3
6. `networkx` v3.4.2
7. `scipy` v1.15.3
8. `node2vec` v0.5.0
9. `umap-learn` v0.5.9
10. `scikit-learn` v1.7.2

## Hardware/OS tested
The program was tested on a single HPC-node running Rocky Linux 9.3 (kernel 5.14).

## Background

Colloidal gels have complex mechanical and rheological properties that emerge from their underlying network structure. In monodisperse systems 
the percolation of a gel network can be predicted from three state variables: the total particle volume fraction, the interparticle attraction 
range, and interparticle attraction strength. We use simulations to investigate how this relationship changes in gels made with two particle sizes 
(size ratio 1:2). We find that bimodal composition changes phase behavior, and that this effect can be mapped onto the monodisperse phase diagram 
using the contact-weighted reduced second virial coefficient; however, the multiscale network structure of these systems is not uniform and still 
depends on the specific details of each bimodal composition.

Here's the abstract: 

*Many industries rely on models of colloidal gelation and phase behavior to predict how processing equipment will handle new product formulations. 
In colloidal systems with nearly uniform particle size it is possible to predict the formation of a gel from only three state variables; however, 
the extension of this model to bimodal and multimodal systems is not well described. Here we present a systematic simulation study of the attractive 
bimodal depletion gel phase diagram. For these systems, we find that gelation is highly dependent on total volume fraction, attraction conditions, 
and bimodal composition. Attraction and composition can be effectively collapsed into a single parameter using the contact-weighted second virial 
coefficient. This successfully projects the bimodal phase diagram onto the established monodisperse phase diagram, including accurately describing 
most bulk assembly properties; however, the underlying network structure of these gels remains highly dependent on the details of bimodal 
composition and attraction.*


## Contributors
This project was in collaboration with the [Soft Matter Engineering Laboratory](https://smel.eng.uci.edu/) at the University of California, Irvine. \
This work was done by [Rob Campbell](https://scholar.google.com/citations?user=i8S54zYAAAAJ&hl=en), Calvin (Ziye) Zhuang, 
[Ali Mohraz](https://scholar.google.com/citations?user=pW80NaAAAAAJ&hl=en), and [Safa Jamali](https://scholar.google.com/citations?user=D1asaYIAAAAJ&hl=en).\
Authors acknowledge support from the National Science Foundation (PMP-2025613) and NASA ROSES FINESST (80NSSC23K0015).
