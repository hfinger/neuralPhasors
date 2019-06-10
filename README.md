# Neural Phasors: A Neural Network And Image Segmentation Toolbox

This repository contains code that was written as part of my PhD thesis at the Institute of Cognitive Science. 
The following is a short list of some of the functionality that is implemented in this toolbox:

* A Sparse Autoencoder for unsupervised learning of features from natural images.
* Stacking several autoencoder layers to create a deep network.
* Code for ZCA-whitening of color images.
* Neural network of Kuramoto oscillators that are used for image segmentation.
* A complex-valued autoencoder implemented in the machine-learning framework Torch.
* Code to extract functional brain connectivity from EEG data using different types of coherence and causality measures.
* Probabilistic fiber tracking using diffusion tensor imaging (DTI) data.
* A GUI for paramter exploration on the Sun Grid Engine.
* Dynamic neural network simulations of Jansen-Rit neural mass models in the cortical connectome.


## Instructions to run the code:

Edit the paths in "matlab/include/dataPaths.m" to match your system setup.

In Matlab run the following script to add all relevant paths to the environment:
'''
addScriptPaths();
'''

There are several scripts in the parameters subfolder that can be used to start computations.


## List of publications that were produced using code in this repository:

* Finger, H., Gast, R., Gerloff, C., Engel,
A. K., & König, P. (in review). Probing
Neural Networks for Dynamic Switches
of Communication Pathways. PLoS Cb.
[Link](http://www.holgerfinger.com/cv/publications/2019_probing.pdf)

* Finger, H. (2017). Information Process-
ing in Neural Networks: Learning of
Structural Connectivity and Dynamics
of Functional Activation. Dissertation.
[Link](http://www.holgerfinger.com/cv/publications/2017_information.pdf)

* Finger H, *Bönstrup M, Cheng B,
Messé A, Hilgetag C, Thomalla G, et al.
(2016). Modeling of Large-Scale Func-
tional Brain Networks Based on Struc-
tural Connectivity from DTI: Compari-
son with EEG Derived Phase Coupling
Networks and Evaluation of Alternative
Methods along the Modeling Path.
PLoS Comput Biol 12(8): e1005025.
doi:10.1371/journal.pcbi.1005025.
[Link](http://www.holgerfinger.com/cv/publications/2016_structurally.pdf)

* Finger, H., & König, P. (2014). Phase
synchrony facilitates binding and
segmentation of natural images in
a coupled neural oscillator network.
Frontiers in computational neurosci-
ence, 7, 195.
[Link](http://www.holgerfinger.com/cv/publications/2014_phase.pdf)

