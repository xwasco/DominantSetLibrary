#Dominant Set Library
===========

The Dominant Set Library (DSLib) is an open source Matlab library implementing the Dominant Set (DS) clustering method. The DS is a graph-based clustering technique rooted in the evolutionary game theory that starts gaining lots of interests in the computer science community. It has been originally introduced in (Dominant Sets and Pairwise Clustering , M.Pavan & M.Pelillo, PAMI 2007) and, thanks to its duality with game theory, has been explored in many directions not only related to clustering problems. For example applications in matching, segmentation, classification, biomedical imaging and network analysis are common in literature using the original approach. In this package we propose not only the original implementation but a still growing collections, as more comprehensive as possible, of methods and hacks from different researchers based on the original core.

This package implements the following papers, the plan is to include more in the near future:

[1] M. Pavan and M. Pelillo. Dominant sets and pairwise clustering. PAMI 2007

If you use the Infection Immunization Dynamics (dynType = 1) please cite this work:

[2] S. Rota Bulo, and I. M. Bomze.  Infection and immunization:  a new class of evolutionarygame dynamics. Games and Economic Behaviour, vol.  71, pp.  193–211, 2011.Special issue in honor of John F. Nash, jr.

If you use the prototype generation method please cite:

[3] S. Vascon et al. Using Dominant Sets for k-NN Prototype Selection. ICIAP 2013

If you use this implementation of the Dominant Set please cite this:
@misc{vascon2017dominantsetlib,
  title = {{A Matlab library for the Dominant Set clustering }},
  author = {Sebastiano Vascon}, 
  url = {https://github.com/xwasco/DominantSetLibrary},
}


## Beginning:
Execute the 'setup.m' script in order to set the path and compile the mex files. In the folder demo there are three examples.

The dominantset.m file perform the clustering. The helper associated to the dominantset.m file (type: help dominantset in your Matlab environment) will show you all the parameters and some hints.

## Requisite:
No particular requisite are needed, the library run on any environment (Windows, Linux, Mac) and it has no dependencies. It has been tested on Matlab 2014b,2016b,2018b in Windows/Linux and MacOSX. Some components are written in C/mex for optimization purpose but equivalent Matlab code are available without the needing of compiling it.

### What’s inside the package ?
In this version of the library we have implemented the following papers: Pavan and Pelillo (2003, 2007); Rota Bulo and Bomze (2011); Vascon et al. (2013) and a still growing set of papers will be added in the future.

### Dynamical systems
The method uses a dynamical system to extract the dominant sets. The library comes with three different type of dynamics, RepDyn, ExpRep-Dyn, InfImm, and it can be easily extended to upcoming new optimization method.

## Pratical usage
In this section we show how to setup a basic running example to cluster a generated synthetic dataset.

Step 1: Generating data and build the affinity matrix A:
```
rng('default'); % For reproducibility
cx = [1 1;10 10 ;20 20]; %center of the clouds of points
npts=100; pts= repmat(cx,npts,1) + randn(npts*size(cx,1),2);
A=pdist(pts); %pairwise Euclidean distances
s=3*var(A); %an euristic to find sigma
A=exp(−A./s); A=squareform(A); %from distance to similarity
A=A.*not(eye(size(A))); %the graph should not have self−loops
```

Step 2: Choosing the evolutionary dynamics and the DS parameters:
```
dynType=1; %0=Replicator Dynamics, 1=InfectionImmunization 2=Exponential replicator dynamics
precision=1e−6; %the precision required from the dynamical system
maxIters=1000; %number of maximum iteration of the dynamical system
x=ones(size(A,1),1)./size(A,1); %starting point of the dynamical system
theta=1e−5; %threshold used to extract the support from x.
```

Step 3: Call the clustering method
```
[C]=dominantset(A,x,theta,precision,maxIters,dynType);
```

alternatively the clustering method can be called providing only the similarity matrix and using the default parameters.
```
[C]=dominantset(A);
```

Step 4: Show the cluster results
```
gscatter(pts(:,1),pts(:,2),C); %show the points colored by cluster
```
The package comes with more complex examples and a deep inline documentation.
