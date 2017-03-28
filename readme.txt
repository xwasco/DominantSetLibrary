Beginning:
Execute the 'setup.m' script in order to set the path and compile the mex files. In folder demo there are three example for the usage of the package.

The dominantset.m file perform the clustering. The helper associated to the dominantset.m file (type: help dominantset in you Matlab environment) will show you all the paramters and some hints.


Requisite
No particular requisite are needed, the library run on any environment (Windows, Linux, Mac) and it has no dependencies. It has been tested on Matlab 2014b and 2016b in Windows/Linux and MacOSX. Some components are written in C for optimization purpose but equivalent Matlab code are available without the needing of compiling them.

What’s inside the package ?
In this version of the library we have implemented the following papers: Pavan and Pelillo (2003, 2007); Rota Bulo and Bomze (2011); Vascon et al. (2013) and a still growing set of papers will be added in the future.

Dynamical systems
The method uses a dynamical system to extract the dominant sets. The library comes with three different type of dynamics, RepDyn, ExpRep-Dyn, InfImm, and it can be easily extended to upcoming new optimization method for the program (1).

Pratical usage
In this section we show how to setup a basic running example to cluster a generated synthetic dataset.

Step 1: Generating data and build the affinity matrix A:

rng('default'); % For reproducibility
cx = [1 1;5 5 ;8 8]; %center of the clouds of points
npts=100; pts= repmat(cx,npts,1) + randn(npts*size(cx,1),2);
A=pdist(pts); %pairwise Euclidean distances
s=3*var(A); %a heuristic to find sigma
A=exp(-A./s); A=squareform(A); %from distance to similarity
A=A.*not(eye(size(A))); %the graph should not have self−loops

Step 2: Choosing the evolutionary dynamics and the DS parameters:

dynType=1; %0=Replicator Dynamics, 1=InfectionImmunization 2=Exponential replicator dynamics
precision=1e-6; %the precision required from the dynamical system
maxIters=1000; %number of maximum iteration of the dynamical system
x=ones(size(A,1),1)./size(A,1); %starting point of the dynamical system
theta=1e-5; %threshold used to extract the support from x.


Step 3: Call the clustering method
[C]=dominantset(A,x,theta,precision,maxIters,dynType);

alternatively the clustering method can be called providing only the similarity matrix and using the default parameters.
[C]=dominantset(A);

Step 4: Show the cluster results
scatter(pts(:,1),pts(:,2),5,C); %show the points colored by cluster

The package comes with more complex examples and a deep inline documentation.

-----

Sebastiano Vascon, 2017
