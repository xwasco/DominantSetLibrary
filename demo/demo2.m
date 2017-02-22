close all
clear all

rng('default'); % For reproducibility
cx = [1 1;10 10 ;20 20]; %center of the clouds of points
npts=100; pts= repmat(cx,npts,1) + randn(npts*size(cx,1),2);
A=pdist(pts); %pairwise Euclidean distances
s=3*var(A); %max(max(A))-min(min(A)); %an euristic to find sigma
A=exp(-A./s); A=squareform(A);%from distance to similarity
A=A.*not(eye(size(A))); %the graph should not have the self-loop

dynType=1; %0=Replicator Dynamics, 1=InfectionImmunization 2=Exponential replicator dynamics
precision=1e-6; %the precision required from the dynamical system
maxIters=1000; %number of maximum iteration of the dynamical system
x=ones(size(A,1))./size(A,1); %starting point of the dynamical system
theta=1e-5; %threshold used to extract the support from x.

[C]=dominantset(A,x,theta,precision,maxIters,dynType);
gscatter(pts(:,1),pts(:,2),C);