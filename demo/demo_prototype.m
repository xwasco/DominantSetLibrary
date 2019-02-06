%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released under: MIT License
% 2019, Sebastiano Vascon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all

clc

addpath(genpath('../'));

rng('default');

%% Create a set of vertices (4 clusters)
numPts=100;
cx = [0 0;20 20 ;0 20; 20 0]; %center of the clouds of points
pts= repmat(cx,numPts,1) + randn(numPts*size(cx,1),2);

%% Create the edge-weighted graph (affinity matrix)
D=pdist(pts); %all the pairwise distances
sigma=3*var(D); %an heuristic to tune sigma
A=squareform(exp(-D./sigma)); %compute the affinity matrix of the graph as a Gaussian Kernel

A=A.*not(eye(size(A))); %no self-loop, set the diagonal to zero IMPORTANT !

%% Plot the clusters obtained by the different dynamical systems.
figure;
subplot(2,2,1);
scatter(pts(:,1),pts(:,2),5);
title('Data Points');

subplot(2,2,2);
[C,stat,S,P]=dominantset(A,[],[],[],[],0,[],pts);
gscatter(pts(:,1),pts(:,2),C);
hold on 
scatter(pts(stat(:,4),1),pts(stat(:,4),2),50,'d');
scatter(pts(stat(:,4),1),pts(stat(:,4),2),50,'d');
scatter(P.avg(:,1),P.avg(:,2),50,'^');
scatter(P.wavg(:,1),P.wavg(:,2),50,'*');
title(['Replicator Dynamics (#' num2str(numel(S)) ' DS)']);

subplot(2,2,4);
[C,stat,S,P]=dominantset(A,[],[],[],[],2,10,pts);
gscatter(pts(:,1),pts(:,2),C);
hold on 
scatter(pts(stat(:,4),1),pts(stat(:,4),2),50,'d');
scatter(P.avg(:,1),P.avg(:,2),50,'^');
scatter(P.wavg(:,1),P.wavg(:,2),50,'*');
title(['Exponential Replicator Dynamics (#' num2str(numel(S)) ' DS)']);

subplot(2,2,3);
[C,stat,S,P]=dominantset(A,[],[],[],[],1,[],pts);
gscatter(pts(:,1),pts(:,2),C);
hold on 
scatter(pts(stat(:,4),1),pts(stat(:,4),2),50,'d');
scatter(P.avg(:,1),P.avg(:,2),50,'^');
scatter(P.wavg(:,1),P.wavg(:,2),50,'*');
title(['Infection Immunization (#' num2str(numel(S)) ' DS)']);

