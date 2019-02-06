function [ C,stat,S,P] = dominantset( A,x0,supportThreshold,precision,maxIters,dynType,kappa,pts)
%DOMINANTSET Extract all the dominant sets in a given weighted graph represented in terms of its weighted adjacency matrix.
%
% This package implements the following papers, the plan is to include more in the near future:
%
% [1] M. Pavan and M. Pelillo. Dominant sets and pairwise clustering. PAMI 2007
%
% If you use the Infection Immunization Dynamics (dynType = 1) please cite this work:
% [2] S. Rota Bulo, and I. M. Bomze.  Infection and immunization:  a new class of evolutionarygame dynamics.
% Games and Economic Behaviour, vol.  71, pp.  193–211, 2011.Special issue in honor of John F. Nash, jr.
%
% If you use the prototype generation method please cite:
% [3] S. Vascon et al. Using Dominant Sets for k-NN Prototype Selection. ICIAP 2013
%
%
% If you use this library please cite the following paper:
%
%
%
% Input:
%   A                   a matrix of pairwise affinities (nxn) where n is the number of elements.
%                       The matrix must NOT be of type sparse if dynType=1;
%   x0                  The initial population state. If omitted or empty then x0 is located in the
%                       baricenter of the n-dimensional symplex.
%   supportThreshold    The thershold used to extract the support from the population vector x (default 1e-4);
%   precision           The maximum population distance (Euclidean) between two
%                       successive steps to consider the dynamics in equilibrium (default 1e-8).
%   maxIters            Maximum number of iterations of the dynamical systems (both Infection Immunization of
%                       Replicator). Default = 1000.
%   dynType             Set to 0 for the Replicator Dynamics, 1 for the InfectionImmunization dynamics and 2 for the Exponential Replicator Dynamics (default 1).
%                       For large dataset (n>1000) the use of dynType=1 is highly reccomended.
%
%   kappa               In case that the dynType=2 then kappa is used. Kappa represents the "accelletarion" of the exponential dynamics (default 1).
%
%   pts                 The nxd matrix of points. When set the algorithm generates a set of prototypes that can be used as substitutes to the
%                       one of [1]. The prototypes has been used in [3].
%
% Output:
%
%   C                   An nx1 vector with the assigments of points to cluster. C(i) returns the cluster assigned to the i-th node.
%                       The value C(i)=0 means that the i-th node share no similarites with any node and thus can be considered as an outlier.
%
%   stat                Is an mx4 matrix (where m is the number of dominant sets) in which each row correspond to a dominant set.
%                       The content of each row is the following [Cohesiveness,Iterations,Precision,CentroidID]:
%                       - Cohesiveness is a value of compactness/goodness of a cluster (higher the value more similar are the elements within it.
%                       - Iterations are the number of iterations needed by the dynamical system to converge.
%                       - Precision correspond the to the difference of two successive steps in the dynamical system.
%                       - CentroidID is the id of the cluster centroid (the element with maximum value in the characteristic vector)
%
%   S                   A vector of structs sthat contains the variables C and Stat organized as follow:
%                       S{i}.index contains the index of the elements in the i-th dominant set.
%                       S{i}.cohesiveness is the cohesiveness of the i-th cluster
%                       S{i}.iters the number of iteration needed to extract the i-th cluster
%                       S{i}.precision the precision reached by the dynamical system (difference of two successive steps).
%                       S{i}.ClusterID the id of the centroid of the i-th cluster
%
%   P                   A vector of structs containing the generated prototype per cluster as in [3]. The set of points (pts) must be also provided.
%                       P{i}.avg is the centroid generated as the mean of points of cluster i.
%                       P{i}.wavg is the centroid computed as the weighted mean of points of cluster i.
%
% MIT License
% 2017, Sebastiano Vascon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<8
    pts=[];
end

if nargin<7 || isempty(kappa)
    kappa=1;
end

if nargin<6 || isempty(dynType)
    dynType=1;
end

if nargin<5 || isempty(maxIters)
    maxIters=max(1000,10^(ceil(log10(size(A,1)))+1));
end

if nargin<4 || isempty(precision)
    precision=1e-8;
end

if nargin<3  || isempty(supportThreshold)
    supportThreshold=max(1e-4,10^-(ceil(log10(size(A,1)))+1));
end

if nargin<2 || isempty(x0)
    x0=ones(size(A,1),1); x0=x0./sum(x0);
end

if dynType==1
    if issparse(A)
        warning('DSLib: Matrix A cannot be of type SPARSE when dynType=1. Matrix A is going to be converted to full.');
        A=full(A);
    end
end

if nargout>=4 && isempty(pts)
    error('DSLib: if the set of prototypes are requested as output it is necessary to provide the points descriptors.');
end

C=zeros(size(A,1),1);
stat=[];

if numel(nonzeros(A))==0 || isempty(A)
    warning('DSLib: The similarity matrix is empty!');
    stat=[0,0,0,0];
end

d=1:size(A,1);
cid=1;
S=[];
P=[];

while ~isempty(A) && numel(nonzeros(A))>0
    x=x0(d); x=x./sum(x);
    
    if dynType==1
        if exist('inImDynC','file')==3
            %if the mex exist use it
            [x,iters,nasherror]=inImDynC(A,x,precision,int32(maxIters));
        else
            [x,iters,nasherror]=inImDynM(A,x,precision,int32(maxIters));
        end
    elseif dynType==2
        [x,iters,nasherror]=ExpRepDyn(A,x,precision,maxIters,kappa);
    else
        [x,iters,nasherror]=RepDyn(A,x,precision,maxIters);
    end
    
    cohes=(x'*A)*x;
    xid=x>supportThreshold;
    
    if sum(xid)==0
        warning('DSLib: The support of vector x is empty, consider lower the support threshold\n');
        break;
    else
        [~,p]=max(x); p=d(p); %get the centroid id
        idx=d(xid);
        
        C(idx)=cid;
        
        if nargout>=2
            stat(cid,:)=[cohes,iters,nasherror,p];
        end
        
        if nargout>=3
            S{cid,1}.index=idx;
            S{cid,1}.cohesiveness=cohes;
            S{cid,1}.iters=iters;
            S{cid,1}.precision=nasherror;
            S{cid,1}.centroid=p;
        end
        
        if nargout>=4
            P.avg(cid,:)=mean(pts(idx,:));
            P.wavg(cid,:)=sum(repmat(x,1,size(pts,2)).*pts(d,:));
        end
        
        %remove the unused edges and nodes, this improve the efficiency
        d(xid)=[];
        A(xid,:)=[];
        A(:,xid)=[];
        cid=cid+1;
    end
end

end
