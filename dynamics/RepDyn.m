function [x,iters,NashError]=RepDyn(A,x,toll,maxiter)
%REPDYN Replicator Dynamcs. 
%
% Input:
%   A           A pairwise nxn similarity matrix (with zero diagonal)
%
%   x           An nx1 vector in the n-dimensional simplex (it should add up to 1).
%
%   toll        The precision required from the dynamical system
%
%   maxIters    The maximum numer of iterations
%
% Output:
%
%   x           The population vector at convergence
%
%   iters       The number of iterations needed to converge
%
%   NashError   The precision reached by the dynamical system
%
% Released under: MIT License
% 2019, Sebastiano Vascon


    NashError=2*toll;
    iters=0;
    while (toll<=NashError && iters<maxiter)
        old_x=x;
        x=x.*(A*x);
        x=x./sum(x);

        NashError= norm(x-old_x);

        iters=iters+1;
    end

end
