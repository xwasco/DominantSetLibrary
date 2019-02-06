function [x,iters,NashError] = inImDynM(A,x,toll,maxIters)
%INIMDYNM Infection-Immunization dynamcs. 
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
%
% If you use the Infection Immunization Dynamics please cite this work:
% [2] S. Rota Bulo, and I. M. Bomze.  Infection and immunization:  a new class of evolutionarygame dynamics.
% Games and Economic Behaviour, vol.  71, pp.  193–211, 2011.Special issue in honor of John F. Nash, jr.
%
% Released under: MIT License
% 2019, Samuel Rota Bulo

    iters=1;
    NashError=2*toll;

    g = A*x;
    while iters<maxIters && NashError>=toll
        r = g - (x'*g);
        NashError=norm(min(x,-r));
        
        i = selectPureStrategy(x,r);
        den = A(i,i) - g(i) - r(i);
        do_remove=0;
        if r(i)>=0
            mu = 1;
            if den<0
                mu = min(mu, -r(i)/den);
                if mu<0
                    mu=0;
                end
            end
        else
            do_remove=1;
            mu = x(i)/(x(i)-1);
            if den<0
                [mu , do_remove] = max([mu -r(i)/den]);
                do_remove=do_remove==1;
            end
        end
        tmp = -x;
        tmp(i) = tmp(i)+1;
        x = mu*tmp + x;
        if(do_remove)
            x(i)=0;
        end;
        x=abs(x)/sum(abs(x));
        g = mu*(A(:,i)-g) + g;
        
        iters=iters+1;
    end
end

function [i] = selectPureStrategy(x,r)
    index=1:length(x);
    mask = x>0;
    masked_index = index(mask);
    [~, i] = max(r);
    [~, j] = min(r(x>0));
    j = masked_index(j);
    if r(i)<-r(j)
        i = j;
    end
end
