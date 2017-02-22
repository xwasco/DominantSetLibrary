function [x,iters,NashError]=ExpRepDyn(A,x,toll,maxiter,kappa)
    NashError=2*toll;
    iters=0;
    while (toll<=NashError && iters<maxiter)
        old_x=x;
        x=x.*exp(kappa*(A*x));
        x=x./sum(x);

        NashError= norm(x-old_x);

        iters=iters+1;
    end

end

