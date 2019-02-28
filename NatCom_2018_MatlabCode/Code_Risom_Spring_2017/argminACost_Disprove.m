%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identifies a dynamics matrix with fixed predicted data via l2-regularized least-squares.
%2 phenotypic states: green, *not* green. 
%Switching specified on or off
%Death gain constraints, or lack thereof, specified via AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AVar - dynamics matrix, n - number of states, 
%noDeath - 'noDeath' phenotype, switching - on or off, AllDeathNoDeathGains_GSKBEZ - 1 or 0, EqualDeathGains_DMSO - 1 or 0 
%variable names reflect those in ConstraintsOnA_Disprove
%Mhat - predicted data, M - original data.

function Ahat = argminACost_Disprove(Mhat, lambda, T, M, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO)

[n, TW] = size(M); W = TW/T;

J_Exists = GetExistingSamples(M);

cvx_begin

    variable AVar(n,n)
    
    minimize( norm(Mhat(:,W+1:TW) - AVar*Mhat(:,1:(T-1)*W),'fro') + lambda*norm(AVar,'fro') + norm(Mhat(:,J_Exists)-M(:,J_Exists),'fro') )
    
    subject to
    
        ConstraintsOnA_Disprove
        
cvx_end

Ahat = AVar;

%X = Mhat(:,1:(T-1)*W);

%XPlus = Mhat(:,W+1:TW);




