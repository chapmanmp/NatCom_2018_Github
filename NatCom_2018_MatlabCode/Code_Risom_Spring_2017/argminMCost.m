%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identifies predicted data with fixed dynamics matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mhat = argminMCost(Ahat, lambda, T, M)

[n,TW] = size(M); W = TW/T;

J_Exists = GetExistingSamples(M);

cvx_begin

    variable MVar(n,TW) nonnegative
    
    minimize( norm(MVar(:,W+1:TW) - Ahat*MVar(:,1:(T-1)*W),'fro') + lambda*norm(Ahat,'fro') + norm(MVar(:,J_Exists)-M(:,J_Exists),'fro') )
        
cvx_end

Mhat = MVar;

%X = M(:,1:(T-1)*W);

%XPlus = M(:,W+1:TW);




