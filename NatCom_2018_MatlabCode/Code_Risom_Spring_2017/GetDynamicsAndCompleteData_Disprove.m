%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uses alternating minimization to identify dynamics in the presence of discarded data.
%2 phenotypic states: green, *not* green.
%Switching specified on or off
%Death gain constraints, or lack thereof, specified via AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ahat, Mhat] = GetDynamicsAndCompleteData_Disprove(M, lambda, T, epsilon, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO)

[n,~] = size(M);

Abar = ones(n,n);

Mbar = GetCoarseEstimate(M,T);

while(true)
   
    AbarPLUS = argminACost_Disprove(Mbar, lambda, T, M, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO); %Identify dynamics matrix with fixed predicted data.
        
    MbarPLUS = argminMCost(AbarPLUS, lambda, T, M); %Identify predicted data with fixed dynamics matrix.
    
    if (norm(AbarPLUS - Abar,'fro')/norm(Abar,'fro') < epsilon)
        break;
    else
        Abar = AbarPLUS;
        Mbar = MbarPLUS; 
    end
     
end

Ahat = AbarPLUS;

Mhat = MbarPLUS;

%Looks good 4/21, 6/8
