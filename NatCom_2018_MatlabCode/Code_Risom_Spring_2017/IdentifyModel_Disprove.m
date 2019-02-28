%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uses alternating minimization wrapped around l2-regularized least-squares to identify optimal dynamics in the presence of discarded data.
%2 phenotypic states: green, *not* green. 
%Switching specified on or off
%Death gain constraints, or lack thereof, specified via AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [AStar, MStar, e_lambda, indexStar] = IdentifyModel_Disprove(M, SetWExists, T, SetLambda, epsilon, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO)

[~,TWTotal] = size(M); WTotal = TWTotal/T;

%Split data set into a test set and a training set. Test set contains wells without missing data for all time points, all conditions of interest.
[MTest, MTrain] = GetTestTrainMatrices(M, SetWExists, WTotal);

%Get optimal regularizer using cross-validation.
[e_lambda, indexStar] = GetBestRegularizer_Disprove(MTest, MTrain, T, SetLambda, epsilon, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO);

%Get optimal dynamics and complete data set.
lambdaStar = SetLambda(indexStar);
[AStar, MStar] = GetDynamicsAndCompleteData_Disprove(M, lambdaStar, T, epsilon, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO);



