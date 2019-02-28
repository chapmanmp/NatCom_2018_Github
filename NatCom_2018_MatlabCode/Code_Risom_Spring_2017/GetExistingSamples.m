%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identifies indices of samples that are completely specified with observations of all data types.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_Exists = GetExistingSamples(M)

[~,TW] = size(M);

J_Empty = GetDiscardedSamples(M); 

J_Exists = setdiff(1:TW,J_Empty);

%VERIFICATION
%This function was tested during the week of 4/12.
