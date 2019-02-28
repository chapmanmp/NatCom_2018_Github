%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identifies indices of samples that lack observations of any data type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J_Empty = GetDiscardedSamples(M)

[n,TW] = size(M); %TW = #Time points * #Wells = #Samples

J_Empty = [];

for j = 1:TW
    if ( M(n,j) == Inf || M(n,j) == 0 || isnan(M(n,j)) )
       J_Empty = [J_Empty j]; 
    end
end

%VERIFICATION
%This function was tested during the week of 4/12, 6/2.





