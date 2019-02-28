%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Completes discarded samples with the average of existing samples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sample index = column index of data matrix M

function M_hat = GetCoarseEstimate(M,T)

[~,TW] = size(M); %TW = #Time points * #Wells = #Samples

W = TW/T;

J_Empty = GetDiscardedSamples(M);

M_hat = M;

for k = 0:T-1
    
    Jk = (1+k*W):(k+1)*W; %Sample indices at time k
    
    Jk_Empty = intersect(Jk, J_Empty); %Sample indices at time k missing measurements.
    
    Jk_Exist = setdiff(Jk, Jk_Empty); %Sample indices at time k with all measurements.
    
    for index = 1:length(Jk_Empty)
        
        j = Jk_Empty(index);
        
        M_hat(:,j) = 1/length(Jk_Exist)*sum(M(:,Jk_Exist),2); %Sum over all columns.
        
    end
    
end

%VERIFICATION
%This function was tested during the week of 4/12.

% %Test code: Do sample indices at time k make sense? Let T = 7.
% W = 15;
% if(sum(1:W == 1:15)==W), display('k = 0 is correct');
% else display('k = 0 is NOT correct');
% end
% 
% if(sum(W+1:2*W == 16:30)==W), display('k = 1 is correct');
% else display('k = 1 is NOT correct');
% end
% 
% if(sum(5*W+1:6*W == 76:90)==W), display('k = 5 is correct');
% else display('k = 5 is NOT correct');
% end
% 
% if(sum(6*W+1:7*W == 91:105)==W), display('k = 6 is correct');
% else display('k = 6 is NOT correct');
% end













