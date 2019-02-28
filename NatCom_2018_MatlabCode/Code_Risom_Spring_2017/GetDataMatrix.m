%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Combines phenotypic state counts (P), percent death (d) into a single matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.

function M = GetDataMatrix(P,d)

[NumPheno,TWTotal] = size(P); %NumPheno = number of phenotypes, TWTotal = #Time points * #Wells = Number of samples.

n = NumPheno + 1; %Number of states.

M = zeros(n,TWTotal);

for i = 1:n
    for j = 1:TWTotal
        if i==n
            M(i,j) = d(j)*sum(P(:,j)); %Sum over column j.
        else
            M(i,j) = (1-d(j))*P(i,j);
        end
    end
end

%VERIFICATION
%This function was tested during the week of 4/12.