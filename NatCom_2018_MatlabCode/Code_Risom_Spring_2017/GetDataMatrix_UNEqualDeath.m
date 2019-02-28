%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Combines phenotypic state counts (P), percent death (d) into a single matrix.
%Assumes phenotypic states: green (row 1 of P), *not* green (row 2 of P).
%All death in 'allDeath' phenotype. 
%No death in 'noDeath' phenotype.
%1 = green, 2 = *not* green.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.

function M = GetDataMatrix_UNEqualDeath(P, d, allDeath, noDeath)

[NumPheno,TWTotal] = size(P); %NumPheno = number of phenotypes, TWTotal = #Time points * #Wells = Number of samples.

n = NumPheno + 1; %Number of states.

M = zeros(n, TWTotal);

for j = 1:TWTotal
    
    M(3,j) = d(j)*sum(P(:,j)); % # dead cells = death percentage * population total
    
    M(allDeath, j) = P(allDeath, j) - M(3,j); % # 'allDeath' phenotype (alive) = # 'allDeath' phenotype counts - # dead cells
    
    M(noDeath, j) = P(noDeath, j); % # 'noDeath' phenotype (alive) = # 'noDeath' phenotype counts
    
end

%Verified 11/4/2016
