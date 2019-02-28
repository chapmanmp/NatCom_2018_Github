%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Margaret P. Chapman
%Date: May 3, 2017

%PURPOSE: Computes fold change (drug vs. DMSO) of phenotype p over time horizon
%INPUT:
    %traj_DMSO(:, k) = [#K14hi live; #K14lo live; #dead] at time k, following DMSO input at time 0
    %traj_drug(:, k) = [#K14hi live; #K14lo live; #dead] at time k, following drug treatment at time 0
    %p = phenotype index for fold change computation
    %equal = designates how death is split between K14hi & K14lo under drug
        %equal = 1 if death should be split evenly between K14hi & K14lo
        %equal = 0 if death should be allocated to one phenotype only
    %allDeath = index of phenotype that absorbs all death; relevant only if equal = 0
%NOTE: 
    % phenotype index : 1 = K14hi, 2 = K14lo
    % #K14hi dead, #K14lo dead depend on how death is allocated to phenotypes under drug 
    % #K14hi = #K14hi live + #K14hi dead
    % #K14lo = #K14lo live + #K14lo dead
    % cell total = #K14hi live + #K14lo live + #dead
    % fraction K14hi = #K14hi / cell total
    % fraction K14lo = #K14lo / cell total
%OUTPUT: FoldChange_drugvsDMSO(k) = [fraction phenotype p, time k, drug] / [fraction phenotype p, time k, DMSO]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function foldChange_drugvsDMSO = ComputeFoldChange(traj_DMSO, traj_drug, p, equal, allDeath)

[ ~, T ] = size( traj_DMSO ); %# time points

foldChange_drugvsDMSO = zeros( T, 1 );

for k = 1 : T
    
    xk_DMSO = traj_DMSO( :, k );
    
    %split death evenly between K14hi and K14lo under DMSO
    livePhenotype_plus_deadPhenotype_DMSO = LivePhenotype_Plus_DeadPhenotype_UnderEqualDeath( xk_DMSO, p );
    
    fraction_phenotype_DMSO = livePhenotype_plus_deadPhenotype_DMSO / sum( xk_DMSO );

    xk_drug = traj_drug( :, k );
    
    %if death is split evenly between K14hi and K14lo under drug
    if equal, livePhenotype_plus_deadPhenotype_drug = LivePhenotype_Plus_DeadPhenotype_UnderEqualDeath( xk_drug, p );
    
    %if not, allocate all death into 'allDeath' phenotype
    else      livePhenotype_plus_deadPhenotype_drug = xk_drug( p ) + ( allDeath == p ) * xk_drug( end ); 
    end
    
    fraction_phenotype_drug = livePhenotype_plus_deadPhenotype_drug / sum( xk_drug );
            
    foldChange_drugvsDMSO(k) = fraction_phenotype_drug / fraction_phenotype_DMSO;
    
end
    
%Tests passed 5/3    
    