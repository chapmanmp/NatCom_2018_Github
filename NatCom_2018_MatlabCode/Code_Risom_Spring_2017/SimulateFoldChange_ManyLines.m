%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Margaret P. Chapman
%Date: May 3, 2017

%PURPOSE: To compute fold change for K14hi & K14lo drug vs. DMSO via dynamics simulated from many replicate initial conditions 
%INPUT:
    %x0_DMSO = initial state, DMSO
    %A_DMSO = LTI dynamics matrix under DMSO
    %x0_AllWells_drug(:, w) = initial state of drug-treated well w
    %A_drug = LTI drug-treated dynamics matrix
    %T = length of simulation horizon
    %equal = designates how death is split between K14hi & K14lo following drug treatment
    %allDeath = index of phenotype that absorbs all death; relevant only if death is split unequally
%OUTPUT: 
    %foldChange_ phenotype p _drugvsDMSO_AllWells( :, w ) = fold change, phenotype p, drug-treated dynamics initialized in well w, DMSO-input dynamics initialized at x0_DMSO
%NOTE: 
    %state vector x = [#K14hi live; #K14lo live; #dead]
    %fold change phenotype p, time k = [fraction phenotype p, time k, drug] / [fraction phenotype p, time k, DMSO]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [foldChange_K14hi_drugvsDMSO_AllWells, foldChange_K14lo_drugvsDMSO_AllWells] = SimulateFoldChange_ManyLines(x0_DMSO, A_DMSO, x0_drug_AllWells, A_drug, T, equal, allDeath)

[ ~, WTotal ] = size( x0_drug_AllWells );

foldChange_K14hi_drugvsDMSO_AllWells = zeros( T, WTotal ); 

foldChange_K14lo_drugvsDMSO_AllWells = zeros( T, WTotal );

%simulate DMSO-input dynamics, initialized at x0_DMSO
traj_DMSO = SimulateLTIdynamics( A_DMSO, x0_DMSO, T );

%for each drug-treated well
for w = 1 : WTotal
    
    %simulate drug-treated dynamics, initialized in well w
    traj_drug_w = SimulateLTIdynamics( A_drug, x0_drug_AllWells( :, w ), T );

    %get fold change K14hi drug (well w) vs. DMSO
    foldChange_K14hi_drugvsDMSO_AllWells( :, w ) = ComputeFoldChange(traj_DMSO, traj_drug_w, 1, equal, allDeath); %phenotype index : 1 = K14hi
    
    %get fold change K14lo drug (well w) vs. DMSO
    foldChange_K14lo_drugvsDMSO_AllWells( :, w ) = ComputeFoldChange(traj_DMSO, traj_drug_w, 2, equal, allDeath); %phenotype index : 2 = K14lo
    
end


 

    
    
    
    
    



