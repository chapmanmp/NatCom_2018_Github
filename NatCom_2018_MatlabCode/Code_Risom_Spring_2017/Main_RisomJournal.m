%AUTHOR: Margaret P. Chapman
%DATE: May 9, 2017

%PURPOSE: Risom journal paper. Tests if: 
    %(1) K14high Darwinian selection or cell-state transition is the stronger driver of Trametinib-induced differentiation state-enrichment
    %(2) K14low Darwinian selection or cell-state transition is the stronger driver of BEZ235-induced differentiation state-enrichment
    
%NOTE: 
    %Dynamical system model, x_{k+1} = A x_k, identified via 15-well time series data
    %x = (# K14high live cells, # K14low live cells, # dead cells)'
    %Does not include statistical analysis.
    %HYPOTHESIS = 1, Darwinian selection
    %HYPOTHESIS = 2, Cell-state transition 
    %Ensure that path variable is properly defined to import data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. SPECIFY HYPOTHESIS

close all; clear all; clc;

HYPOTHESIS = 2;

if HYPOTHESIS == 1 %Darwinian selection
    equal = 0;                                  %Training data: death NOT distributed equally (GSK, BEZ)
    AllDeathNoDeathGains_GSKBEZ = 1;            %Model ID: all death in one gain, no death in other gain (GSK, BEZ)
    EqualDeathGains_DMSO = 1;                   %Model ID: equal death gains (DMSO)
    switching = 0;                              %Model ID: NO switching (GSK, BEZ, DMSO)
    
elseif HYPOTHESIS == 2 %Cell-state transition
    equal = 1;                                  %Training data: death distributed equally (GSK, BEZ)
    AllDeathNoDeathGains_GSKBEZ = 0;            %Model ID: no special death gain constraints (GSK, BEZ)
    EqualDeathGains_DMSO = 0;                   %Model ID: no special death gain constraints (DMSO)
    switching = 1;                              %Model ID: YES switching (GSK, BEZ, DMSO)
end

%Function arguments for unequal death distribution (GSK, BEZ)
allDeath_GSK = 2; %'allDeath' phenotype = 2, all death in K14low (GSK)
allDeath_BEZ = 1; %'allDeath' phenotype = 1, all death in K14high (BEZ)
%% 1. CONSTRUCTION OF TRAINING DATA MATRICES

path = 'C:\Users\chapmanm\Documents\MATLAB\BiologicalNetworks\Phenotypic State Modeling\D15_Summer2016to2017\Risom_Spring_2017\Code_Risom_Spring_2017\Phenotype, Death Assay Raw Counts Wells Matched, June 16 2016, TR.xlsx';
%Missing data replaced with 0s in spreadsheet. Combination data corrected.
PhenotypeObs_Col1 = 'E'; PhenotypeObs_ColLast = 'H'; DeathObs_TotalCol = 'J'; DeathObs_DeadCol = 'K'; T = 7; NumConditions = 4; WTotal = 15;

%Import experimentally observed phenotypic state counts (P) and death fraction (d) (DMSO, GSK, BEZ)
[P_DMSO, d_DMSO, P_GSK, d_GSK, P_BEZ, d_BEZ, ~, ~] = GetRawData(path, T, WTotal, NumConditions, PhenotypeObs_Col1, PhenotypeObs_ColLast, DeathObs_TotalCol, DeathObs_DeadCol);

%Get phenotype decomposition: K14hi, K14low
P_DMSO = Get_GNG_Decomposition(P_DMSO); %row 1 : K14hi, row 2 : K14low 
P_GSK = Get_GNG_Decomposition(P_GSK); 
P_BEZ = Get_GNG_Decomposition(P_BEZ);

%Get training data matrix for specified death distribution
M_DMSO = GetDataMatrix( P_DMSO, d_DMSO );                         %death is distributed equally (DMSO)
if equal                                                          %if death is distributed equally (GSK, BEZ)
    M_GSK = GetDataMatrix( P_GSK, d_GSK );
    M_BEZ = GetDataMatrix( P_BEZ, d_BEZ );
else                                                              %if death is NOT distributed equally (GSK, BEZ)
    M_GSK = GetDataMatrix_UNEqualDeath( P_GSK, d_GSK, allDeath_GSK, 1 ); %put all death in K14low (GSK)    
    M_BEZ = GetDataMatrix_UNEqualDeath( P_BEZ, d_BEZ, allDeath_BEZ, 2 ); %put all death in K14hi (BEZ)    
end

%% 2. MODEL IDENTIFICATION: x = (# K14hi live, # K14low live, # dead)'

%Identify wells for which no observations have been discarded for DMSO, GSK, BEZ for all time
SetWExists = GetCompleteWells_Disprove(M_DMSO, M_GSK, M_BEZ, WTotal); 

%Set regularizer range
NSetLambda = 6; SetLambda = logspace(0, 3, NSetLambda);

%Set stopping criterion for alternating minimization
epsilon = 10^(-6);

%ContraintsOnA_Disprove.m: turns switching on or off, applies special death gain constraints if specified 
%Optimize dynamics for specified hypothesis
[AStar_DMSO, ~, e_lambda_DMSO, indexStar_DMSO] = IdentifyModel_Disprove(M_DMSO, SetWExists, T, SetLambda, epsilon, 0, switching, 0, EqualDeathGains_DMSO);
[AStar_GSK, ~, e_lambda_GSK, indexStar_GSK] = IdentifyModel_Disprove(M_GSK, SetWExists, T, SetLambda, epsilon, 1, switching, AllDeathNoDeathGains_GSKBEZ, 0);
[AStar_BEZ, ~, e_lambda_BEZ, indexStar_BEZ] = IdentifyModel_Disprove(M_BEZ, SetWExists, T, SetLambda, epsilon, 2, switching, AllDeathNoDeathGains_GSKBEZ, 0);

%Plot generalization error versus regularization parameter, lambda.
PlotGeneralizationErrorvsReg_Disprove(e_lambda_DMSO, indexStar_DMSO, e_lambda_GSK, indexStar_GSK, e_lambda_BEZ, indexStar_BEZ, SetLambda);

%% 3. ANALYSIS - Plot fold change

x0_DMSO_avg = mean( M_DMSO( :, 1:WTotal ), 2 ); %average state at time 0 under DMSO

%Get fold change using SIMULATED values of #K14high (live & dead), #K14low (live & dead)
%GSK vs. DMSO
[FoldChange_K14hi_GSKvsDMSO_AllWells, FoldChange_K14lo_GSKvsDMSO_AllWells] = SimulateFoldChange_ManyLines(x0_DMSO_avg, AStar_DMSO, M_GSK( :, 1:WTotal ), AStar_GSK, T, equal, allDeath_GSK);
K14hi_GSKvsDMSO_Time72_Avg = mean( FoldChange_K14hi_GSKvsDMSO_AllWells(end,:) );
K14lo_GSKvsDMSO_Time72_Avg = mean( FoldChange_K14lo_GSKvsDMSO_AllWells(end,:) );

%BEZ vs. DMSO
[FoldChange_K14hi_BEZvsDMSO_AllWells, FoldChange_K14lo_BEZvsDMSO_AllWells] = SimulateFoldChange_ManyLines(x0_DMSO_avg, AStar_DMSO, M_BEZ( :, 1:WTotal ), AStar_BEZ, T, equal, allDeath_BEZ);
K14hi_BEZvsDMSO_Time72_Avg = mean( FoldChange_K14hi_BEZvsDMSO_AllWells(end,:) );
K14lo_BEZvsDMSO_Time72_Avg = mean( FoldChange_K14lo_BEZvsDMSO_AllWells(end,:) );

%Compare simulated vs observed fold change.
RealFoldChangeObservations; %Get fold change using EXPERIMENTALLY OBSERVED values of #K14high (live & dead), #K14low (live & dead)

plotFoldChange(FoldChange_K14hi_GSKvsDMSO_AllWells, FoldChange_K14lo_GSKvsDMSO_AllWells, FoldChange_K14high_GSKvsControl_Real, STD_K14high_GSKvsControl_Real, ...
                                                                                                  FoldChange_K14low_GSKvsControl_Real, STD_K14low_GSKvsControl_Real, 'Trametinib vs. DMSO', 4.6); 

plotFoldChange(FoldChange_K14hi_BEZvsDMSO_AllWells, FoldChange_K14lo_BEZvsDMSO_AllWells, FoldChange_K14high_BEZvsControl_Real, STD_K14high_BEZvsControl_Real, ...
                                                                                                  FoldChange_K14low_BEZvsControl_Real, STD_K14low_BEZvsControl_Real, 'BEZ235 vs. DMSO', 2.0); 


