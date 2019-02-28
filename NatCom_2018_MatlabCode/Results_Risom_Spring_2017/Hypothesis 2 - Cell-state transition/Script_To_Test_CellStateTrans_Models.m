%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Margaret P. Chapman
% DATE: September 20, 2017
% PURPOSE: Test dynamics matrices identified on 15-well data under cell state transition hypothesis
% NOTE: Test data is 4-well data collected on different day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clearvars; clc;

%% Declare constants

SHEET_NAME = 'Summary';                         % missing data marked with 'na' for 'not available'
COL_ROW_1 = 'E5';                               % 1st col-row
N_AGENT = 4;                                    % # agents in data spreadsheet
N_T = 6;                                        % # time points in 4-well data set (0 h, 12 h, ..., 60 h)
EVEN_DEATH = 2;                                 % death distributed evenly across phenotypic states
N_W = 4;                                        % # wells in 4-well data set

%% Get 4-well phenotype & death timeseries raw data

[ P4W, d4W ] = SortBySeries( SortRawData( 'Timeseries_Raw_4wells.xlsx', SHEET_NAME, COL_ROW_1, 'K124', N_AGENT, N_T, N_W ) );

%% Get dynamics identified on 15-well data under cell state transition hypothesis

load('CellStateTransition_allResults_May9.mat'); % -> AStar_BEZ, AStar_DMSO, AStar_GSK from 15-well data

% Put dynamics into cell (1: DMSO, 2: GSK, 3: BEZ)
A_STAR = cell( N_AGENT-1, 1 ); A_STAR{1} = AStar_DMSO; A_STAR{2} = AStar_GSK; A_STAR{3} = AStar_BEZ; 

%% Compute test and predicted trajectories

M = cell( N_AGENT-1, 1 ); M_HAT = cell( N_AGENT-1, 1 );

for a = 1 : N_AGENT-1 % COMBO dynamics were not identified under cell state transition hypothesis
    
    % format raw phenotype timeseries : K14hi, K14low
    P4W_a = GetK14( P4W{a} );
    
    % 4-well data, even death
    M{a} = GetCellCtsMATRIX( GetCellCtsCELL( P4W_a, d4W{a}, EVEN_DEATH, 1:N_W ), 1:N_W ); % all 4 wells
    
    % prediction using A*(15-well), ( 12h, ..., 72h )
    M_HAT{a} = A_STAR{a} * M{a}; 
               
end

%% Compare test vs. predicted trajectories

SeeDynFit( M, M_HAT, N_W, {'DMSO', 'Tram.', 'BEZ'}, {'K14+ live', 'K14- live', 'dying/dead'} );


