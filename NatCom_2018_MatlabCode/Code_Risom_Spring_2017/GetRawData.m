%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PURPOSE: Import raw phenotypic state counts (P) and percent death (d) for each condition.

%INPUTS: path of spreadsheet, T = Number of time points, W = Number of wells, C = Number of conditions,
%Phenotype_Col1 = column assoc. with phenotype 1, Phenotype_ColLast = column assoc. with final phenotype, 
%DeathExpt_TotalCol = column assoc. with total counts from dead expt. DeathExpt_DeadCol = column assoc. with counts of dead cells

%OUTPUTS: P(i,:) = number of cells in phenotypic state i, d(j) = percent death for sample j. 
%Column j in...
                %1:W = well 1 to well W at time 0,
                %W+1:2W = well 1 to well W at time 1, ...
                %(T-2)W+1:(T-1)W = well 1 to well W at time T-1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_Control, d_Control, P_GSK, d_GSK, P_BEZ, d_BEZ, P_COMBO, d_COMBO] = GetRawData(path, T, W, C, Phenotype_Col1, Phenotype_ColLast, DeathExpt_TotalCol, DeathExpt_DeadCol)

Control_PhenotypeExpt = []; GSK_PhenotypeExpt = []; BEZ_PhenotypeExpt = []; COMBO_PhenotypeExpt = [];
Control_DeadExpt = []; GSK_DeadExpt = []; BEZ_DeadExpt = []; COMBO_DeadExpt = [];

LastRow = 4; %Initial condition, row number of "Phenotype" "Death" titles.

for time = 1:T %0 hr to 72 hr
    for condition = 1:C 
        FirstRow = LastRow + 2; LastRow = FirstRow + W - 1;
        NextPhenotypeData = xlsread(path,'Summary',[Phenotype_Col1, num2str(FirstRow),':', Phenotype_ColLast, num2str(LastRow)]);
        NextDeadRawData = xlsread(path,'Summary',[DeathExpt_TotalCol, num2str(FirstRow),':', DeathExpt_DeadCol, num2str(LastRow)]);
        NextDeadPercentage = NextDeadRawData(:,2)./NextDeadRawData(:,1);
        if condition == 1
            Control_PhenotypeExpt = [Control_PhenotypeExpt; NextPhenotypeData];
            Control_DeadExpt = [Control_DeadExpt; NextDeadPercentage];
        elseif condition == 2
            GSK_PhenotypeExpt = [GSK_PhenotypeExpt; NextPhenotypeData];
            GSK_DeadExpt = [GSK_DeadExpt; NextDeadPercentage]; 
        elseif condition == 3
            BEZ_PhenotypeExpt = [BEZ_PhenotypeExpt; NextPhenotypeData];
            BEZ_DeadExpt = [BEZ_DeadExpt; NextDeadPercentage];
        elseif condition == 4
            COMBO_PhenotypeExpt = [COMBO_PhenotypeExpt; NextPhenotypeData];
            COMBO_DeadExpt = [COMBO_DeadExpt; NextDeadPercentage];
        end
    end
end

%Each column is associated with a time point and well.
P_Control = Control_PhenotypeExpt';
d_Control = Control_DeadExpt';
P_GSK =  GSK_PhenotypeExpt';
d_GSK = GSK_DeadExpt';
P_BEZ = BEZ_PhenotypeExpt';
d_BEZ = BEZ_DeadExpt';
P_COMBO = COMBO_PhenotypeExpt';
d_COMBO = COMBO_DeadExpt';

%VERIFICATION
%This function was tested during the week of 4/12, 6/2.





                