%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Margaret P. Chapman
%Date: May 2, 2017

%PURPOSE: Computes # live + # dead cells for phenotype p assuming equal death proportion for K14hi and K14lo
%INPUT:
    %x = (#K14hi live cells, #K14lo live cells, #dead cells)^T
    %p = phenotype index : 1 = K14hi, 2 = K14lo
%OUTPUT: # phenotype p live + # phenotype p dead
%NOTE: For equal allocation of death,
    % dead fraction = (#dead cells)/(#K14hi live + #K14lo live + #dead)
    % [#K14hi live] = (1 - dead fraction) * [#K14hi]
    % [#K14lo live] = (1 - dead fraction) * [#K14lo]
    % [#K14hi] := [#K14hi live] / (1 - dead fraction)
    % [#K14lo] := [#K14lo live] / (1 - dead fraction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NumPhenotypeLivePlusDead = LivePhenotype_Plus_DeadPhenotype_UnderEqualDeath( x, p )

NumDead = x( end );
    
NumTotal = sum( x );
    
DeadFraction = NumDead / NumTotal;

NumPhenotypeLive = x( p );

NumPhenotypeLivePlusDead = NumPhenotypeLive / ( 1 - DeadFraction );
