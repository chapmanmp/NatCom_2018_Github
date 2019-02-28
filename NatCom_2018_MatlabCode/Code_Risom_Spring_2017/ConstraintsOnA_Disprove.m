%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constraints on dynamics to ensure a realistic biological system.
%Used to identify a dynamics matrix with fixed predicted data.
%2 phenotypic states: green (1), *not* green (2). 
%Death gain constraints, or lack thereof, specified via AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO
%Switching specified on or off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AVar - dynamics matrix variable that is undergoing optimization. 
%n - number of states.
%AllDeathNoDeathGains_GSKBEZ = 1 -> all death in one gain, no death in other gain; noDeath - 'noDeath' phenotype (1 = green or 2 = *not* green).
%EqualDeathGains_DMSO = 1 -> equal death gains
%AllDeathNoDeathGains_GSKBEZ = 0 & EqualDeathGains_DMSO = 0 -> no special contraints on death gains
%switching: 1 = on, 0 = off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AVar(:,n) == [zeros(n-1,1); 1]; %Constraint on last col of dynamics matrix, n scalar equality contraints

ones(1,n)*AVar(:, 1:n-1) >= 1; %Cell division gains must be >= 1, n-1 scalar inequality constraints
       
diag( AVar( 1:n-1 , 1:n-1 ) ) >= 0; %first n-1 diagonal terms must be >= 0, n-1 scalar inequality constraints
              
for j = 1:n-1
    for i = 1:n
        if i ~= j
            AVar(i,j) >= 0;
            AVar(i,j) <= 1;
        end                                          %Nondiagonal terms must be between 0 and 1.
    end
    
    if j ~= 1                                       
        sum(AVar(:,1)) == sum(AVar(:,j));            %Set all cell division gains equivalent due to EdU+ analysis.
    end
end

if AllDeathNoDeathGains_GSKBEZ                      %If all death in one gain, no death in other gain (GSK, BEZ)
    AVar(n, noDeath) == 0;                              %enforce no death in 'noDeath' phenotype. set death gain of 'noDeath' = 0.
elseif EqualDeathGains_DMSO                         %If equal death gains (DMSO), set death gains equal.
    AVar(n, 1) == AVar(n, 2);
else                                                %If not, then no special constraints on death gains
end

%If NO switching, turn switching off
if ~switching
    AVar(1,2) == 0;
    AVar(2,1) == 0; 
end

% CVX User guide 2008, p. 18, Sec. 3.3 Constraints
% "These equality and inequality operators work for arrays. When both sides of
% the constraint are arrays of the same size, the constraint is imposed elementwise. For
% example, if a and b are m×n matrices, then a<=b is interpreted by cvx as mn (scalar)
% inequalities, i.e., each entry of a must be less than or equal to the corresponding entry
% of b. cvx also handles cases where one side is a scalar and the other is an array. This
% is interpreted as a constraint for each element of the array, with the (same) scalar
% appearing on the other side. As an example, if a is an m × n matrix, then a>=0 is
% interpreted as mn inequalities: each element of the matrix must be nonnegative."


 %GSK vs. Control: all death in *not* green, p_1D = 0, 'noDeath' phenotype = green = 1
 %BEZ vs. Control: all death in green, p_2D = 0, 'noDeath' phenotype = *not* green = 2
