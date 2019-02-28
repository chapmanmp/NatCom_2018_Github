%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Margaret P. Chapman
%Date: May 3, 2017

%PURPOSE: Propogates linear time-invariant dynamics, x_{k+1} = A x_k, forward in time
%INPUT:
    %A = dynamics matrix
    %x0 = initial condition
    %T = length of simulation horizon
%OUTPUT:
    %traj_forward(:, 1) = x0, traj_forward(:, 2) = x1, ..., traj_forward(:, T) = x_{T-1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajectory = SimulateLTIdynamics( A, x0, T )

n = length( x0 ); %number of states

trajectory = zeros( n, T );

trajectory( :, 1 ) = x0;

for k = 1 : T - 1,
    
    trajectory( :, k+1 ) = A * trajectory( :, k );
    
end
