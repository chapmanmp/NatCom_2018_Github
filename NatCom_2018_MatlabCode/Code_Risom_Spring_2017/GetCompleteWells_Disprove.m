%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identifies indices of wells that have existing samples for Control, GSK, and BEZ for all time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SetWExists = GetCompleteWells_Disprove(M_Control, M_GSK, M_BEZ, WTotal)

J_Empty_Control = GetDiscardedSamples(M_Control);

J_Empty_GSK = GetDiscardedSamples(M_GSK);

J_Empty_BEZ = GetDiscardedSamples(M_BEZ);

SetWEmpty = union( mod( J_Empty_BEZ, WTotal ), union( mod( J_Empty_Control, WTotal ), mod( J_Empty_GSK, WTotal ) ) );
%Set of incomplete wells: at least one observation has been discarded for a condition, a time point.

for i = 1:length(SetWEmpty)
    
    if SetWEmpty(i) == 0, SetWEmpty(i) = WTotal; end
    
end

SetWExists = setdiff( 1:WTotal, SetWEmpty );

%Verified 11/4/2016
