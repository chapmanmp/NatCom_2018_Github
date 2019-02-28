%Changes phenotype distribution green(1)/red(2)/blue(3)/not(4) to green(1)/not green(2).
%P(i,j) = # cells in phenotype i from sample j

function P_GnotG = Get_GNG_Decomposition(P)

Green = P(1, :); %row 1

NotGreen = sum( P( 2:end, :) ); %row 2 + row 3 + row 4

P_GnotG = [Green; NotGreen];