%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Selects a regularizer via cross-validation for model identification.
%2 phenotypic states: green, *not* green. 
%Switching specified on or off
%Death gain constraints, or lack thereof, specified via AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e_lambda, indexStar] = GetBestRegularizer_Disprove(MTest, MTrain, T, SetLambda, epsilon, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO)

NSetLambda = length(SetLambda);

[XTest, XTestPlus] = ShiftTimeHorizon(MTest, T);

e_lambda = zeros(1, NSetLambda);

for i = 1:NSetLambda
    lambda = SetLambda(i);
    [Ahatlambda, ~] = GetDynamicsAndCompleteData_Disprove(MTrain, lambda, T, epsilon, noDeath, switching, AllDeathNoDeathGains_GSKBEZ, EqualDeathGains_DMSO);
    e_lambda(i) = norm(XTestPlus - Ahatlambda*XTest,'fro');
    pause(0.5); %Pauses execution for 0.5 sec so that Control-c can stop code, if needed.
end

%Find largest lambda with minimum error rounded to 1 decimal place.

ehat_lambda = roundn(e_lambda,-1); %Round error to tenth decimal place.
indexStar = find(ehat_lambda == min(ehat_lambda),1,'last'); %Locates highest index corresponding to minimum rounded error.
%lambdaStar = SetLambda(indexStar);

% %Test script: Is the highest index associated with minimum rounded error found?
% Errorhat = [10.5, 6.7, 11.3, 4.6, 5.7, 4.6, 7.2];
% %(Errorhat == min(ErrorHat)) should output [0 0 0 1 0 1 0]
% indexStar = find(Errorhat == min(Errorhat),1,'last');
% if (indexStar == 6), display('Correct index found.');
% else display ('Incorrect index found.');
% end

%Test passes 4/21






