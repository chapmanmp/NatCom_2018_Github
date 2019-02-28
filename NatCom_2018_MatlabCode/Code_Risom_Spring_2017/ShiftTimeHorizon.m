%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Forms data matrices of earlier and later observations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%'M','X','XPlus' contain data for a specific condition.

function [X, XPlus] = ShiftTimeHorizon(M,T)

[~,TW] = size(M);

W = TW/T;

X = M(:,1:(T-1)*W);

XPlus = M(:,W+1:TW);

%VERIFICATION
%This function was tested during the week of 4/12.

