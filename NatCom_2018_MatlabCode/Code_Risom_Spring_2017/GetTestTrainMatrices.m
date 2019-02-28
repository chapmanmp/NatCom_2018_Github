%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Splits original data into test and train matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MTest, MTrain] = GetTestTrainMatrices(M,SetWExists,WTotal)

[n,TWTotal] = size(M); T = TWTotal/WTotal;

WTest = length(SetWExists); %Number of test wells
WTrain = WTotal - WTest; %Number of train wells

MTest = zeros(n,T*WTest); MTrain = zeros(n,T*WTrain);

for k = 0:T-1
    wTest = 0; 
    wTrain = 0;
    for w = 1:WTotal
        if (ismember(w,SetWExists))
            wTest = wTest + 1;
            MTest(:,wTest + k*WTest) = M(:,w + k*WTotal);
        else
            wTrain = wTrain + 1;
            MTrain(:,wTrain + k*WTrain) = M(:,w + k*WTotal);
        end
    end
end

%VERIFICATION
%This function was tested during the week of 4/12.




