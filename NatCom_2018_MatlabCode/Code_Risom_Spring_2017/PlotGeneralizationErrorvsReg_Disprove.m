%Plots generalization error versus regularization parameter, lambda, to verify the choice of regularizer range.
%Optimal lambda should be a minimum argument of generalization error.

function PlotGeneralizationErrorvsReg_Disprove(e_lambda_Control, indexStar_Control, e_lambda_GSK, indexStar_GSK, e_lambda_BEZ, indexStar_BEZ, SetLambda)

NConditions = 3;

figure

FigureSettings;

for c = 1 : NConditions
    
    subplot(1, NConditions, c)
    
    if c == 1
        e_lambda = e_lambda_Control;
        indexStar = indexStar_Control;
    elseif c == 2
        e_lambda = e_lambda_GSK;
        indexStar = indexStar_GSK;
    else
        e_lambda = e_lambda_BEZ;
        indexStar = indexStar_BEZ;
    end
    
    semilogx(SetLambda, e_lambda); hold on
    semilogx(SetLambda(indexStar), e_lambda(indexStar),'*r');
    
    if c == 1
        title('DMSO');
        ylabel('Generalization error');
    elseif c == 2
        title('Trametinib');    
    else
        title('BEZ235');
        legend('e_\lambda','(\lambda^*, e_\lambda^*)');
    end
    xlabel('\lambda');
    
end
        