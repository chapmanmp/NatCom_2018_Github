function plotFoldChange(FoldChange_K14high_TherapyvsControl_AllWells, FoldChange_K14low_TherapyvsControl_AllWells, FoldChange_K14high_TherapyvsControl_Real, STD_K14high_TherapyvsControl_Real, ...
                                                                                                                   FoldChange_K14low_TherapyvsControl_Real, STD_K14low_TherapyvsControl_Real, plotTitle, YMAX) 
[T, WTotal] = size(FoldChange_K14high_TherapyvsControl_AllWells);

figure
FigureSettings

TimeHorizon = 0: 12: (T-1)*12; 

%Simulated fold change (#K14high living, #K14low living)
for w = 1 : WTotal
    
    plot(TimeHorizon, FoldChange_K14high_TherapyvsControl_AllWells(:, w), 'g'); hold on 

    plot(TimeHorizon, FoldChange_K14low_TherapyvsControl_AllWells(:, w), 'Color', [.6 .6 .6]); hold on %[red, green, blue]
      
end

%Observed fold change (#K14high living and dead, #K14low living and dead)
errorbar(TimeHorizon, FoldChange_K14high_TherapyvsControl_Real, STD_K14high_TherapyvsControl_Real, 'Color', [0 .5 0]); hold on

errorbar(TimeHorizon, FoldChange_K14low_TherapyvsControl_Real, STD_K14low_TherapyvsControl_Real, 'k'); 

legend('K14 high', 'K14 low')
title(plotTitle);
xlabel('Time (hr)')
ylabel('Simulated fold change');
axis([0 72 0 YMAX]);


