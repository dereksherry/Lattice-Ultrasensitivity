%% This script makes figures helping to justify the approximation to the mf we use to show the diverging susceptibility.
clear
close all

%% Make Figure
myFig = figure('Position',[0,0,1200,960]);

myColorOrder = colororder;
myColorOrder(2,:) = myColorOrder(2,:)*1.15;
colororder(myColorOrder)

load ../p57_data.mat
mArray = mArray/2;

%% Useful Equations
f  = @(P,L,m) P.n*(-P.ar*(m-P.mr0) + log((1+L./abs(P.Ki))./(1+L/abs(P.Ka))));
p  = @(P,L,m) 1./(1+exp(f(P,L,m)));
g  = @(P,m) P.as*(m-P.ms0);
c  = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));
ActSS = @(P,L,m) (1 - c(P,L,m) + 2*P.k/(1+exp(g(P,m))) + c(P,L,m)*2*P.k/(1+exp(-g(P,m))) ...
    - sqrt((-1 + c(P,L,m) - 2*P.k/(1+exp(g(P,m))) - c(P,L,m)*2*P.k/(1+exp(-g(P,m)))).^2 - 4*(1-c(P,L,m)).*c(P,L,m)*2*P.k/(1+exp(-g(P,m)))))...
    ./(2*(1-c(P,L,m)));
%% Compare Gain of Approximate vs Exact MF
% for large k
axLargek = axes(myFig);
kIndex = 3;
plot(mArray,hillcoefsApprox(kIndex,:),mArray,hillcoefsFull(kIndex,:),'LineWidth',2)
hold on
plot([params.ms0,params.ms0]/2,get(gca,'YLim'),'k','LineWidth',3)
legend('MF Approx.','MF')

xlabel('Methylation Level')
ylabel('Hill Coefficient')
set(gca,'OuterPosition',[0.5,0.5,0.5,0.5])
text(1,3.5,'{\itk} = 0.1')

% for medium k
axLargek = axes(myFig);
kIndex = 2;
plot(mArray,hillcoefsApprox(kIndex,:),mArray,hillcoefsFull(kIndex,:),'LineWidth',2)
hold on
plot([params.ms0,params.ms0]/2,get(gca,'YLim'),'k','LineWidth',3)
legend('MF Approx.','MF')
text(1,20,'{\itk} = 0.01')

xlabel('Methylation Level')
ylabel('Hill Coefficient')
set(gca,'OuterPosition',[0,0,0.5,0.5])

% for small k
axLargek = axes(myFig);
kIndex = 1;
plot(mArray,hillcoefsApprox(kIndex,:),mArray,hillcoefsFull(kIndex,:),'LineWidth',2)
hold on
plot([params.ms0,params.ms0]/2,get(gca,'YLim'),'k','LineWidth',3)
legend('MF Approx.','MF')
text(1,200,'{\itk} = 0.001')


xlabel('Methylation Level')
ylabel('Hill Coefficient')
set(gca,'OuterPosition',[0.5,0,0.5,0.5])


%% Compare DRCs
load ../p57_DRCdata.mat
kIndex = 3;
axDRCCompare = axes(myFig);
hold on

set(gca,'OuterPosition',[0,0.5,0.5,0.5],'XScale','log')
myColor = [0,0,0];
plot(0,0,'k--',0,0,'k','lineWidth',2)
colororder(gca,'gem12')
set(gca,'ColorOrderIndex',1)
for i = 1:1:9
    plot(L,drcFull{i},'lineWidth',2)
end
set(gca,'ColorOrderIndex',1)
for i = 1:1:9
    plot(L,drcApprox{i},'--','lineWidth',2)
end
plot(L,ActSS(params,L,params.ms0),'k','LineWidth',4)
legend('MF Approx.','MF','location','west')
xlabel('Ligand Concentration')
ylabel('Activity')
xlim([10^-5,10^2])

%% Finishing Touches
textH = findall(gcf,'Type','text');
set(textH,'FontSize',23)

legendH = findall(gcf,'Type','legend');
set(legendH,'FontSize',20)
titleH = findall(gcf,'Type','Title');
set(titleH,'FontSize',20)
axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2,'FontSize',20,'LabelFontSizeMultiplier',1.2);


% Make Figure Labels
figLabelSize = 30;
annotation('textbox',[0 1 0 0], ...
    'String','A','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[0.5 1 0 0], ...
    'String','B','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[0 0.5 0 0], ...
    'String','C','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[0.5 0.5 0 0], ...
    'String','D','EdgeColor','none','FontSize',figLabelSize)
%% Save Figures
exportgraphics(gcf,'Plots/FigSIMFApprox.pdf','ContentType','vector')