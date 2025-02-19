clear
close all
%% Make Figure
myFig = figure('Position',[0,0,1800,960]);
myColorOrder = colororder;
myColorOrder(2,:) = myColorOrder(2,:)*1.15;
colororder(myColorOrder)

%% Load Data
load('Fig3DataA.mat')
load('Fig3DataB.mat')
load('Fig3DataC.mat')
load('Fig3DataD.mat')
load('Fig3DataE.mat')
load('Fig3DataF.mat')


% load('defaultParams.mat')

%% 3A: Plot PSD vs k
lorentzFunction = @(P,x) P(1)./(1+((x)/P(2)).^2); % define lorentz function
axPSDvsk = axes(myFig);
for i = 1:4
    loglog(axPSDvsk,Fig3DataA.freq,Fig3DataA.psdk{i},'linewidth',2) % plot psd from sim
    hold on
    loglog(axPSDvsk,Fig3DataA.freq,lorentzFunction(Fig3DataA.lorentzFitVals{i},Fig3DataA.freq),'k','linewidth',2) % plot lorentz fit of psd
end

xlabel('Frequency (1/s)')
ylabel('Power Spectral Density')
set(gca,'XScale','log','YScale','log','FontSize',16,'OuterPosition',[0,1/2,1/3,0.45])
legend('{\itk} = 0.0001','','{\itk} = 0.001','','{\itk} = 0.01','','{\itk} = 0.1','location','southwest')
text(0.1,8,'Mean Activity = 0.5','FontSize',14)

%% 3B: Plot stats vs k
axStatsVsk = axes(myFig);
yyaxis left
loglog(axStatsVsk,Fig3DataB.k,Fig3DataB.variance,'o-','linewidth',2)
xlabel('{\itk}')
ylabel('Activity Variance (<A^2> - <A>^2)')
text(10^-4,0.8*10^-2,'Mean Activity = 0.5','FontSize',14)

yyaxis right
loglog(axStatsVsk,Fig3DataB.k,Fig3DataB.lorentzFitVals,'o-','linewidth',2)
ylabel('Characteristic Frequency (1/s)','Rotation',-90)
xlim([Fig3DataB.k(1),Fig3DataB.k(end)])
set(axStatsVsk,'XTick',[10^-5,10^-4,10^-3,10^-2,10^-1,1],'OuterPosition',[0,0,1/3,0.45])
% fontsize(axStatsVsk,23,'pixels')

%% 3C: Plot Variance vs L
load Fig1-3_constrainedShimizuParams.mat
M = fitParamsStruct.ms0;
% defaultParams = fitParamsStruct;

f  = @(P,L,m) P.n*(-P.ar*(m-P.mr0) + log((1+L./abs(P.Ki))./(1+L/abs(P.Ka))));
p  = @(P,L,m) 1./(1+exp(f(P,L,m)));
g  = @(P,m) P.as*(m-P.ms0);
c  = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));

h = @(P,m) g(P,m) - P.n*(-P.ar*(m-P.mr0));
Q = @(P,m,c) (c.*exp(-h(P,m))).^(1/P.n);
LTocFunc = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));
cToLFunc = @(P,m,c) (Q(P,m,c)-1) ./ (1/P.Ka - Q(P,m,c)/P.Ki);

L = cToLFunc(fitParamsStruct,M,Fig3DataC.c);
varianceVals = Fig3DataC.variance;

% remove unrealizable L1 vals
LLim1 = find(L > 0,1);
LLim2 = LLim1 + find(L(LLim1+1:end) < 0,1)-1;
if isempty(LLim2)
    L = L(LLim1:end);
    varianceVals = varianceVals(LLim1:end);
else
    L = L(LLim1:LLim2);
    varianceVals = varianceVals(LLim1:LLim2);
end

% make plot
ax1VarVsL = axes(myFig);
semilogx(ax1VarVsL,L,varianceVals,'o-','linewidth',2)
xlabel('Ligand Concentration [MeAsp] (mM)')
ylabel('Activity Variance (<A^2> - <A>^2)')
text(0.01,10^-2,'{\itk} = 0.01','Fontsize',14)
set(gca,'XTick',[10^-3,10^-2,10^-1,1],'OuterPosition',[1/3,0.5,1/3,0.45])
xlim([min(L),max(L)])

ax1VarVsL.Box = 'off';

ax2VarVsL = axes(myFig);
ax2VarVsL.XAxisLocation = 'top';
ax2VarVsL.YAxisLocation = 'right';

% cVals = [10,3,1,0.3,0.2];
cVals = [ceilS(LTocFunc(fitParamsStruct,min(L),M),1), round(LTocFunc(fitParamsStruct,min(L)*8,M),1,'significant'),1, ...
    round(LTocFunc(fitParamsStruct,max(L)/6,M),1,'significant'),floorS(LTocFunc(fitParamsStruct,max(L),M),1)]; 
LVals = cToLFunc(fitParamsStruct,M,cVals);
% cStrings = {'10','3','1','0.3','0.2'};
for i = 1:length(cVals)
    cStrings{i} = num2str(cVals(i));
end

set(ax2VarVsL,'XTick',LVals,'XTickLabel',cStrings,'xscale','log')
xlim([min(L),max(L)])
xlabel('Control Parameter, {\itc}')
ax2VarVsL.Color = 'none';
ax2VarVsL.YTickLabel = [];
set(ax2VarVsL,'XMinorTick','off','Position',ax1VarVsL.Position)

%% 3D: Plot tuned, untuned, and adapted PSD
axTunedVsAdapt = axes(myFig);
hold on
set(gca,'ColorOrderIndex',6)
loglog(Fig3DataD.freq,Fig3DataD.untunedPSD,'linewidth',2)
hold on
loglog(Fig3DataD.freq,Fig3DataD.tunedPSD,'linewidth',2)
xlabel('Frequency (1/s)')
ylabel('Power Spectral Density')
set(axTunedVsAdapt,'fontsize',15,'OuterPosition',[1/3,0,1/3,0.45],'XScale','log','YScale','log','Box','on')
xlim([10^-4,10])
ylim([10^-5,10])
loglog(Fig3DataD.adaptingFreq,Fig3DataD.adaptingPSD,'linewidth',2)
legend(axTunedVsAdapt,'Untuned','Tuned','Adapting','location','northeast')
legend('boxoff')

%% 3E: Plot Adapting PSDs
axAdaptingPSD = axes(myFig);
loglog(Fig3DataD.adaptingFreq, Fig3DataE.slowAdaptPSD,...
    Fig3DataD.adaptingFreq,Fig3DataE.mediumAdaptPSD,...
    Fig3DataD.adaptingFreq,Fig3DataE.fastAdaptPSD,'linewidth',2)
xlabel('Frequency (1/s)')
ylabel('Power Spectral Density')
xlim([10^-4,10])
ylim([10^-5,10])
set(gca,'fontsize',15)
legend('{\itk_B/k} = 0.016','{\itk_B/k} = 0.16','{\itk_B/k} = 16','fontsize',11)
legend('boxoff')
set(gca,'OuterPosition',[2/3,0.5,1/3,0.45])

%% 3F: Plot representative time series
% Make Big Figure
axSeries(1) = axes(myFig,'OuterPosition',[2/3,0,1/3,0.45]);
% set(gca,'FontSize',14,'Box','off','XColor','none','YColor','none','Color','none')
set(gca,'FontSize',20,'YColor','none','Color','none')
xlabel('Time (s)','Fontsize',20,'Color',[.15 .15 .15])
ylabel('Mean Activity ({\itA})','Fontsize',20,'Color',[.15 .15 .15])
xlim([0,20])
boundaries = axSeries(1).Position;

% Make Smaller Plots
axSeries(2) = axes(myFig);
plot(Fig3DataF.slowAdapt_ts, Fig3DataF.slowAdapt_As,'linewidth',2);
set(gca,'fontsize',15)
xticklabels({})
xlim([500,520])
ylim([-0.1,1.1])
text(519,0.1,'{\itk_B/k} = 0.016','HorizontalAlignment','right')

axSeries(3) = axes(myFig);
plot(Fig3DataF.mediumAdapt_ts, Fig3DataF.mediumAdapt_As,'Color',[0.8500 0.3250 0.0980]*1.15,'linewidth',2);
set(gca,'fontsize',15)
xticklabels({})
xlim([500,520])
ylim([-0.1,1.1])
text(519,0.1,'{\itk_B/k} = 0.16','HorizontalAlignment','right')


axSeries(4) = axes(myFig);
plot(Fig3DataF.fastAdapt_ts-500, Fig3DataF.fastAdapt_As,'color',[0.9290 0.6940 0.1250],'linewidth',2);
set(gca,'fontsize',15)
xticklabels({})
xlim([0,20])
ylim([-0.1,1.1])
text(19,0.1,'{\itk_B/k} = 16','HorizontalAlignment','right')


% add text labels


for i = 2:4
    axSeries(i).Position(1) = boundaries(1);
    axSeries(i).Position(2) = boundaries(2) + boundaries(4) - boundaries(4)*(i-1)/3;
    axSeries(i).Position(3) = boundaries(3);
    axSeries(i).Position(4) = boundaries(4)/3.1;
end
%% Final Touches
% fontsize(myFig,23,'pixels')
% 
% axesH = findall(gcf,'Type','axes');
% set(axesH,'LineWidth',2);

textH = findall(gcf,'Type','text');
set(textH,'FontSize',23)

legendH = findall(gcf,'Type','legend');
set(legendH,'FontSize',20)
titleH = findall(gcf,'Type','Title');
set(titleH,'FontSize',20)
axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2,'FontSize',20,'LabelFontSizeMultiplier',1.2);


% Fix Axes
set(ax2VarVsL,'Position',ax1VarVsL.Position)

% add figure labels
figLabelSize = 30;
annotation('textbox',[0 1 0 0], ...
    'String','A','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[0 0.5 0 0], ...
    'String','B','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[1/3 1 0 0], ...
    'String','C','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[1/3 0.5 0 0], ...
    'String','D','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[2/3 1 0 0], ...
    'String','E','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[2/3 0.5 0 0], ...
    'String','F','EdgeColor','none','FontSize',figLabelSize)

% exportgraphics(gcf,'Plots/Fig3Complete.pdf','ContentType','vector')

%% Functions
function flooredVal = floorS(val, nS)
    pw = ceil(log10(val)); % Find order of magnitude of val.
    res = 10^(pw-nS); % Resolution to round to.
    flooredVal = ceil(val/res)*res; % < change floor() to ceil(), for ceiling equivalent.
end

function ceiledVal = ceilS(val, nS)
    pw = ceil(log10(val)); % Find order of magnitude of val.
    res = 10^(pw-nS); % Resolution to round to.
    ceiledVal = floor(val/res)*res; % < change floor() to ceil(), for ceiling equivalent.
end