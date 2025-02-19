% This script generate figure 2

clear
close all
%% Make Figure
myFig = figure('Position',[300,300,1200,480]);
myColorOrder = colororder;
myColorOrder(2,:) = myColorOrder(2,:)*1.15;
colororder(myColorOrder)

%% Define Useful Functions
f  = @(P,L,m) P.n*(-P.ar*(m-P.mr0) + log((1+L./abs(P.Ki))./(1+L/abs(P.Ka))));
p  = @(P,L,m) 1./(1+exp(f(P,L,m)));
g  = @(P,m) P.as*(m-P.ms0);
c  = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));
ActSS = @(P,L,m) (1 - c(P,L,m) + 2*P.k/(1+exp(g(P,m))) + c(P,L,m)*2*P.k/(1+exp(-g(P,m))) ...
    - sqrt((-1 + c(P,L,m) - 2*P.k/(1+exp(g(P,m))) - c(P,L,m)*2*P.k/(1+exp(-g(P,m)))).^2 - 4*(1-c(P,L,m)).*c(P,L,m)*2*P.k/(1+exp(-g(P,m)))))...
    ./(2*(1-c(P,L,m)));

h = @(P,m) g(P,m) - P.n*(-P.ar*(m-P.mr0));
Q = @(P,m,c) (c.*exp(-h(P,m))).^(1/P.n);
LTocFunc = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));
cToLFunc = @(P,m,c) (Q(P,m,c)-1) ./ (1/P.Ka - Q(P,m,c)/P.Ki);

%% 2A: Adapt Time Series Panel
% load p10b_Fig2AdaptPanelDataLong4.mat AsFinal tsFinal MsFinal% load data
load Fig2DataForAdaptTimeSeries.mat AsFinal tsFinal MsFinal% load data


% select portion of time series to plot
As = AsFinal(4000:14000);
Ms = MsFinal(4000:14000);
ts = tsFinal(4000:14000)-tsFinal(4000);

axTimeSeries = axes(myFig); % create axis
% yyaxis left
plot(ts,As,'LineWidth',2.5) % plot data

% make data pretty
xlabel('Time (s)')
ylabel('Mean Activity ({\itA})')
xlim([min(ts),max(ts)])
ylim([0.1,0.7])
% ylim([0,1])

% add annotations for location of changing ligand
arrow1 = annotation("textarrow");
set(arrow1,'Parent',gca,'Position',[10,0.55,0,-0.1],'String',{'Ligand','Added'})
arrow2 = annotation('textarrow');
set(arrow2,'Parent',gca,'Position',[60,0.25,0,0.1],'String',{'Ligand','Removed'})
set(gca,'OuterPosition',[0,0,0.5,.9])

%% 2B: adaptation tuning
load Fig1-3_constrainedShimizuParams.mat

cArrayBase = 10.^(-3:0.00001:3); % array of c for plotting mf
M = fitParamsStruct.ms0; % methylation used to convert between c and L
LArrayBase = cToLFunc(fitParamsStruct,M,cArrayBase); % convert c to L for plotting mf data

% remove values of c and L corresponding to non-realizable values
LLim1 = find(LArrayBase > 0,1);
LLim2 = find(LArrayBase(LLim1+1:end) < 0,1) + LLim1 - 1;
if isempty(LLim2)
    LArrayBase = LArrayBase(LLim1:end);
    cArrayBase = cArrayBase(LLim1:end);
else
    LArrayBase = LArrayBase(LLim1:LLim2);
    cArrayBase = cArrayBase(LLim1:LLim2);
end

load('Fig1-2_stochasticDRC.mat') % load data for plotting

L = cToLFunc(fitParamsStruct,M,c); %  convert c to L for plotting simulation data

for i = 1:length(kArray)
    param{i} = fitParamsStruct;
    param{i}.k = kArray(i);
end

xMin = 10^-2;
xMax = 3;

cVals = [ceilS(LTocFunc(fitParamsStruct,xMin,M),1), round(LTocFunc(fitParamsStruct,10^-1.5,M),1,'significant'),1, ...
    round(LTocFunc(fitParamsStruct,xMax/4,M),1,'significant'),floorS(LTocFunc(fitParamsStruct,xMax,M),1)]; 
LVals = cToLFunc(fitParamsStruct,M,cVals);
for i = 1:length(cVals)
    cStrings{i} = num2str(cVals(i));
end



useLigand = 1;

kIndex = 3; %set k value for drc in plot: kIndex = 3 -> k = 0.01
ax1Adapt = axes(myFig);% make axis

% plot dose response curves
semilogx(ax1Adapt,L(kIndex,:),drc(kIndex,:),'lineWidth',2.5) % plot simulation data
set(gca,'ColorOrderIndex',1)
hold on
semilogx(ax1Adapt,LArrayBase,ActSS(param{kIndex},LArrayBase,M),'--','linewidth',2) % plot mf data
semilogx(ax1Adapt,LArrayBase,0.4*ones(length(LArrayBase),1),'k--','linewidth',2.5) % plot adaptation target
text(1.5*10^-2,0.5,{'Adaptation', 'Target'},'Color','k') % label adaptation target
xlim([xMin,xMax]) % set xaxis limits

% get indices corresponding to 90% and 10% activity
topThreshold = 0.8;
bottomThreshold = 0.2;
topIndex = find(drc(kIndex,:) > topThreshold, 1);
bottomIndex = find(drc(kIndex,:) > bottomThreshold, 1);

% get axis limits
xlims = xlim(ax1Adapt);
ylims = ylim(ax1Adapt);

% get ligand, activity, and control parameter values corresponding to top
% and bottom activity threshholds
x1 = L(kIndex,topIndex);
x2 = L(kIndex,bottomIndex);
y1 = bottomThreshold;
y2 = topThreshold;
c1 = LTocFunc(fitParamsStruct,x1,M);
c2 = LTocFunc(fitParamsStruct,x2,M);

% create shaded boxes
fill([xlims(1),xlims(1),xlims(2),xlims(2)],[y1,y2,y2,y1],[0,1,0]/2,'FaceAlpha',0.1,'EdgeColor','none') % wide tuning in Activity space
fill([x1,x1,x2,x2],[ylims(1),ylims(2),ylims(2),ylims(1)],[1,0,0]/2,'FaceAlpha',0.2,'EdgeColor','none') % narrow tuning in ligand space
fill([x1,x1,x2,x2],[y1,y2,y2,y1],[1,1,0],'FaceAlpha',0.4,'EdgeColor','none') % overlap between areas

set(gca,'OuterPosition',[0.5,0,0.5,.9],'XScale','log') % set plot size and location

% create axis labels
xlabel('Ligand Concentration [MeAsp] (mM)')
ylabel('Mean Activity ({\itA})')

fontsize(myFig,23,'pixels') % set font sizes

% convert ligand and activity levels to figure position units
[fx1,fy1] = axxy2figxy(x1,y1);
[fx2,fy2] = axxy2figxy(x2,y2);

% set colors for annotations
green = [34,139,34]/256/1.7;
red = [178,24,43]/256;

% create annotation for activity tuning
verticalArrow = annotation('doublearrow');
set(verticalArrow,'Position',[fx1+0.02,fy1,0,fy2-fy1],...
    'LineWidth',2,'Head1Width',15,'Head1Length',10,'Head2Width',15,'Head2Length',10,'Color',green)
text(0.2,0.5,{'Tune Activity',['within ' num2str(100*(topThreshold-bottomThreshold)) '%']},'Color',green)

% create annotation for ligand/control parameter tuning
horizontalArrowLeft = annotation('arrow');
horizontalArrowRight = annotation('arrow');
set(horizontalArrowLeft,'Position',[fx1-0.03,fy2+0.01,0.03,0],...
    'LineWidth',2,'HeadWidth',15,'HeadLength',10,'Color',red)
set(horizontalArrowRight,'Position',[fx2+0.03,fy2+0.01,-0.03,0],...
    'LineWidth',2,'HeadWidth',15,'HeadLength',10,'Color',red)
cTuningWidth = round(100*abs((c1-c2)));
text(0.3,topThreshold+.00,{'Tune {\itc}',['within ' num2str(cTuningWidth) '%']},'Color',red)


% Add Control parameter axis  and label
ax1Adapt.Box = 'off';
ax2Adapt = axes(myFig);
ax2Adapt.XAxisLocation = 'top';
ax2Adapt.YAxisLocation = 'right';
set(ax2Adapt,'XTick',LVals,'XTickLabel',cStrings,'xscale','log')
xlim([xMin,xMax])
xlabel('Control Parameter, {\itc}')
ax2Adapt.Color = 'none';
ax2Adapt.YTickLabel = [];


%% Final touches
% fontsize(myFig,23,'pixels') % set font sizes

textH = findall(gcf,'Type','text');
set(textH,'FontSize',23)

legendH = findall(gcf,'Type','legend');
set(legendH,'FontSize',20)
titleH = findall(gcf,'Type','Title');
set(titleH,'FontSize',20)
axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2,'FontSize',20,'LabelFontSizeMultiplier',1.2);

% axesH = findall(gcf,'Type','axes');
% set(axesH,'LineWidth',2);

% make sure all double x-axes align
if useLigand == 1
    set(ax2Adapt,'XMinorTick','off','Position',ax1Adapt.Position)
end
% set(ax2DRC,'XMinorTick','off','Position',ax1DRC.Position)


% make figure labels
figLabelSize = 30;
annotation('textbox',[0 1 0 0], ...
    'String','A','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[0.5 1 0 0], ...
    'String','B','EdgeColor','none','FontSize',figLabelSize)

% saveas(gcf,'Plots/p10b_Fig2Complete.png') % save plot
% exportgraphics(gcf,'Plots/Fig2Complete.pdf','ContentType','vector')

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