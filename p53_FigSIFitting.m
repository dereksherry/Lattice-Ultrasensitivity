clear
close all
load ../p54_constrainedShimizuParams.mat

%% Load Data
data  = readmatrix('../../Code For LU Paper/Shimizu Raw','NumHeaderLines',2);
% data  = table2array(data);
xdata = data(:,1:2:17);
ydata = data(:,2:2:18);
ydata = ydata./max(ydata,[],'all');

colors = 'ymcbrgkbr';

% EEEEx     = xdata(:,1);
% wtx       = xdata(:,2);
% QEEEx     = xdata(:,3);
% QEQEx     = xdata(:,4);
% QEQQx     = xdata(:,5);
% QQQQx     = xdata(:,6);
% QEmQQx    = xdata(:,7);
% QEmQEmx   = xdata(:,8);
% EmEmEmEmx = xdata(:,9);
% EEEEy     = ydata(:,1);
% wty       = ydata(:,2);
% QEEEy     = ydata(:,3);
% QEQEy     = ydata(:,4);
% QEQQy     = ydata(:,5);
% QQQQy     = ydata(:,6);
% QEmQQy    = ydata(:,7);
% QEmQEmy   = ydata(:,8);
% EmEmEmEmy = ydata(:,9);

%% Nice Equations
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

%% Plot Stuff
figure('Position',[200,200,600,500])
colororder(gca,'gem12')
scatter(xdata(:,[1,3:9])*10^3,ydata(:,[1,3:9]),'filled','SizeData',100)
hold on
set(gca,'ColorOrderIndex',1,'XScale','log')
for i = 0:7
    plotHandle = plot(0,0,'o-','LineWidth',2);
    set(plotHandle, 'MarkerFaceColor', get(plotHandle,'Color')); 
end
set(gca,'ColorOrderIndex',1)

L = 10.^(-4:0.01:2);
for m = [0:6,8]
    semilogx(L*10^3,ActSS(fitParamsStruct,L,m),'LineWidth',2)
end

xlabel('[MeAsp] (\muM)')
ylabel('Activity')
legend('','','','','','','','','EEEE','QEEE','QEQE','QEQQ','QQQQ','QEmQQ','QEmQEm','EmEmEmEm', ...
    'NumColumns',4,'Location','northoutside');
xlim([10^-1,10^5])
box on

%% Final Stuff
textH = findall(gcf,'Type','text');
set(textH,'FontSize',23)

legendH = findall(gcf,'Type','legend');
set(legendH,'FontSize',20)
titleH = findall(gcf,'Type','Title');
set(titleH,'FontSize',20)
axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2,'FontSize',20,'LabelFontSizeMultiplier',1.2);

exportgraphics(gcf,'Plots/FigSFitting.pdf','ContentType','vector')
