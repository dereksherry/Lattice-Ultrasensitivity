% This script generates figure 1
clear
close all

%% Set up the figure
myFig = figure(1);
% myFig.Position(3) = 1800;
% myFig.Position(4) = 960;
myFig.Position(3) = 1600;
myFig.Position(4) = 960;
myColorOrder = colororder;
myColorOrder(2,:) = myColorOrder(2,:)*1.15;
colororder(myColorOrder);

%% Define Subfigure Positions
myFig.Position(3) = 1360;
myFig.Position(4) = 960;

schematicPosition = [0.0,0.5,0.5,0.5];
dynamicsPosition = [0.5,0.4,0.5,0.6];
drcPosition = [0,0,0.4,0.45];
aminPosition = [0.4,0,0.3,0.4];
frankPosition = [0.7,0,0.3,0.4];

%% 1A: Plot schematics Key
axKey = axes(myFig);
set(axKey,'Box','off','XColor','none','YColor','none','Position',schematicPosition,'Color','none')%58->5

xlim([-5.9,5.7])
ylim([-5.8,5.3])
pbaspect([1,1,1])
receptorGrey1 = [1,1,1]*0.95;
receptorGrey2 = [1,1,1]*0.7;
nameArray = {'Position','Curvature'};

% Draw Molecular CSC
CSCScale = 1.2;
CSCCenterX = -3/CSCScale;
CSCCenterY = 1.5/CSCScale;
receptorSize = 0.8;

r1 = drawShadedRectangle(axKey,'Position',CSCScale*[CSCCenterX+1,CSCCenterY-receptorSize/2*cos(pi/12)/cos(pi/4),receptorSize,receptorSize],'none',receptorGrey1,receptorGrey2,-pi/6);
r2 = drawShadedRectangle(axKey,'Position',CSCScale*[CSCCenterX+1,CSCCenterY+receptorSize/2*cos(pi/12)/cos(pi/4),receptorSize,receptorSize],'none',receptorGrey1,receptorGrey2,pi/6);
r3 = drawShadedRectangle(axKey,'Position',CSCScale*[CSCCenterX+1+receptorSize*(1/2+1/2*cos(pi/12)/cos(pi/4)),CSCCenterY,receptorSize,receptorSize],'none',receptorGrey1,receptorGrey2);
r4 = drawShadedRectangle(axKey,'Position',CSCScale*[CSCCenterX-1,CSCCenterY-receptorSize/2*cos(pi/12)/cos(pi/4),receptorSize,receptorSize],'none',receptorGrey1,receptorGrey2,pi/6);
r5 = drawShadedRectangle(axKey,'Position',CSCScale*[CSCCenterX-1,CSCCenterY+receptorSize/2*cos(pi/12)/cos(pi/4),receptorSize,receptorSize],'none',receptorGrey1,receptorGrey2,-pi/6);
r6 = drawShadedRectangle(axKey,'Position',CSCScale*[CSCCenterX-1-receptorSize*(1+cos(pi/12)/cos(pi/4))/2,CSCCenterY,receptorSize,receptorSize],'none',receptorGrey1,receptorGrey2);

kinaseGrey1 = [1,1,1]*0.85;
kinaseGrey2 = [1,1,1]*0.4;

displaceX = 0.6;
displaceY = 1.3;
kinaseConnect1 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX+0.3,CSCCenterY+0.6,0.1,1],[0.5*0.8,0.5*0.8]},'slopeDown',kinaseGrey1,kinaseGrey2);
kinaseConnect2 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX-0.3,CSCCenterY-0.6,0.1,1],[0.5*0.8,0.5*0.8]},'slopeDown',kinaseGrey1,kinaseGrey2);
kinase1 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX,CSCCenterY,1,1],[1,1]},'none',kinaseGrey1,kinaseGrey2);
kinase2 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX+displaceX,CSCCenterY+displaceY,1.1,0.6],[1,1]},'none',kinaseGrey1,kinaseGrey2,pi/6);
kinase3 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX-displaceX,CSCCenterY-displaceY,1.1,0.6],[1,1]},'none',kinaseGrey1,kinaseGrey2,pi/6);
CheW1 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX+displaceX,CSCCenterY-displaceY,1.0,0.6],[0.5*0.8,0.5*0.8]},'none',[1,1,1],kinaseGrey2,-pi/6);
CheW2 = drawShadedRectangle(axKey,nameArray,{CSCScale*[CSCCenterX-displaceX,CSCCenterY+displaceY,1.0,0.6],[0.5*0.8,0.5*0.8]},'none',[1,1,1],kinaseGrey2,-pi/6);


% Make legend
kinaseCoordinates = [3.5, CSCCenterY*CSCScale];
receptorCoordinates = [kinaseCoordinates(1), CSCCenterY*CSCScale];

drawShadedRectangle(axKey,nameArray,{[receptorCoordinates(1),receptorCoordinates(2),4,4],[0,0]},'none',receptorGrey1,receptorGrey2);
text(receptorCoordinates(1),receptorCoordinates(2)+1.5,'Receptor','FontSize',20,'HorizontalAlignment','center')

drawShadedRectangle(axKey,nameArray,{[kinaseCoordinates(1),kinaseCoordinates(2),2,2],[1,1]},'none',kinaseGrey1,kinaseGrey2);
text(kinaseCoordinates(1),kinaseCoordinates(2),'Kinase','FontSize',20,'HorizontalAlignment','center')

text(0,CSCCenterY*CSCScale+2.8,'Core Signaling Unit','HorizontalAlignment','center')

equalSignX = 0.6;
equalSignTop = drawShadedRectangle(axKey,'Position',[equalSignX,CSCCenterY*CSCScale+0.2,1,0.1],'none',[0,0,0],[0,0,0]);
equalSignBottom = drawShadedRectangle(axKey,'Position',[equalSignX,CSCCenterY*CSCScale-0.2,1,0.1],'none',[0,0,0],[0,0,0]);

% Draw LU CSCs
offColor = 0.8;
LUCSCHeight = -3.5;

receptorBlue1 = [offColor,offColor,1]*0.9; %modifier for greyscale
receptorBlue2 = [offColor,offColor,1]*0.8;
receptorRed1 = [1,offColor,offColor];
receptorRed2 = [1,offColor,offColor]*0.8;
llR = drawShadedRectangle(axKey,'Position',[-4.7,LUCSCHeight,2,2],'none',receptorBlue1,receptorBlue2);
lmR = drawShadedRectangle(axKey,'Position',[-1.7,LUCSCHeight,2,2],'none',receptorRed1,receptorRed2);
rmR = drawShadedRectangle(axKey,'Position',[1.7,LUCSCHeight,2,2],'none',receptorGrey1,receptorGrey2);
rrR = drawShadedRectangle(axKey,'Position',[4.7,LUCSCHeight,2,2],'none',receptorGrey1,receptorGrey2);

kinaseBlue1 = [0.3,0.3,1]*0.8; %modifier for greyscale
kinaseBlue2 = [0.3,0.3,1]*0.7;
kinaseRed1 = [1,0.3,0.3];
kinaseRed2 = [1,0.3,0.3]*0.7;
llK = drawShadedRectangle(axKey,nameArray,{[-4.7,LUCSCHeight,1,1],1},'none',kinaseGrey1,kinaseGrey2);
lmK = drawShadedRectangle(axKey,nameArray,{[-1.7,LUCSCHeight,1,1],1},'none',kinaseGrey1,kinaseGrey2);
rmK = drawShadedRectangle(axKey,nameArray,{[1.7,LUCSCHeight,1,1],1},'none',kinaseBlue1,kinaseBlue2);
rrK = drawShadedRectangle(axKey,nameArray,{[4.7,LUCSCHeight,1,1],1},'none',kinaseRed1,kinaseRed2);


% Make State Labels
text(-3.2,LUCSCHeight+1.5,'Receptor States','HorizontalAlignment','center')
text(3.2,LUCSCHeight+1.5,'Kinase States','HorizontalAlignment','center')
text(-4.7,LUCSCHeight-1.5,'Inactivating','HorizontalAlignment','center')
text(-1.7,LUCSCHeight-1.5,'Activating','HorizontalAlignment','center')
text(1.7,LUCSCHeight-1.5,'Inactive','HorizontalAlignment','center')
text(4.7,LUCSCHeight-1.5,'Active','HorizontalAlignment','center')


%% 1B: Plot Dynamics
axDynamics = axes(myFig);
set(gca,'Box','off','XColor','none','YColor','none','Color','none','Position',dynamicsPosition)%.5->.58
pbaspect([2,1.53,1])
xlim([0,20])
ylim([-0.3,14])

receptorArray = [0 0 0 0; 1 0 0 1; 0 1 0 0; 1 1 0 1]';
kinaseArray1 = [1 1 1 0; 1 0 0 0; 0 1 0 0; 1 1 1 0]';
kinaseArray2 = [0 1 0 0; 1 0 0 0; 0 1 0 0; 1 1 1 1]';
latSize = length(receptorArray);

receptorRed = [1,0.8,0.8];
receptorRed = receptorRed1;
receptorBlue = [0.8,0.8,1];
receptorBlue = receptorBlue1;
kinaseRed = kinaseRed1;
kinaseBlue = kinaseBlue1;

% Make Lattices
boxLineWidth = 5;
spontColor = [34,139,34]/256*1.8;
inactColor = [0,0,1];
actColor = [1,0,0];

gridScale = 2;

xShift = 0;
yShift = 2;
for i = 1:latSize
    for j = 1:latSize
        drawShadedRectangle(axDynamics,'Position',gridScale*[i+xShift,yShift+latSize-j+1,1,1],'none',...
            receptorRed*receptorArray(i,j)+receptorBlue*(1-receptorArray(i,j)),0);
        drawShadedRectangle(axDynamics,{'Position','Curvature'},{gridScale*[i+xShift,yShift+latSize-j+1,0.5,0.5],1,},'none',...
            kinaseRed*kinaseArray1(i,j)+kinaseBlue*(1-kinaseArray1(i,j)),0);
    end
end
drawShadedRectangle(axDynamics,{'Position','LineWidth','EdgeColor'},{gridScale*[1+xShift,latSize+yShift,1,1],boxLineWidth,spontColor});
drawShadedRectangle(axDynamics,{'Position','LineWidth','EdgeColor'},{gridScale*[3+xShift,latSize+yShift-0.5,1,2],boxLineWidth,inactColor});
drawShadedRectangle(axDynamics,{'Position','LineWidth','EdgeColor'},{gridScale*[3.5+xShift,1+yShift,2,1],boxLineWidth,actColor});


xShift = 5;
for i = 1:latSize
    for j = 1:latSize
        drawShadedRectangle(axDynamics,'Position',gridScale*[i+xShift,yShift+latSize-j+1,1,1],'none',...
            receptorRed*receptorArray(i,j)+receptorBlue*(1-receptorArray(i,j)),0);
        drawShadedRectangle(axDynamics,{'Position','Curvature'},{gridScale*[i+xShift,yShift+latSize-j+1,0.5,0.5],1,},'none',...
            kinaseRed*kinaseArray2(i,j)+kinaseBlue*(1-kinaseArray2(i,j)),0);
    end
end
drawShadedRectangle(axDynamics,{'Position','LineWidth','EdgeColor'},{gridScale*[1+xShift,latSize+yShift,1,1],boxLineWidth,spontColor});
drawShadedRectangle(axDynamics,{'Position','LineWidth','EdgeColor'},{gridScale*[3+xShift,latSize+yShift-0.5,1,2],boxLineWidth,inactColor});
drawShadedRectangle(axDynamics,{'Position','LineWidth','EdgeColor'},{gridScale*[3.5+xShift,1+yShift,2,1],boxLineWidth,actColor});

rightArrow = annotation('arrow');
rightArrow.Parent = axDynamics;
set(rightArrow,'Position', gridScale*[xShift-0.375,yShift+2.5,0.75,0],...
    'LineWidth',4,'HeadWidth',30,'HeadLength',20)

% Make Labels
% Spreading Inactivity
inactSpread = annotation('textbox');
set(inactSpread,'Parent',axDynamics,'String',{'Spreading','Inactivity'},'HorizontalAlignment','center','FontSize',20,...
    'LineWidth',boxLineWidth,'EdgeColor',inactColor,'Color',inactColor,'Position',gridScale*[1,yShift-.75,2,1],'VerticalAlignment','middle')
% Spreading Activity
actSpread = annotation('textbox');
set(actSpread,'Parent',axDynamics,'String',{'Spreading','Activity'},'HorizontalAlignment','center','FontSize',20,...
    'LineWidth',boxLineWidth,'EdgeColor',actColor,'Color',actColor,'Position',gridScale*[xShift-1.5,yShift-.75,2,1],'VerticalAlignment','middle')

% Spontaneous Switching
spontSwitchText = text(gridScale*(2),gridScale*(yShift-1.5),{'Spontaneous','Switching'},'HorizontalAlignment','center',...
    'FontSize',20,'Color',spontColor);
spontSwitchBox = rectangle(axDynamics,'LineWidth',boxLineWidth,'EdgeColor',spontColor,'Position',gridScale*[3.5,yShift-2,1,1]);

% Fast and Slow Processes
text(gridScale*(xShift+2.5),gridScale*(yShift-0.25),...
    'Fast Processes','FontSize',20,'HorizontalAlignment','center')
text(gridScale*(xShift+2.5),gridScale*(yShift-1.5),'Slow Process','FontSize',20,'HorizontalAlignment','center')

text(gridScale*(latSize+1)/2,gridScale*(latSize+yShift+0.8),'{\it t = t_0}','FontSize',20,'HorizontalAlignment','center')
text(gridScale*((latSize+1)/2 + xShift),gridScale*(latSize+yShift+0.8),'{\it t = t_0 + {\Delta}t}','FontSize',20,'HorizontalAlignment','center')

% dividers
hold on
plot([1.5,18.5],gridScale*[yShift+0.375,yShift+0.375],'k')
plot([1.5,18.5],gridScale*[yShift-0.875,yShift-0.875],'k')
plot([1.5,18.5],gridScale*[yShift-2.125,yShift-2.125],'k')

%% Define Useful Functions
f  = @(P,L,m) P.n*(-P.ar*(m-P.mr0) + log((1+L./abs(P.Ki))./(1+L/abs(P.Ka))));
p  = @(P,L,m) 1./(1+exp(f(P,L,m)));
g  = @(P,m) P.as*(m-P.ms0);
c  = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));
ActSS = @(P,L,m) (1 - c(P,L,m) + 2*P.k/(1+exp(g(P,m))) + c(P,L,m)*2*P.k/(1+exp(-g(P,m))) ...
    - sqrt((-1 + c(P,L,m) - 2*P.k/(1+exp(g(P,m))) - c(P,L,m)*2*P.k/(1+exp(-g(P,m)))).^2 - 4*(1-c(P,L,m)).*c(P,L,m)*2*P.k/(1+exp(-g(P,m)))))...
    ./(2*(1-c(P,L,m)));

fAsp  = @(P,L,m) P.n*(-P.ar*(m-P.mr0) + log((1+L./abs(P.KiAsp))./(1+L/abs(P.KaAsp))));
pAsp  = @(P,L,m) 1./(1+exp(fAsp(P,L,m)));
gAsp  = @(P,m) P.as*(m-P.ms0);
cAsp  = @(P,L,m) pAsp(P,L,m)./(1-pAsp(P,L,m)).*exp(gAsp(P,m));
ActSSAsp = @(P,L,m) (1 - cAsp(P,L,m) + 2*P.k/(1+exp(gAsp(P,m))) + cAsp(P,L,m)*2*P.k/(1+exp(-gAsp(P,m))) ...
    - sqrt((-1 + cAsp(P,L,m) - 2*P.k/(1+exp(gAsp(P,m))) - cAsp(P,L,m)*2*P.k/(1+exp(-gAsp(P,m)))).^2 - 4*(1-cAsp(P,L,m)).*cAsp(P,L,m)*2*P.k/(1+exp(-gAsp(P,m)))))...
    ./(2*(1-cAsp(P,L,m)));
bind = @(P,L,m) L./(P.KaAsp + L).*pAsp(P,L,m) + L./(P.KiAsp + L).*(1-pAsp(P,L,m));

h = @(P,m) g(P,m) - P.n*(-P.ar*(m-P.mr0));
Q = @(P,m,c) (c.*exp(-h(P,m))).^(1/P.n);
LTocFunc = @(P,L,m) p(P,L,m)./(1-p(P,L,m)).*exp(g(P,m));
cToLFunc = @(P,m,c) (Q(P,m,c)-1) ./ (1/P.Ka - Q(P,m,c)/P.Ki);

% Load data
load ../p53_Fig1Data.mat
% save('p53_Fig1Data.mat','indexValsAmin','indexValsFrank','mWeightsAmin','mWeightsFrank','xvalsAmin','xvalsFrank','yvalsAmin','yvalsFrank')

% define plot colors and sizes
markerSize = 10;
lineWidth = 2;
blue = [0 0.4470 0.7410];
red = myColorOrder(2,:);
grey = [0.5,0.5,0.5];
white = [1,1,1];

% process amin data
b0x = xvalsAmin(1:indexValsAmin(1))*10^3;
b0y = yvalsAmin(1:indexValsAmin(1));
b3x = xvalsAmin(indexValsAmin(1)+1:indexValsAmin(2))*10^3;
b3y = yvalsAmin(indexValsAmin(1)+1:indexValsAmin(2));
k0x = xvalsAmin(indexValsAmin(2)+1:indexValsAmin(3))*10^3;
k0y = yvalsAmin(indexValsAmin(2)+1:indexValsAmin(3));
k3x = xvalsAmin(indexValsAmin(5)+1:indexValsAmin(6))*10^3;
k3y = yvalsAmin(indexValsAmin(5)+1:indexValsAmin(6));

% process Frank data
wtx{1} = xvalsFrank(1:indexValsFrank(1))*10^3;
wtx{2} = xvalsFrank(indexValsFrank(3)+1:indexValsFrank(4))*10^3;
wty{1} = yvalsFrank(1:indexValsFrank(1))*mWeightsFrank(1);
wty{2} = yvalsFrank(indexValsFrank(3)+1:indexValsFrank(4))*mWeightsFrank(4);
x2x{1} = xvalsFrank(indexValsFrank(1)+1:indexValsFrank(2))*10^3;
x2x{2} = xvalsFrank(indexValsFrank(2)+1:indexValsFrank(3))*10^3;
x2y{1} = yvalsFrank(indexValsFrank(1)+1:indexValsFrank(2))*mWeightsFrank(2);
x2y{2} = yvalsFrank(indexValsFrank(2)+1:indexValsFrank(3))*mWeightsFrank(3);

%% Plot 1D: difference in amplification and adaptation between receptors and kinases
L = 10.^(-4:0.01:4); % define ligand values for hill curve inputs

% generate curves from fit
load ../p54_constrainedAminParams.mat
fitParamsStruct.KaAsp = fitParamsStruct.KaAsp*10^3;
fitParamsStruct.KiAsp = fitParamsStruct.KiAsp*10^3;
bindRes{1} = bind(fitParamsStruct,L,0)*fitParamsStruct.bindingScaleFactor/mWeightsAmin(1);
bindRes{2} = bind(fitParamsStruct,L,3)*fitParamsStruct.bindingScaleFactor/mWeightsAmin(2);
actRes{1} = ActSSAsp(fitParamsStruct,L,0)*fitParamsStruct.actScaleFactor/mWeightsAmin(3);
actRes{2} = ActSSAsp(fitParamsStruct,L,3)*fitParamsStruct.actScaleFactor/mWeightsAmin(6);

axAmin = axes(myFig); % make axes
axAmin.OuterPosition = aminPosition; %define axes position

% plot activity and make it pretty
semilogx(L,actRes{1}/max(actRes{1}),'-','LineWidth',lineWidth,'Color',blue)
hold on
semilogx(L,actRes{2}/max(actRes{2}),'-','LineWidth',lineWidth,'Color',blue)
semilogx(k0x,k0y/max(actRes{1}),'o','MarkerFaceColor',blue,'LineWidth',lineWidth,'MarkerSize',markerSize,'Color',blue)
semilogx(k3x,k3y/max(actRes{2}),'o','MarkerFaceColor',white,'LineWidth',lineWidth,'MarkerSize',markerSize,'Color',blue)

% plot binding and make it pretty
semilogx(L,1-bindRes{1}/max(bindRes{1}),'-','LineWidth',lineWidth,'Color',red)
semilogx(L,1-bindRes{2}/max(bindRes{2}),'-','LineWidth',lineWidth,'Color',red)
semilogx(b0x,1-b0y/max(bindRes{1}),'o','MarkerFaceColor',red,'LineWidth',lineWidth,'MarkerSize',markerSize,'Color',red)
semilogx(b3x,1-b3y/max(bindRes{2}),'o','MarkerFaceColor',white,'LineWidth',lineWidth,'MarkerSize',markerSize,'Color',red)

ylabel('Normalized Values')
xlabel('Ligand [Asp] ({\mu}M)')
ylim([-0.15,1.3])

xlim([0.3*10^-3,3*10^3])
set(gca,'FontSize',18,'XTick',10.^(-3:2:3))

legend('\color[rgb]{0, 0.4470, 0.7410} Kinase','','','',...
    '\color[rgb]{0.9775,0.3737,0.1127} Receptor','location','north','Orientation','horizontal')

% legend(['\color[rgb]{0, 0.4470, 0.7410} Kinase'  newline  '\color[rgb]{0, 0.4470, 0.7410} Activity'],'','','',...
%     ['\color[rgb]{0.9775,0.3737,0.1127} Receptor' newline '\color[rgb]{0.9775,0.3737,0.1127} Vacancy'],'location','north','Orientation','horizontal')
legend('boxoff')

% add legend for different methylation levels
axAminMethLabels = axes('Position',get(gca,'Position'));
semilogx(0,0,'o-','MarkerFaceColor',grey,'MarkerEdgeColor',grey,'LineWidth',2,'Color',grey,'MarkerSize',10)
hold on
semilogx(0,0,'o-','MarkerFaceColor',white,'MarkerEdgeColor',grey,'LineWidth',2,'Color',grey,'MarkerSize',10)
lgH = legend('{\itm} = 0','{\itm} = 1.5','location','southwest');
title(lgH,'Methylation','FontWeight','normal')
legend('boxoff')
set(gca,'FontSize',21,'Visible','off')
axAminMethLabels.Color = 'none';



%% 1E: lesser adaptation in CheW-X2
L = 10.^(-2:0.01:4); % define ligand values for hill curve inputs

normalized = 1;

load ../p54_constrainedFrankParams.mat
fitParamsStruct.Ka = fitParamsStruct.Ka*10^3;
fitParamsStruct.Ki = fitParamsStruct.Ki*10^3;


% make figure axis
axX2 = axes(myFig);
axX2.OuterPosition = frankPosition;

% plot fits
if normalized == 0
    semilogx(L,ActSS(fitParamsStruct,L,1),'LineWidth',2,'Color',blue);
    hold on
    semilogx(L,p(fitParamsStruct,L,1),'LineWidth',lineWidth,'Color',red);
    semilogx(L,p(fitParamsStruct,L,4),'LineWidth',lineWidth,'Color',red)
    semilogx(L,ActSS(fitParamsStruct,L,4),'LineWidth',lineWidth,'Color',blue)
elseif normalized == 1
    semilogx(L,ActSS(fitParamsStruct,L,1)/max(ActSS(fitParamsStruct,L,1)),'LineWidth',2,'Color',blue);
    hold on
    semilogx(L,p(fitParamsStruct,L,1)/max(p(fitParamsStruct,L,1)),'LineWidth',lineWidth,'Color',red);
    semilogx(L,p(fitParamsStruct,L,4)/max(p(fitParamsStruct,L,4)),'LineWidth',lineWidth,'Color',red)
    semilogx(L,ActSS(fitParamsStruct,L,4)/max(ActSS(fitParamsStruct,L,4)),'LineWidth',lineWidth,'Color',blue)
end

% plot data
if normalized == 0
    semilogx(wtx{1},wty{1},'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',blue,'Color',blue)
    semilogx(x2x{1},x2y{1},'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',red,'Color',red)
    semilogx(x2x{2},x2y{2},'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',white,'Color',red)
    semilogx(wtx{2},wty{2},'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',white,'Color',blue)
elseif normalized == 1
    semilogx(wtx{1},wty{1}/max(ActSS(fitParamsStruct,L,1)),'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',blue,'Color',blue)
    semilogx(x2x{1},x2y{1}/max(p(fitParamsStruct,L,1)),'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',red,'Color',red)
    semilogx(x2x{2},x2y{2}/max(p(fitParamsStruct,L,4)),'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',white,'Color',red)
    semilogx(wtx{2},wty{2}/max(ActSS(fitParamsStruct,L,4)),'o','LineWidth',lineWidth,'MarkerSize',markerSize,'MarkerFaceColor',white,'Color',blue)
end

% make pretty
legend('\color[rgb]{0, 0.4470, 0.7410} WT','\color[rgb]{0.9775,0.3737,0.1127} CheW-X2','location','north','Orientation','horizontal')
legend('boxoff')
xlabel('Ligand [MeAsp] ({\mu}M)')
if normalized == 0
    ylabel('Kinase Activity')
elseif normalized == 1
    ylabel('Normalized Activity')
end
set(gca,'xlim',[min(L),max(L)],'ylim',[-0.15,1.3],'FontSize',22,'XTick',10.^(-2:2:4))

% add legend for different methylation levels
axX2meth = axes('Position',get(gca,'Position'));
semilogx(0,0,'o-','MarkerFaceColor',grey,'MarkerEdgeColor',grey,'LineWidth',2,'Color',grey,'MarkerSize',markerSize)
hold on
semilogx(0,0,'o-','MarkerFaceColor',white,'MarkerEdgeColor',grey,'LineWidth',2,'Color',grey,'MarkerSize',markerSize)
if normalized == 0
    lgH = legend('{\itm} = 0.5','{\itm} = 2','location','west');
elseif normalized == 1
    lgH = legend('{\itm} = 0.5','{\itm} = 2','location','southwest');
end
title(lgH,'Methylation','FontWeight','normal')
legend('boxoff')
set(gca,'FontSize',21,'Visible','off')
axX2meth.Color = 'none';

%% 1C: Amplification Schematic
load ../p54_constrainedShimizuParams.mat % load parameters for c<->L conversions

% define c and L arrays for plotting
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


%% The actual figures
load('p10b_stochasticDRC.mat','c','drc','kArray') % load data for plotting

% Define complete set of params for each drc curve
for i = 1:length(kArray)
    param{i} = fitParamsStruct;
    param{i}.k = kArray(i);
end

L = cToLFunc(fitParamsStruct,M,c); %  convert c to L for plotting simulation data

% Setup Axis
ax1DRC = axes(myFig); % make axes
hold on
set(gca,'XScale','log','OuterPosition',drcPosition)

% Plot MF Data
set(gca,'ColorOrderIndex',3)
semilogx(LArrayBase,ActSS(param{4},LArrayBase,M),'--','LineWidth',2)
set(gca,'ColorOrderIndex',1)
semilogx(LArrayBase,ActSS(param{3},LArrayBase,M),'--',LArrayBase,ActSS(param{2},LArrayBase,M),'--','linewidth',2) % plot mf results

% Plot Sim Data
set(gca,'ColorOrderIndex',3)
semilogx(L(4,:)',drc(4,:)','linewidth',2.5)
set(gca,'ColorOrderIndex',1)
semilogx(ax1DRC,L(3:-1:2,:)',drc(3:-1:2,:)','linewidth',2.5) % plot simulation data

% make plot pretty
xlabel('Ligand Concentration [MeAsp] (mM)')
ylabel('Mean Activity ({\itA})')
legend(['MF: {\itk} = ' num2str(kArray(4))],['MF: {\itk} = ' num2str(kArray(3))],['MF: {\itk} = ' num2str(kArray(2))],...
    ['sim: {\itk} = ' num2str(kArray(4))], ['sim: {\itk} = ' num2str(kArray(3))],['sim: {\itk} = ' num2str(kArray(2))],...
    'location','northeast')
legend('boxoff')
xMin = 10^-2;
xMax = 3;
xlim([xMin,xMax])
ax1DRC.Box = 'off';
set(ax1DRC,'XTick',10.^(-3:2))

% make second axis to display control parameter values
ax2DRC = axes(myFig);
ax2DRC.XAxisLocation = 'top';
ax2DRC.YAxisLocation = 'right';

% set control parameter values to display
% cVals = [12,10,3,1,0.3,0.15]; 
cVals = [ceilS(LTocFunc(fitParamsStruct,xMin,M),1), round(LTocFunc(fitParamsStruct,10^-1.5,M),1,'significant'),1, ...
    round(LTocFunc(fitParamsStruct,xMax/4,M),1,'significant'),floorS(LTocFunc(fitParamsStruct,xMax,M),1)]; 
LVals = cToLFunc(fitParamsStruct,M,cVals);
for i = 1:length(cVals)
    cStrings{i} = num2str(cVals(i));
end

% display control parameter values on top axis
set(ax2DRC,'XTick',LVals,'XTickLabel',cStrings,'xscale','log')

% make things pretty
xlim([xMin,xMax])
xlabel('Control Parameter, {\itc}')
ax2DRC.Color = 'none';
ax2DRC.YTickLabel = [];

%% Final Touches
% fontsize(myFig,23,'pixels') %set font sizes
textH = findall(gcf,'Type','text');
set(textH,'FontSize',23)

legendH = findall(gcf,'Type','legend');
set(legendH,'FontSize',20)
titleH = findall(gcf,'Type','Title');
set(titleH,'FontSize',20)
axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2,'FontSize',20,'LabelFontSizeMultiplier',1.2);
% get(lgH,'Position');
% lgH.Position(1) = lgH.Position(1)*1.01;
% lgH.Position(2) = lgH.Position(2)*0.99;

% make sure all double x-axes align
set(ax2DRC,'XMinorTick','off','Position',ax1DRC.Position)

% add figure labels
figLabelSize = 30;
annotation('textbox',[schematicPosition(1) schematicPosition(2)+schematicPosition(4) 0 0], ...
    'String','A','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[dynamicsPosition(1) dynamicsPosition(2)+dynamicsPosition(4) 0 0], ...
    'String','B','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[drcPosition(1) drcPosition(2)+drcPosition(4) 0 0], ...
    'String','C','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[aminPosition(1) aminPosition(2)+aminPosition(4) 0 0], ...
    'String','D','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[frankPosition(1) frankPosition(2)+frankPosition(4) 0 0], ...
    'String','E','EdgeColor','none','FontSize',figLabelSize)
% annotation('textbox',[0.5 0.55 0 0], ...
%     'String','F','EdgeColor','none','FontSize',figLabelSize)
% annotation('textbox',[0.79 0.55 0 0], ...
%     'String','G','EdgeColor','none','FontSize',figLabelSize)

%% 1e Insets
% general Setup
load ../Fig1Lats.mat
gridScale = 1.5;
latSize = length(receptorArray{1});
circleScale = 0.4;
xShift = -0.45;
yShift = -0.45;
insetSizeeAxXY = 12.5;

% Define Inset Postions
ogPosition = get(ax1DRC,'Position');
insetHeight = ogPosition(4)*0.4;
insetWidth = insetHeight * myFig.Position(4)/myFig.Position(3);
actPositionMid = [ogPosition(1)+insetWidth/8,ogPosition(2)+insetHeight/8,insetWidth,insetHeight];
actPositionHigh = [ogPosition(1)+insetWidth/8,ogPosition(2)+insetHeight*10/8,insetWidth,insetHeight];
actPositionLow = [ogPosition(1)+ogPosition(3)-insetWidth*9/8,ogPosition(2)+insetHeight/8,insetWidth,insetHeight];

% Define Arrow Positions
lowArrowPosition = [0.7,0.25,-0.15,-0.25+mean(kinaseArray{3},'all')];
midArrowPosition = [0.3,0.25,0.15,-0.25+0.5];
highArrowPosition = [0.3,0.75,0.12,-0.75+mean(kinaseArray{2},'all')];

%% Add High Act lattice
latticeIndex = 2;

axHighAct = axes(myFig);
set(gca,'Box','off','XColor','none','YColor','none','Color','none')
pbaspect([1,1,1])
set(gca,'Position',actPositionHigh)
xlim([0,insetSizeeAxXY])
ylim([0,insetSizeeAxXY])

for i = 1:latSize
    for j = latSize:-1:1
        drawShadedRectangle(axHighAct,'Position',gridScale*[i+xShift,j+yShift,1,1],'none',...
            receptorRed*receptorArray{latticeIndex}(i,j)+receptorBlue*(1-receptorArray{latticeIndex}(i,j)),0);
        drawShadedRectangle(axHighAct,{'Position','Curvature'},{gridScale*[i+xShift,j+yShift,circleScale,circleScale],1,},'none',...
            kinaseRed*kinaseArray{latticeIndex}(i,j)+kinaseBlue*(1-kinaseArray{latticeIndex}(i,j)),0);
    end
end

%% Add Mid Act lattice
latticeIndex = 1;

axMidAct = axes(myFig);
set(gca,'Box','off','XColor','none','YColor','none','Color','none')
ogPosition = get(ax1DRC,'Position');
set(gca,'Position',actPositionMid)
pbaspect([1,1,1])
xlim([0,insetSizeeAxXY])
ylim([0,insetSizeeAxXY])

for i = 1:latSize
    for j = latSize:-1:1
        drawShadedRectangle(axMidAct,'Position',gridScale*[i+xShift,j+yShift,1,1],'none',...
            receptorRed*receptorArray{latticeIndex}(i,j)+receptorBlue*(1-receptorArray{latticeIndex}(i,j)),0);
        drawShadedRectangle(axMidAct,{'Position','Curvature'},{gridScale*[i+xShift,j+yShift,circleScale,circleScale],1,},'none',...
            kinaseRed*kinaseArray{latticeIndex}(i,j)+kinaseBlue*(1-kinaseArray{latticeIndex}(i,j)),0);
    end
end
%% add low act lattice
latticeIndex = 3;
axLowAct = axes(myFig);
set(gca,'Box','off','XColor','none','YColor','none','Color','none')
ogPosition = get(ax1DRC,'Position');
set(gca,'Position',actPositionLow)
pbaspect([1,1,1])
xlim([0,insetSizeeAxXY])
ylim([0,insetSizeeAxXY])

for i = 1:latSize
    for j = latSize:-1:1
        drawShadedRectangle(axLowAct,'Position',gridScale*[i+xShift,j+yShift,1,1],'none',...
            receptorRed*receptorArray{latticeIndex}(i,j)+receptorBlue*(1-receptorArray{latticeIndex}(i,j)),0);
        drawShadedRectangle(axLowAct,{'Position','Curvature'},{gridScale*[i+xShift,j+yShift,circleScale,circleScale],1,},'none',...
            kinaseRed*kinaseArray{latticeIndex}(i,j)+kinaseBlue*(1-kinaseArray{latticeIndex}(i,j)),0);
    end
end

%% Add Arrows
axDRCArrows = axes('Position',get(ax1DRC,'Position'));
set(axDRCArrows,'Color','none','YTicklabel',[],'XTicklabel',[],'XTick',[],'YTick',[]);

arrowLow = annotation("arrow");
set(arrowLow,'Parent',axDRCArrows,'Position',lowArrowPosition);
arrowMid = annotation("arrow");
set(arrowMid,'Parent',axDRCArrows,'Position',midArrowPosition);
arrowHigh = annotation("arrow");
set(arrowHigh,'Parent',axDRCArrows,'Position',highArrowPosition);


%% save figure
exportgraphics(gcf,'Plots/Fig1Complete.pdf','ContentType','vector')

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