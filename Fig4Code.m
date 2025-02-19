% This script generate figure 4.
close all
clear

nHistBins = 20; % set consistent bin size for histograms
%% Make Initial Figure

f = figure('Position',[1000,500,1800,800]);

% define big axis to help with sizing and positioning of smaller subplots
for i = 1:3
    ax(i) = axes(f);
    ax(i).OuterPosition = [(i-1)/3, 0, 1/3, 1];
    set(gca,'Box','off','XColor','none','YColor','none','Color','none')
    % xlabel('Time (s)','Fontsize',20,'Color',[.15 .15 .15])
    % ylabel('Activity','Fontsize',20,'Color',[.15 .15 .15])
    boundaries(i,:) = ax(i).Position;
end
%% 1A: switching w/ k
panelIndex = 1;
numPlots = 5;
load('Fig4Data.mat') % load data

% for-loop to make individual subplots
for i = 1:numPlots
    ax1(i) = axes(f); % make axis
    plot(Fig4Data.kSeries.ts{i},Fig4Data.kSeries.As{i},'linewidth',2) % plot time series data

    % make plot pretty
    set(gca,'fontsize',14)
    ylim([-0.1,1.1])
    xlim([0,2000])

    ax2(i) = axes(f); % make second axis
    histogram(Fig4Data.kSeries.As{i},nHistBins,'Orientation','horizontal') % plot histogram

    % make plot pretty
    % set(ax2(i),'Box','off','XColor','none','Ycolor','none','FontSize',14) 
    set(ax2(i),'XTickLabel',[],'YTickLabel',[],'XTick',[])
    title(['{\itk} = ' num2str(Fig4Data.kSeries.param{i}.k)],'Rotation',270,'Units','normalized','Position',[1.05,0.5],'FontWeight','normal')
    if i ~= numPlots
        ax1(i).XTickLabel = [];
    end
end

% set positions for all subplots
for i = 1:numPlots
    ax1(i).Position(1) = boundaries(panelIndex,1);
    ax1(i).Position(2) = boundaries(panelIndex,2) + boundaries(panelIndex,4)*(numPlots-i)/5;
    ax1(i).Position(3) = boundaries(panelIndex,3)*0.8;
    ax1(i).Position(4) = boundaries(panelIndex,4)/6;
    ax2(i).Position(1) = ax1(i).Position(1) + ax1(i).Position(3) + 0.005;
    ax2(i).Position(3) = boundaries(panelIndex,3)*0.2;
    ax2(i).Position([2,4]) = ax1(i).Position([2,4]);
end

% make histogram and time series axis align
linkaxes([ax1 ax2], 'y')

%% 4b: switching with L
clear ax1 ax2
panelIndex = 2;
numPlots = 4;
% load('p10b_switchingWithLatSize_small.mat') % load data

% for-loop to make individual subplots
for i = 1:numPlots
    ax1(i) = axes(f); % make axis
    plot(Fig4Data.LatSizeSeries.ts{i},Fig4Data.LatSizeSeries.As{i},'linewidth',2) % plot time series data

    % make plot pretty
    set(gca,'fontsize',14)
    ylim([-0.1,1.1])
    xlim([0,50])

    ax2(i) = axes(f);% make second axis
    histogram(Fig4Data.LatSizeSeries.As{i},nHistBins,'Orientation','horizontal') % plot histogram

    % make plot pretty
    % set(ax2(i),'Box','off','XColor','none','Ycolor','none','FontSize',14)
    set(ax2(i),'XTickLabel',[],'YTickLabel',[],'XTick',[])
    title(['{\itL} = ' num2str(length(Fig4Data.LatSizeSeries.param{i}.A0))],'Rotation',270,'Units','normalized','Position',[1.05,0.5],'FontWeight','normal')
    if i ~= numPlots
        ax1(i).XTickLabel = [];
    end
end

% set positions for all subplots
for i = 1:numPlots
    ax1(i).Position(1) = boundaries(panelIndex,1);
    ax1(i).Position(2) = boundaries(panelIndex,2) + boundaries(panelIndex,4)*(numPlots-i)/numPlots;
    ax1(i).Position(3) = boundaries(panelIndex,3)*0.8;
    ax1(i).Position(4) = boundaries(panelIndex,4)/(numPlots+1);
    ax2(i).Position(1) = ax1(i).Position(1) + ax1(i).Position(3) + 0.005;
    ax2(i).Position(3) = boundaries(panelIndex,3)*0.2;
    ax2(i).Position([2,4]) = ax1(i).Position([2,4]);
end

% make histogram and time series y-axes align
linkaxes([ax1 ax2], 'y')

%% Make plot with external noise
clear ax1 ax2
numPlots = 5;
% load('p10b_switchingWithExtFluct_small.mat') % load data

% for-loop to make individual subplots
for i = 1:numPlots
    ax1(i) = axes(f); % make axis
    plot(Fig4Data.ExtFluct.ts{i},Fig4Data.ExtFluct.As{i},'linewidth',2) % plot time series data

    % make plots pretty
    set(gca,'fontsize',14)
    ylim([-0.1,1.1])
    xlim([0,2000])

    ax2(i) = axes(f); % make second axis
    histogram(Fig4Data.ExtFluct.As{i},nHistBins,'Orientation','horizontal') % plot histogram

    % make plots pretty
    % set(ax2(i),'Box','off','XColor','none','Ycolor','none','FontSize',14)
    set(ax2(i),'XTickLabel',[],'YTickLabel',[],'XTick',[])
    if i == 1
        title('\sigma^2_{{\itx}_{OU}} = 0','Rotation',270,'Units','normalized','Position',[1.05,0.5],'FontWeight','normal')
    elseif i == 5
        title('\sigma^2_{{\itx}_{OU}} = 1','Rotation',270,'Units','normalized','Position',[1.05,0.5],'FontWeight','normal')
    else
        title(['\sigma^2_{{\itx}_{OU}} = 10^{' num2str(-5 + i) '}'],'Rotation',270,'Units','normalized','Position',[1.05,0.5],'FontWeight','normal')
    end
    if i ~= numPlots
        ax1(i).XTickLabel = [];
    end
end

% set positions for all subplots
for i = 1:numPlots
    ax1(i).Position(1) = boundaries(3,1);
    ax1(i).Position(2) = boundaries(3,2) + boundaries(3,4)*(numPlots-i)/numPlots;
    ax1(i).Position(3) = boundaries(3,3)*0.8;
    ax1(i).Position(4) = boundaries(3,4)/(numPlots+1);
    ax2(i).Position(1) = ax1(i).Position(1) + ax1(i).Position(3) + 0.005;
    ax2(i).Position(3) = boundaries(3,3)*0.2;
    ax2(i).Position([2,4]) = ax1(i).Position([2,4]);
end

% make histogram and time series y-axes align
linkaxes([ax1 ax2], 'y')

%% Make Titles and Figure Labels
fontsize(f,20,'pixels') % set fontSizes

for i = 1:3
    xlabel(ax(i), 'Time (s)','Fontsize',23,'Color',[.15 .15 .15])
    ylabel(ax(i), 'Mean Activity ({\itA})','Fontsize',23,'Color',[.15 .15 .15])
end

axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2);

% make titles
ax(1).Title.String = {'Decreasing {\itk}', '({\itL} = 16)'};
ax(1).Title.FontSize = 23;
ax(2).Title.String = {'Decreasing Lattice Size', '({\itk} = 0.01)'};
ax(2).Title.FontSize = 23;
ax(3).Title.String = {'Increasing Extrinsic Noise', '({\itL} = 16, {\itk} = 0.01)'};
ax(3).Title.FontSize = 23;

% make figure labels
figLabelSize = 30;
annotation('textbox',[0 1 0 0], ...
    'String','A','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[0.33 1 0 0], ...
    'String','B','EdgeColor','none','FontSize',figLabelSize)
annotation('textbox',[.66 1 0 0], ...
    'String','C','EdgeColor','none','FontSize',figLabelSize)

% exportgraphics(f,'Plots/Fig4Complete.pdf','ContentType','vector')