clear

load('FigSIResampling.mat')
[As,ts] = resample(sim.results.As,sim.results.ts/100,100);

f = figure(1);
clf
f.Position = [500,500,1200,480];
plot(sim.results.ts/100,sim.results.As,'linewidth',2)
hold on
plot(ts,As,'linewidth',2)
xlabel('Time (s)')
ylabel('Activity')
legend('Original Time Points','Resampled Time Points')
ylim([0.2,0.8])


fontsize(f,23,'pixels') % set font sizes

axesH = findall(gcf,'Type','axes');
set(axesH,'LineWidth',2);

% exportgraphics(gcf,'Plots/p10bFigSI1Complete.pdf','ContentType','vector')
