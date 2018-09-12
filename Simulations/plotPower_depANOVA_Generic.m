resDir = '<PATH_TO_SIM_RESULTS>'; % insert your own path

model_str = 'S'; % compare savePrefix in Sim_depANOVA_Generic

height = 1.2; % width/golden ratio
width = 2.4/1.618;

scale = 300; % 3.13 inch

xpos = 50;
ypos = 300;

figure; % Create Figure
axes('FontName','Helvetica') % Set axis font style
box('on'); % Define box around whole figure
set(gcf,'Position',[xpos ypos scale*width scale*height]) % Set figure format

subplot(2,2,1)
hold on
cd(resDir)
load(['2way_dep_a_' model_str '_gauss_ftest.mat'])
plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
load(['2way_dep_a_' model_str '_gauss_exact.mat'])
plot(param,p_val,'g-<','color',[42 41 112]./255,'linewidth',1.2)
load(['2way_dep_a_' model_str '_gauss_raw.mat'])
plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
load(['2way_dep_a_' model_str '_gauss_res.mat'])
plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)

set(gcf,'color','w')

xlabel('\theta_A','FontSize',10,'FontName','Helvetica')
title('Gauss. errors','FontSize',11,'FontName','Helvetica','FontWeight','bold','FontAngle','italic')

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'YTick'       , [0.05 0.2:0.2:1], ...   %'YTickLabel'       , [], ... 
  'LineWidth'   , 1         );
  box off
  xlim([min(param) max(param)])

subplot(2,2,2)
hold on
load(['2way_dep_a_' model_str '_exp_ftest.mat'])
plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
load(['2way_dep_a_' model_str '_exp_exact.mat'])
plot(param,p_val,'g-<','color',[42 41 112]./255,'linewidth',1.2)
load(['2way_dep_a_' model_str '_exp_raw.mat'])
plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
load(['2way_dep_a_' model_str '_exp_res.mat'])
plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)

xlabel('\theta_A','FontSize',10,'FontName','Helvetica')
title('Exp. errors','FontSize',11,'FontName','Helvetica','FontWeight' , 'bold','FontAngle','italic')


hleg1 = legend('F-test','exact','raw','res');

% Set Legend Properties
set(hleg1,'Location','SouthEast')
set(hleg1,'FontName','Helvetica')
set(hleg1,'FontSize',9)

set(hleg1,'box','off')

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'YTick'       , [0.05 0.2:0.2:1], ...   %'YTickLabel'       , [], ... 
  'YTickLabel'       , [], ... 
'LineWidth'   , 1         );
box off
xlim([min(param) max(param)])

%% Second row: interaction factor
cd(resDir)

subplot(2,2,3)
hold on
load(['2way_dep_iaxb_' model_str '_gauss_ftest.mat'])
plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
load(['2way_dep_iaxb_' model_str '_gauss_raw.mat'])
plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
load(['2way_dep_iaxb_' model_str '_gauss_res.mat'])
plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)

xlabel('\theta_{AB}','FontSize',10,'FontName','Helvetica')


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'YTick'       , [0.05 0.2:0.2:1], ...   %'YTickLabel'       , [], ... 
'LineWidth'   , 1         );
box off
xlim([min(param) max(param(1:end))])

subplot(2,2,4)
hold on
load(['2way_dep_iaxb_' model_str '_exp_ftest.mat'])
plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
load(['2way_dep_iaxb_' model_str '_exp_raw.mat'])
plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
load(['2way_dep_iaxb_' model_str '_exp_res.mat'])
plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)

xlabel('\theta_{AB}','FontSize',10,'FontName','Helvetica')

hleg1 = legend('F-test','raw','res');

set(hleg1,'Location','SouthEast')
set(hleg1,'FontName','Helvetica')
set(hleg1,'FontSize',9)

set(hleg1,'box','off')

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'YTick'       , [0.05 0.2:0.2:1], ... 
  'YTickLabel'       , [], ... 
'LineWidth'   , 1         );
box off
xlim([min(param) max(param(1:end))])

fig = gcf;
style = hgexport('factorystyle');
style.Bounds = 'loose';
hgexport(fig,'Figure.eps',style,'applystyle', true);
drawnow;

print('-depsc2 ', '-tiff', ['depAnova_' model_str '.eps'])


