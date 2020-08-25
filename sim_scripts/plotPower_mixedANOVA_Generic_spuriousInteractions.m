stemFolder = '/data/pt_user-helbling_ticket017439/helbling/permANOVA/permANOVA_Github/';
resDir = [stemFolder 'SimData/ResultsMixedANOVA'];

% compare savePrefix in Sim_mixedANOVA_Generic
model_str = 'S';

%% ADD TOTUNRES
totunres_flag = 1;
iaxb_flag = 0;
height = 0.4;
nmb_rows = 1;

width = 2.4/1.618;
scale = 300; % 3.13 inch

xpos = 50;
ypos = 300;

figure; % Create Figure
axes('FontName','Helvetica') % Set axis font style
box('on'); % Define box around whole figure
set(gcf,'Position',[xpos ypos scale*width scale*height]) % Set figure format


%% Spurious interactions
subplot(1,2,1)
hold on
load(['2way_mix_iaxb_' model_str '_gauss_ftest.mat'])
plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2);
load(['2way_mix_iaxb_' model_str '_gauss_raw.mat'])
plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2);
load(['2way_mix_iaxb_' model_str '_gauss_res.mat'])
plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2);
if totunres_flag
    load(['2way_mix_iaxb_' model_str '_gauss_totunres.mat'])
    plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
end


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
title('Gauss. errors','FontSize',11,'FontName','Helvetica','FontWeight','bold','FontAngle','italic')

subplot(1,2,2)
hold on
load(['2way_mix_iaxb_' model_str '_exp_ftest.mat'])
plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2);
load(['2way_mix_iaxb_' model_str '_exp_raw.mat'])
plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2);
load(['2way_mix_iaxb_' model_str '_exp_res.mat'])
plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2);
if totunres_flag
    load(['2way_mix_iaxb_' model_str '_exp_totunres.mat'])
    plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
end

xlabel('\theta_{AB}','FontSize',10,'FontName','Helvetica')

% Create Legend
hleg1 = legend('F-test','raw','res');

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
  'YTick'       , [0.05 0.2:0.2:1], ... 
  'YTickLabel'       , [], ... 
'LineWidth'   , 1         );
box off
xlim([min(param) max(param(1:end))])
title('Exp. errors','FontSize',11,'FontName','Helvetica','FontWeight' , 'bold','FontAngle','italic')

fig = gcf;
style = hgexport('factorystyle');
style.Bounds = 'loose';
hgexport(fig,'Figure.eps',style,'applystyle', true);
drawnow;

print('-depsc2', '-tiff', ['mixedAnova_' model_str '_spurious_interactions.eps'])



