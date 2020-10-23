% stemFolder = '<INSERT_YOUR_OWN_WORKING_DIR';
stemFolder = '/data/pt_np-helbling/permANOVA/';
addpath([stemFolder 'stat_util/'])
resDir = [stemFolder 'SimData/ResultsMixedANOVA/'];

model = 'S';

totunres_flag = 1; % plot power
effect = {'confMain'}; % 'mainInt', 'spurInt', 'confMain'

for eff = 1:length(effect)
    
    if strcmp(effect{eff},'mainInt')
        height = 1.8;
        nmb_rows = 3;
    elseif strcmp(effect{eff},'confMain')
        height = 1.2;
        nmb_rows = 2;
    elseif strcmp(effect{eff},'spurInt')
        height = 0.6;
        nmb_rows = 1;
    end
    
    width = 2.4/1.618;
    scale = 300;
    
    xpos = 50;
    ypos = 300;
    
    figure; % Create Figure
    axes('FontName','Helvetica') % Set axis font style
    box('on'); % Define box around whole figure
    set(gcf,'Position',[xpos ypos scale*width scale*height]) % Set figure format
    
    if strcmp(effect{eff},'mainInt')
        model_str = ['a_' model];
    elseif strcmp(effect{eff},'confMain')
        model_str = ['a_' model '_AB'];
    elseif strcmp(effect{eff},'spurInt')
        model_str = ['iaxb_' model];
    end
    
    % open file to write alpha errors to
    fileID = fopen(['mixAnov_' model '_' effect{eff} '.txt'],'w');
    fprintf(fileID,'Model: %s\n',model_str);
    
    %% First row: between factor
    
    subplot(nmb_rows,2,1)
    hold on
    cd(resDir)
    if strcmp(effect{eff},'mainInt')
        fprintf(fileID,'\nMain effects A\n');
    elseif strcmp(effect{eff},'confMain')
        fprintf(fileID,'\nMain effects A - confounding interaction\n');
    elseif strcmp(effect{eff},'spurInt')
        fprintf(fileID,'\nSpurious interaction\n');
    end
    load(['2way_mix_' model_str '_gauss_ftest.mat'])
    plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
    fprintf(fileID,'gauss ftest: %1.3f\n',p_val(1));
    load(['2way_mix_' model_str '_gauss_exact.mat'])
    plot(param,p_val,'g-<','color',[42 41 112]./255,'linewidth',1.2)
    fprintf(fileID,'gauss exact: %1.3f\n',p_val(1));
    load(['2way_mix_' model_str '_gauss_raw.mat'])
    plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
    fprintf(fileID,'gauss raw: %1.3f\n',p_val(1));
    load(['2way_mix_' model_str '_gauss_res.mat'])
    plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)
    fprintf(fileID,'gauss res: %1.3f\n',p_val(1));
    if totunres_flag
        load(['2way_mix_' model_str '_gauss_totunres.mat'])
        plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
        fprintf(fileID,'gauss totunres: %1.3f\n',p_val(1));
    end
    
    set(gcf,'color','w')
    xlabel('\theta_A','FontSize',10,'FontName','Helvetica')
    title('Gauss. errors','FontSize',11,'FontName','Helvetica','FontWeight','bold','FontAngle','italic')
    ylabel('Power','FontSize',10,'FontName','Helvetica')
    
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'      , ...
        'YGrid'       , 'off'      , ...
        'YTick'       , [0.05 0.2:0.2:1], ...
        'LineWidth'   , 1         );
    box off
    xlim([min(param) max(param)])
    
    subplot(nmb_rows,2,2)
    hold on
    load(['2way_mix_' model_str '_exp_ftest.mat'])
    plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
    fprintf(fileID,'exp ftest: %1.3f\n',p_val(1));
    load(['2way_mix_' model_str '_exp_exact.mat'])
    plot(param,p_val,'g-<','color',[42 41 112]./255,'linewidth',1.2)
    fprintf(fileID,'exp exact: %1.3f\n',p_val(1));
    load(['2way_mix_' model_str '_exp_raw.mat'])
    plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
    fprintf(fileID,'exp raw: %1.3f\n',p_val(1));
    load(['2way_mix_' model_str '_exp_res.mat'])
    plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)
    fprintf(fileID,'exp raw: %1.3f\n',p_val(1));
    if totunres_flag
        load(['2way_mix_' model_str '_exp_totunres.mat'])
        plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
        fprintf(fileID,'exp raw: %1.3f\n',p_val(1));
    end
    
    
    set(gcf,'color','w')
    if strcmp(effect{eff},'mainInt')||strcmp(effect{eff},'confMain')
        xlabel('\theta_A','FontSize',10,'FontName','Helvetica')
    elseif strcmp(effect,'spurInt')
        xlabel('\theta_{AB}','FontSize',10,'FontName','Helvetica')
    end
    title('Exp. errors','FontSize',11,'FontName','Helvetica','FontWeight' , 'bold','FontAngle','italic')
    
    % Create Legend
    hleg1 = legend('F-test','exact','raw','res','totunres');
    
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
    xlim([min(param) max(param)])
    
    %% Second row: within factor
    if strcmp(effect{eff},'mainInt')||strcmp(effect{eff},'confMain')
        if strcmp(effect{eff},'mainInt')
            fprintf(fileID,'\nMain effects B\n');
        elseif strcmp(effect{eff},'confMain')
            fprintf(fileID,'\nMain effects B - confounding interaction\n');
        end
        
        if strcmp(effect{eff},'mainInt')
            model_str = ['b_' model];
        elseif strcmp(effect{eff},'confMain')
            model_str = ['b_' model '_AB'];
        end
        subplot(nmb_rows,2,3)
        fprintf(fileID,'\nMain effects B\n');
        hold on
        % cd([stemFolder 'ResultsMixedANOVA'])
        load(['2way_mix_' model_str '_gauss_ftest.mat'])
        plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
        fprintf(fileID,'gauss ftest: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_gauss_exact.mat'])
        plot(param,p_val,'g-<','color',[42 41 112]./255,'linewidth',1.2)
        fprintf(fileID,'gauss exact: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_gauss_raw.mat'])
        plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
        fprintf(fileID,'gauss raw: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_gauss_res.mat'])
        plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)
        fprintf(fileID,'gauss res: %1.3f\n',p_val(1));
        if totunres_flag
            load(['2way_mix_' model_str '_gauss_totunres.mat'])
            plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
            fprintf(fileID,'gauss totunres: %1.3f\n',p_val(1));
        end
        
        set(gcf,'color','w')
        xlabel('\theta_B','FontSize',10,'FontName','Helvetica')
        % title('Gauss. errors','FontSize',11,'FontName','Helvetica','FontWeight','bold','FontAngle','italic')
        ylabel('Power','FontSize',10,'FontName','Helvetica')
        
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.02 .02] , ...
            'XMinorTick'  , 'off'      , ...
            'YMinorTick'  , 'off'      , ...
            'YGrid'       , 'off'      , ...
            'YTick'       , [0.05 0.2:0.2:1], ...
            'LineWidth'   , 1         );
        box off
        xlim([min(param) max(param)])
        
        subplot(nmb_rows,2,4)
        hold on
        load(['2way_mix_' model_str '_exp_ftest.mat'])
        plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2)
        fprintf(fileID,'exp ftest: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_exp_exact.mat'])
        plot(param,p_val,'g-<','color',[42 41 112]./255,'linewidth',1.2)
        fprintf(fileID,'exp exact: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_exp_raw.mat'])
        plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2)
        fprintf(fileID,'exp raw: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_exp_res.mat'])
        plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2)
        fprintf(fileID,'exp res: %1.3f\n',p_val(1));
        if totunres_flag
            load(['2way_mix_' model_str '_exp_totunres.mat'])
            plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
            fprintf(fileID,'exp totunres: %1.3f\n',p_val(1));
        end
        
        xlabel('\theta_B','FontSize',10,'FontName','Helvetica')
        % title('Exp. errors','FontSize',11,'FontName','Helvetica','FontWeight' , 'bold','FontAngle','italic')
        
        % Create Legend
        hleg1 = legend('F-test','exact','raw','res','totunres');
        
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
        xlim([min(param) max(param)])
    end
    %% Third row: interaction factor
    if strcmp(effect{eff},'mainInt')
        fprintf(fileID,'\nInteraction\n');
        model_str = ['iaxb_' model '_AB'];
        subplot(nmb_rows,2,5)
        hold on
        load(['2way_mix_' model_str '_gauss_ftest.mat'])
        plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2);
        fprintf(fileID,'gauss ftest: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_gauss_raw.mat'])
        plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2);
        fprintf(fileID,'gauss raw: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_gauss_res.mat'])
        plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2);
        fprintf(fileID,'gauss res: %1.3f\n',p_val(1));
        if totunres_flag
            load(['2way_mix_' model_str '_gauss_totunres.mat'])
            plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
            fprintf(fileID,'gauss totunres: %1.3f\n',p_val(1));
        end
        
        
        xlabel('\theta_{AB}','FontSize',10,'FontName','Helvetica')
        ylabel('Power','FontSize',10,'FontName','Helvetica')
        
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
        
        subplot(nmb_rows,2,6)
        hold on
        load(['2way_mix_' model_str '_exp_ftest.mat'])
        plot(param,p_val,'--d','color',[0.4 0.4 0.4],'linewidth',1.2);
        fprintf(fileID,'exp ftest: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_exp_raw.mat'])
        plot(param,p_val,'-s','color',[178 220 38]./255, 'linewidth',1.2);
        fprintf(fileID,'exp raw: %1.3f\n',p_val(1));
        load(['2way_mix_' model_str '_exp_res.mat'])
        plot(param,p_val,'o-','color',[226 40 12]./255, 'linewidth',1.2);
        fprintf(fileID,'exp res: %1.3f\n',p_val(1));
        if totunres_flag
            load(['2way_mix_' model_str '_exp_totunres.mat'])
            plot(param,p_val,'o-','color',[153 153 255]./255, 'linewidth',1.2)
            fprintf(fileID,'exp totunres: %1.3f\n',p_val(1));
        end
        
        xlabel('\theta_{AB}','FontSize',10,'FontName','Helvetica')
        
        % Create Legend
        hleg1 = legend('F-test','raw','res','totunres');
        
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
    end
%%
fig = gcf;
style = hgexport('factorystyle');
style.Bounds = 'loose';
hgexport(fig,'Figure.eps',style,'applystyle', true);
drawnow;

% print('-depsc2', '-tiff', ['mixedAnova_' model_str '.eps'])
print( ['mixAnov' '_' model '_' effect{eff}], '-depsc', '-tiff')

fclose(fileID);
end