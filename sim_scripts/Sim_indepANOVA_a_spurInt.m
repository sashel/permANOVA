%% Test for spurious interactions in an independent, balanced 2-way ANOVA without interaction with increasing main effect A
% Uses FtSimLink_indepANOVA, anova2_cell_mod
addpath('<PATH_TO_ADAPTED_FIELDTRIPTOOLBOX>')
ft_defaults

stemFolder = '<INSERT_YOUR_OWN_WORKING_DIR';
addpath([stemFolder 'stat_util/']) 
resDir = [stemFolder 'SimData/ResultsIndepANOVA/']; 


savePrefix = '2way_indep_spurInt';

if ~exist(resDir,'dir')
    mkdir(resDir);
end;

Rep = 1000; % number of Monte Carlo simulations
 
a = 3; % number of levels in factor A
b = 4; % number of levels in factor B
n = 5; % number of repetitions in each factor level combination

B = [50 -50 20 -20];
B = reshape(repmat(B,a,n),a,b,n);

design = [repmat([1,2,3],1,b*n);repmat([1 1 1 2 2 2 3 3 3 4 4 4],1,n)]; % design matrix

load([stemFolder filesep 'dummy'],'data');
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg_neighb, data);

Method_List = {'ftest','raw','res'};
Error_List = {'exp','gauss'};
for ee = 1:length(Error_List)
    figure
    hold on
    for mm = 1:length(Method_List)
        
        method = Method_List{mm};
        err_dist = Error_List{ee};
        
        % same as for Sim_indepANOVA_a
        if strcmp(err_dist,'exp')
            sc = [10000, 40, 20, 12, 10, 8, 6, 4, 3];
        elseif strcmp(err_dist,'gauss')
            sc = [100000, 200, 100, 70, 50, 40, 35, 30, 25]*3;
        end
        
        saveStr = [savePrefix,'_',err_dist,'_',method];

        for r = 1:length(sc)
            % gradually increase effecr size of factor A
            if r ==1
                A = [0; 0; 0];
            else
                A = [50; 0; -50]/sc(r);
            end
            m_A = mean(A);
            par = sqrt(sum(((A - m_A).^2)./a));
            A = reshape(repmat(A,b,n),a,b,n);
           
                for j = 1:Rep
                    
                    % put the model together & add errors
                    if strcmp(err_dist,'exp')
                        y = A + B + (exprnd(1,a,b,n)).^3;
                    elseif strcmp(err_dist,'gauss')
                        y = A + B + randn(a,b,n);
                    end
                    
                    switch method
                        case 'ftest' % standard F-test
                            
                            % prepare data input for anova2_cell_mod
                            c{1,1} = squeeze(y(1,1,:))';
                            c{1,2} = squeeze(y(1,2,:))';
                            c{1,3} = squeeze(y(1,3,:))';
                            c{1,4} = squeeze(y(1,4,:))';
                            c{2,1} = squeeze(y(2,1,:))';
                            c{2,2} = squeeze(y(2,2,:))';
                            c{2,3} = squeeze(y(2,3,:))';
                            c{2,4} = squeeze(y(2,4,:))';
                            c{3,1} = squeeze(y(3,1,:))';
                            c{3,2} = squeeze(y(3,2,:))';
                            c{3,3} = squeeze(y(3,3,:))';
                            c{3,4} = squeeze(y(3,4,:))';

                            
                            [FA, FB, FI, dfa, dfb, dfi] =  anova2_cell_mod(c); % adapted from the resampling statistic toolbox
                            res_I(j) = 1 - fcdf(FI, dfi(1), dfi(2));
                            
                        case 'raw' % permutation of raw data
                            fac = 'iaxb';
                            c = reshape(y,1,a*b*n);
                            exact = 'no';
                            statfun = 'indepAnova2way';
                            % permutation ANOVA is called here
                            stat = FtSimLink_indepANOVA(data,neighbours,c,design,statfun,fac,exact);
                            res_I(j) = stat.prob;
                            
                        case 'res' % permutation of residuals                    
                            fac = 'iaxb';
                            c = reshape(y,1,a*b*n);
                            ncond_a = a;
                            ncond_b = b;
                            
                            for nfac_a = 1:ncond_a
                                for nfac_b = 1:ncond_b
                                    idx_ab = design(1,:) == nfac_a & design(2,:) == nfac_b;
                                    anovaIn{nfac_a,nfac_b} = c(idx_ab);
                                end
                            end
                            
                            for ii = 1:size(anovaIn,1)
                                tmp = zeros(size(anovaIn{1,1}));
                                for jj = 1:size(anovaIn,2)
                                    tmp(:,:,jj) = anovaIn{ii,jj};
                                end
                                amean(ii) = squeeze(mean(mean(tmp)));
                            end
                            
                            for jj = 1:size(anovaIn,2)
                                tmp = zeros(size(anovaIn{1,1}));
                                for ii = 1:size(anovaIn,1)
                                    tmp(:,:,ii) = [(anovaIn{ii,jj})];
                                end
                                bmean(jj) = squeeze(mean(mean(tmp)));
                            end
                            
                            for ii = 1:size(anovaIn,1)
                                tmp = zeros(size(anovaIn{1,1}));
                                for jj = 1:size(anovaIn,2)
                                    tmp(:,:,ii,jj) = anovaIn{ii,jj};
                                end
                                
                            end
                            totmean = squeeze(mean(mean(mean(tmp))));
                            
                            for ii = 1:size(anovaIn,1)
                                for jj = 1:size(anovaIn,2)
                                    idx_ab = design(1,:) == ii & design(2,:) == jj;
                                    c_new(idx_ab) = c(idx_ab) - bmean(jj) - amean(ii) + totmean;
                                    
                                end
                            end
                            
                            exact = 'no';
                            statfun = 'indepAnova2way';
                            stat = FtSimLink_indepANOVA(data,neighbours,c_new,design,statfun,fac,exact);
                            res_I(j) = stat.prob;
                            
                    end
                    
                end
                p_val(r) = length(find(res_I <= 0.05))/Rep;
            paramA(r) = par;
        end
        
        switch method
            case 'ftest'
                plot(paramA,p_val,':d','color',[0.4 0.4 0.4],'linewidth',2)
            case 'raw'
                plot(paramA,p_val,'g:d','linewidth',2)
            case 'res'
                plot(paramA,p_val,'m:s','linewidth',2)
        end
        
        set(gcf,'color','w')
        xlabel('\theta_A','FontSize',10,'FontName','Helvetica')
        ylabel('Power','FontSize',10,'FontName','Helvetica')
        if strcmp(err_dist,'gauss')
            title('Gauss. errors','FontSize',11,'FontName','Helvetica','FontWeight','bold','FontAngle','italic')
        elseif strcmp(err_dist,'exp')
            title('Exp. errors','FontSize',11,'FontName','Helvetica','FontWeight','bold','FontAngle','italic')
        end
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
        xlim([min(paramA) max(paramA)])

        cd(resDir)
        save(saveStr,'paramA','p_val*','A','B','sc','method')
        clear paramA p_val* c
        cd(stemFolder)
    end
end


