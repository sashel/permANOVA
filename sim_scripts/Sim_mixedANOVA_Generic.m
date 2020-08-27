% Monte Carlo simulations mixed 2-way ANOVA
% Params: fac, Rep

% stemFolder = '<INSERT_YOUR_OWN_WORKING_DIR';
stemFolder = '/data/pt_np-helbling/permANOVA/';
addpath([stemFolder 'stat_util/']) 
resDir = [stemFolder 'SimData/ResultsMixedANOVA/']; 


if ~exist(resDir,'dir')
    mkdir(resDir);
end

fac = 'b'; % factor or interaction of interest: 'a', 'b' and 'iaxb' for main factors A and B, and the interaction.

Int_flag = false; % Interaction between the two factors?

% subjInt_flag = false; % Interaction between within factors and subjects?
% subjInt_scale = 100;

Rep = 1000;

Method_List = {'ftest','exact','raw','res','totunres'}; % Note: not all tests are meaningful for each simulation. E.g. it doesn't make sense to calculate an exact test for the interaction. 
Error_List = {'exp','gauss'};

%% construct save string prefix
savePrefix = sprintf('2way_mix_%s_S',fac);
% if subjInt_flag
%     savePrefix = [savePrefix sprintf('_subjInt_%d',subjInt_scale)];
% end

if Int_flag
    savePrefix = [savePrefix '_AB'];
end

a = 3; % number of levels group factor A
b = 4; % number of levels within factor B
n = a*b; % number of UOs (e.g. subjects)

A = [50; 0; -50]; % see Anderson et al., 2001
B = [50 -50 20 -20]; % 

B = repmat(B,1,n);

A = [repmat(50,1,b*n/a) zeros(1,b*n/a) repmat(-50,1,b*n/a)];  % see Anderson et al., 2001

if Int_flag
    Int = A.*B;
end

load([stemFolder filesep 'dummy'],'data');
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg_neighb, data);

for ee = 1:length(Error_List)
    figure
    hold on
    for mm = 1:length(Method_List)
        
        method = Method_List{mm};
        err_dist = Error_List{ee};
        
        if strcmp(err_dist,'exp')
            if strcmp(fac, 'a')
                sc = [100000, 230,140, 90, 65, 48, 34, 26, 20]/10;
            elseif strcmp(fac, 'b')
                sc = [10000, 280, 160, 100, 60, 40, 30, 20, 15]/10;
            elseif strcmp(fac,'iaxb')
                sc = [100000, 400, 200, 130, 80, 55, 40, 30, 24]*3;
            end
            
        elseif strcmp(err_dist,'gauss')
            if strcmp(fac, 'a')
                sc = [100000, 140, 80,56,44 ,36,  30, 25, 20]*3;
            elseif strcmp(fac, 'b')
                sc = [100000, 140, 80,56,44 ,36,  30, 25, 20]*3;
            elseif strcmp(fac,'iaxb')
                sc = [100000, 200, 100, 70, 50, 40, 34, 28, 24]*100;
            end
        end
            
        saveStr = [savePrefix,'_',err_dist,'_',method];

        for r = 1:length(sc)
            fprintf('Effect size %d out of %d \n',r,length(sc))
            if strcmp(fac,'a')||(strcmp(fac,'iaxb')&&~Int_flag) % spurious a TAKE CARE OF THAT LATER
                % gradual increase of main factor
                if r == 1
                    A = [repmat(50,1,b*n/a) zeros(1,b*n/a) repmat(-50,1,b*n/a)]/sc(r); % Group
                    A = zeros(size(A));
                else   
                    A = [repmat(50,1,b*n/a) zeros(1,b*n/a) repmat(-50,1,b*n/a)]/sc(r); % Group
                end
                m_A = mean(A);
                par = sqrt(sum(((A - m_A).^2)./a));
                %                 if subjInt_flag  % TAKE CARE OF THAT LATER
                %                     AS = A.*S/subjInt_scale;
                %                     if Int_flag
                %                         ABS = A.*B.*S/subjInt_scale;
                %                     end
                %                 end
                
            elseif strcmp(fac,'b') ||(strcmp(fac,'iaxb')&&~Int_flag) % spurious b, i.e. both, factor A and B increase in strength, but there's no interaction
                if r == 1
                    B = [50 -50 20 -20];
                    B = repmat(B,1,n);
                    B = zeros(size(B));
                else
                    B = [50 -50 20 -20];
                    B = repmat(B,1,n)/sc(r);
                end 
                m_B = mean(B);
                par = sqrt(sum(((B - m_B).^2)./b));  
                 %                 if subjInt_flag
                %                     BS = B.*S/subjInt_scale;
                %                     if Int_flag
                %                         ABS = A.*B.*S/subjInt_scale;
                %                     end
                %                 end
                
            elseif strcmp(fac,'iaxb') && Int_flag
                % gradual increase of interaction factor
                if r ==1
                    AB = zeros(size(Int));
                else
                    AB = Int/sc(r);
                end
                m_AB = mean(mean(AB));
                par = sqrt(sum(sum(((AB - m_AB).^2)./(a*b))));
%                 if subjInt_flag
%                     ABS = AB.*S/subjInt_scale;
%                 end
            end
            res = zeros(1,Rep);
            
            
            for j = 1:Rep
                % put the model together & add errors
                if Int_flag
                    if ~strcmp(fac,'iaxb') % if fac = 'iaxb', AB has been constructed and scaled above
                    AB = A.*B;
                    end
                    y = A + B + AB;
%                     if subjInt_flag % TAKE CARE OF THAT LATER
%                         y = y + AS + BS + ABS;
%                     end
                else
                    y = A + B;
%                     if subjInt_flag
%                         y = y + AS + BS;
%                     end
                end
                % put the model together & add errors
                if strcmp(err_dist,'exp')
                    y = y + (exprnd(1,1,b*n)).^3;
                elseif strcmp(err_dist,'gauss')
                    y = y + randn(1,b*n);
                end
                
                X = [y;[ones(1,b*n/a) ones(1,b*n/a)*2 ones(1,b*n/a)*3];repmat(1:4,1,n); ...
                    [ones(1,b) ones(1,b)*2 ones(1,b)*3 ones(1,b)*4 ones(1,b)*5 ones(1,b)*6 ones(1,b)*7 ones(1,b)*8 ...
                    ones(1,b)*9 ones(1,b)*10 ones(1,b)*11 ones(1,b)*12]];
                
                switch method   
                    case 'ftest'
                        [SSQs, DFs, MSQs, Fs, Ps] = mixed_b1_w1_anova(X',1);
                        if strcmp(fac,'a')
                            res(j) = Ps{1};
                        elseif strcmp(fac,'b')
                            res(j) = Ps{3};
                        elseif strcmp(fac,'iaxb')
                            res(j) = Ps{4};
                        end
                        
                    case 'raw'
                        c = X(1,:);
                        design = X(2:end,:);
                        exact = 'no';
                        statfun = 'mixedAnova_1b_1w';
                        stat = FtSimLink_mixed_1b_1w_ANOVA(data,neighbours,c,design,statfun,fac,exact);
                        res(j) = stat.prob;
                        
                    case 'res'
                        c = X(1,:);
                        design = X(2:end,:);
                        ncond_a = a;
                        ncond_b = b;
                        
                        for nfac_a = 1:ncond_a
                            for nfac_b = 1:ncond_b
                                idx_ab = X(2,:) == nfac_a & X(3,:) == nfac_b;
                                anovaIn{nfac_a,nfac_b} = c(:,idx_ab);
                            end
                        end
                        
                         if strcmp(fac,'a')||strcmp(fac,'iaxb')
                            for jj = 1:size(anovaIn,2)
                                tmp = zeros(size(anovaIn{1,1}));
                                for ii = 1:size(anovaIn,1)
                                    tmp(:,:,ii) = anovaIn{ii,jj};
                                end
                                bmean(jj) = squeeze(mean(mean(tmp)));
                            end
                        end
                        
                        if strcmp(fac,'b')||strcmp(fac,'iaxb')
                            for jj = 1:size(anovaIn,2)
                                tmp = zeros(size(anovaIn{1,1}));
                                for ii = 1:size(anovaIn,1)
                                    tmp(:,:,ii) = anovaIn{ii,jj};
                                end
                                amean(jj) = squeeze(mean(mean(tmp)));
                            end
                        end
                        
                        if strcmp(fac,'a')
                            for ii = 1:size(anovaIn,1)
                                for jj = 1:size(anovaIn,2)
                                    idx_ab = design(1,:) == ii & design(2,:) == jj;
                                    c_new(idx_ab) = c(idx_ab) - bmean(jj); 
                                end
                            end
                        end
                        
                         if strcmp(fac,'b')
                            for ii = 1:size(anovaIn,1)
                                for jj = 1:size(anovaIn,2)
                                    idx_ab = design(1,:) == ii & design(2,:) == jj;
                                    c_new(idx_ab) = c(idx_ab) - amean(ii); 
                                end
                            end
                        end
                        
                        if strcmp(fac,'iaxb')
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
                        end
                        
                        
                        exact = 'no';
                        statfun = 'mixedAnova_1b_1w';
                        stat = FtSimLink_mixed_1b_1w_ANOVA(data,neighbours,c_new,design,statfun,fac,exact);
                        res(j) = stat.prob;
                        
                    case 'totunres' % ignore subject units, just out of curiosity
                        c = X(1,:);
                        design = X(2:end,:);
                        exact = 'totunres';
                        statfun = 'mixedAnova_1b_1w';
                        % permutation ANOVA is called here
                        stat = FtSimLink_mixed_1b_1w_ANOVA(data,neighbours,c,design,statfun,fac,exact);
                        res(j) = stat.prob;
                        
                    case 'exact' % 
                        if strcmp(fac,'iaxb')
                        warning('No exact test for the interaction effect in a mixed-design ANOVA, using permutation scheme for exact group main effect');  
                        c = X(1,:);
                        design = X(2:end,:);
                        exact = 'yes_a';
                        statfun = 'mixedAnova_1b_1w';
                        stat = FtSimLink_mixed_1b_1w_ANOVA(data,neighbours,c,design,statfun,fac,exact);
                        res(j) = stat.prob;    
                        else
                        c = X(1,:);
                        design = X(2:end,:);
                        exact = 'yes';
                        statfun = 'mixedAnova_1b_1w';
                        stat = FtSimLink_mixed_1b_1w_ANOVA(data,neighbours,c,design,statfun,fac,exact);
                        res(j) = stat.prob;    
                        end
                end
                
            end
            p_val(r) = length(find(res <= 0.05))/Rep;
            param(r) = par;
        end

        switch method
            case 'ftest'
                plot(param,p_val,':d','color',[0.4 0.4 0.4],'linewidth',2)
            case 'raw'
                plot(param,p_val,'g:d','linewidth',2)
            case 'exact'
                plot(param,p_val,':*','linewidth',2)
            case 'res'
                plot(param,p_val,'m:s','linewidth',2)
            case 'totunres'
                plot(param,p_val,'c:<','linewidth',2)
        end
        
        set(gcf,'color','w')
        if strcmp(fac,'a')
            xlabel('\theta_A','FontSize',10,'FontName','Helvetica')
        elseif strcmp(fac,'b')
            xlabel('\theta_B','FontSize',10,'FontName','Helvetica')
        elseif strcmp(fac,'iaxb')
            xlabel('\theta_AB','FontSize',10,'FontName','Helvetica')
        end
        
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
        xlim([min(param) max(param)])
        
        cd(resDir)
        save(saveStr,'param','p_val','A*','B*','sc','method','res')
        clear param p_val c
        cd(stemFolder)  
    end
end


