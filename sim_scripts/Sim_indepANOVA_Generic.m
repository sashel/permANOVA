% Monte Carlo simulations independent 2-way ANOVA
% Uses FtSimLink, anova2_cell_mod
% Params: fac, Rep, int_flag, unbalanced_flag

% stemFolder = '<INSERT_YOUR_OWN_WORKING_DIR>';
stemFolder = '/data/pt_np-helbling/permANOVA/';
addpath([stemFolder 'stat_util/'])
resDir = [stemFolder 'SimData/ResultsIndepANOVA/'];

if ~exist(resDir,'dir')
    mkdir(resDir);
end

fac = 'a'; % factor or interaction of interest: 'a' and 'iaxb' for main factor A and the interaction.
% Note: 'b' not tested yet, as 'a'and 'b' both are within-subject factors
Rep = 1000; % reduce number of repetitions REP when testing the code
int_flag = false; % Interaction between the two within factors?
unbalanced_flag = false;

Method_List = {'ftest','exact','raw','res','totunres'}; % Note: not all tests are meaningful for each simulation. E.g. it doesn't make sense to calculate an exact test for the interaction. See
Error_List = {'exp','gauss'};

%% construct save string prefix
savePrefix = sprintf('2way_indep_%s',fac);

if unbalanced_flag
    savePrefix = [savePrefix '_unbalanced'];
end

if int_flag
    savePrefix = [savePrefix '_AB'];
end

%% build data and design matrices
a = 3;
b = 4;
n = 5;

A = [50; 0; -50]; % see Anderson et al., 2001
B = [50 -50 20 -20]; %
if int_flag
    Int = A*B;
end
A = reshape(repmat(A,b,n),a,b,n);
B = reshape(repmat(B,a,n),a,b,n);

design = [repmat([1,2,3],1,b*n);repmat([1 1 1 2 2 2 3 3 3 4 4 4],1,n)];
if unbalanced_flag
    design = design(:,1:end-1);
end

%% define neighbourhood structure for clustering (Fieldtrip functions expects it as an input even if we are performing cluster statistics on a single sensor)
load([stemFolder 'dummy'],'data');
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg_neighb, data);

%% loop across error types and permutation strategies
for ee = 1:length(Error_List)
    figure
    hold on
    err_dist = Error_List{ee};
    % define effect sizes such that they nicely cover the full range of the
    % power curves
    if strcmp(err_dist,'exp')
        if strcmp(fac, 'a')
            sc = [10000, 40, 20, 12, 10, 8, 6, 4, 3];
        elseif strcmp(fac,'iaxb')
            sc = [100000, 200, 100, 80, 60, 50, 40, 35, 30]*3;
        end
    elseif strcmp(err_dist,'gauss')
        if strcmp(fac, 'a')
            sc = [100000, 200, 100, 70, 50, 40, 35, 30, 25]*3;
        elseif strcmp(fac,'iaxb')
            sc = [100000, 200, 100, 70, 50, 40, 30, 25, 20]*100;
        end
    end
    
    for mm = 1:length(Method_List)
        method = Method_List{mm};
        saveStr = [savePrefix,'_',err_dist,'_',method]; % append save string
        p_val = zeros(1,length(sc)); % permutation p-value across effect sizes
        params = zeros(1,length(sc));
        
        for r = 1:length(sc) % loop through increasing effect sizes
            fprintf('Effect size %d out of %d \n',r,length(sc))
            if strcmp(fac,'a')||(strcmp(fac,'iaxb')&&~int_flag) % increase effect size of factor A in the absence of an interaction
                % gradual increase of main factor
                if r ==1
                    A = [0; 0; 0];
                else
                    A = [50; 0; -50]/sc(r);
                end
                m_A = mean(A);
                par = sqrt(sum(((A - m_A).^2)./a));
                A = reshape(repmat(A,b,n),a,b,n);
                if int_flag
                    AB = A.*B;
                end
            elseif strcmp(fac,'iaxb') && int_flag % increase interaction effect
                if r ==1
                    AB = zeros(size(Int));
                else
                    AB = Int/sc(r);
                end
                m_AB = mean(mean(AB));
                par = sqrt(sum(sum(((AB - m_AB).^2)./(a*b))));
                AB = reshape(repmat(AB,1,n),a,b,n);
            end
            
            for j = 1:Rep              
                % put the model together & add errors
                if int_flag
                    y = A + B + AB;
                else
                    y = A + B;
                end
                
                % add error terms
                if strcmp(err_dist,'exp')
                    y = y + (exprnd(1,a,b,n)).^3;
                elseif strcmp(err_dist,'gauss')
                    y = y + randn(a,b,n);
                end
                
                c = reshape(y,1,numel(y));
                if unbalanced_flag
                    c = c(1:end-1);
                end
                
                switch method
                    case 'ftest'                       
                        [p,tbl,stats,terms] = anovan(c, {design(1,:),design(2,:)},'model',2,'display','off');
                        
                        if strcmp(fac,'a')
                            res(j) = 1 - fcdf(tbl{2,6}, tbl{2,3}, tbl{5,3});
                        elseif strcmp(fac,'b')
                            res(j) = 1 - fcdf(tbl{3,6}, tbl{3,3}, tbl{5,3});
                        elseif strcmp(fac,'iaxb')
                            res(j) = 1 - fcdf(tbl{4,6}, tbl{4,3}, tbl{5,3});
                        end
                        
                    case 'raw'
                        exact = 'no';
                        statfun = 'indepAnova2way';
                        % call permutation ANOVA via linking function
                        stat = FtSimLink_indepANOVA(data,neighbours,c,design,statfun,fac,exact);
                        res(j) = stat.prob;
                        
                    case 'exact'
                        exact = 'yes';
                        statfun = 'indepAnova2way';
                        stat = FtSimLink_indepANOVA(data,neighbours,c,design,statfun,fac,exact);
                        res(j) = stat.prob;
                        
                    case 'res'
                        ncond_a = a;
                        ncond_b = b;
                        
                        for nfac_a = 1:ncond_a
                            for nfac_b = 1:ncond_b
                                idx_ab = design(1,:) == nfac_a & design(2,:) == nfac_b;
                                anovaIn{nfac_a,nfac_b} = c(idx_ab);
                            end
                        end
                        
                        if strcmp(fac,'a')||strcmp(fac,'iaxb')
                            for jj = 1:size(anovaIn,2)
                                tmp = zeros(max(max(cellfun('length',anovaIn))),size(anovaIn,1))*NaN;
                                for ii = 1:size(anovaIn,1)
                                    tmp(1:length(anovaIn{ii,jj}),ii) = (anovaIn{ii,jj});
                                end
                                bmean(jj) = squeeze(nanmean(nanmean(tmp)));
                            end
                        end
                        
                        if strcmp(fac,'iaxb')
                            for ii = 1:size(anovaIn,1)
                                tmp = zeros(max(max(cellfun('length',anovaIn))),size(anovaIn,2))*NaN;
                                for jj = 1:size(anovaIn,2)
                                    tmp(1:length(anovaIn{ii,jj}),jj) = (anovaIn{ii,jj});
                                end
                                amean(ii) = squeeze(nanmean(nanmean(tmp)));
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
                        
                        if strcmp(fac,'iaxb')
                            tmp = zeros(max(max(cellfun('length',anovaIn))),size(anovaIn,1),size(anovaIn,2))*NaN;
                            for ii = 1:size(anovaIn,1)
                                for jj = 1:size(anovaIn,2)
                                    tmp(1:length(anovaIn{ii,jj}),ii,jj) = anovaIn{ii,jj};
                                end
                                
                            end
                            totmean = squeeze(nanmean(nanmean(nanmean(tmp))));
                            
                            for ii = 1:size(anovaIn,1)
                                for jj = 1:size(anovaIn,2)
                                    idx_ab = design(1,:) == ii & design(2,:) == jj;
                                    c_new(idx_ab) = c(idx_ab) - bmean(jj) - amean(ii) + totmean;
                                    
                                end
                            end
                        end
                        
                        exact = 'no';
                        statfun = 'indepAnova2way';
                        stat = FtSimLink_indepANOVA(data,neighbours,c_new,design,statfun,fac,exact);
                        res(j) = stat.prob;
                        
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
        end
        title(method)
        cd(resDir)
        save(saveStr,'param','p_val*','A','B','sc','method')
        clear param p_val* c
    end
end


