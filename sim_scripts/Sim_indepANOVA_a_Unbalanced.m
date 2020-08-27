% Monte Carlo simulations independent 2-way ANOVA
% FtSimLink, anova2_cell_mod

% stemFolder = '<INSERT_YOUR_OWN_WORKING_DIR';
stemFolder = '/data/pt_np-helbling/permANOVA/';
addpath([stemFolder 'stat_util/']) 
resDir = [stemFolder 'SimData/ResultsIndepANOVA/']; 

savePrefix = '2way_indep_a_unbalanced';

if ~exist(resDir,'dir')
    mkdir(resDir);
end

Rep = 5;

a = 3;
b = 4;
n = 5;

B = [50 -50 20 -20];
B = reshape(repmat(B,a,n),a,b,n);

design = [repmat([1,2,3],1,b*n);repmat([1 1 1 2 2 2 3 3 3 4 4 4],1,n)];
design = design(:,1:end-1);

load([stemFolder 'dummy'],'data');
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg_neighb, data);

Method_List = {'ftest','raw','exact','res'};
Error_List = {'exp','gauss'};
for ee = 1:length(Error_List)
    figure
    hold on
    for mm = 1:length(Method_List)
        
        method = Method_List{mm};
        err_dist = Error_List{ee};
          
        if strcmp(err_dist,'exp')
            sc = [10000, 40, 20, 12, 10, 8, 6, 4, 3];
        elseif strcmp(err_dist,'gauss')
            sc = [100000, 200, 100, 70, 50, 40, 35, 30, 25]*3; 
        end
        
        
        
        saveStr = [savePrefix,'_',err_dist,'_',method];
        
        for r = 1:length(sc)
            % gradual increase of main factor
            if r ==1
                A = [0; 0;0];
            else
                A = [50; 0;-50]/sc(r);
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
                        case 'ftest'
                            c = reshape(y,1,numel(y));
                            c = c(1:end-1);
                            [p,tbl,stats,terms] = anovan(c, {design(1,:),design(2,:)},'model',2,'display','off');
    
                            res_A(j) = 1 - fcdf(tbl{2,6}, tbl{2,3}, tbl{5,3});
                            res_B(j) = 1 - fcdf(tbl{3,6}, tbl{3,3}, tbl{5,3});
                            res_I(j) = 1 - fcdf(tbl{4,6}, tbl{4,3}, tbl{5,3});
                            
                            
                        case 'raw'
                            fac = 'a';
                            c = reshape(y,1,numel(y));
                            c = c(1:end-1);
                            exact = 'no';
                            statfun = 'indepAnova2way';
                            % permutation ANOVA is called here
                            stat = FtSimLink_indepANOVA(data,neighbours,c,design,statfun,fac,exact);
                            res_A(j) = stat.prob;
                                                        
                        case 'exact'
                            fac = 'a';
                            c = reshape(y,1,a*b*n);
                            c = c(1:end-1);
                            exact = 'yes';
                            statfun = 'indepAnova2way';
                            stat = FtSimLink_indepANOVA(data,neighbours,c,design,statfun,fac,exact);
                            res_A(j) = stat.prob;
                            
                            
                        case 'res'          
                            fac = 'a';
                            c = reshape(y,1,numel(y)); 
                            c = c(1:end-1);
                            ncond_a = a;
                            ncond_b = b;
                            
                            for nfac_a = 1:ncond_a
                                for nfac_b = 1:ncond_b
                                    idx_ab = design(1,:) == nfac_a & design(2,:) == nfac_b;
                                    anovaIn{nfac_a,nfac_b} = c(idx_ab);
                                end
                            end
%                             idx_ab = design(1,:) == ncond_a & design(2,:) == ncond_b;
%                                     anovaIn{ncond_a,ncond_b} = ones(size(c(idx_ab)))*NaN;
%                             
                            for jj = 1:size(anovaIn,2)
                                tmp = zeros(1,max(max(cellfun('length',anovaIn))),size(anovaIn,1));
                                for ii = 1:size(anovaIn,1)
                                    tmp(:,1:length(anovaIn{ii,jj}),ii) = (anovaIn{ii,jj});
                                end
                                bmean(jj) = squeeze(nanmean(nanmean(tmp)));
                            end
                            
                            
                            for ii = 1:size(anovaIn,1)
                                for jj = 1:size(anovaIn,2)
                                    idx_ab = design(1,:) == ii & design(2,:) == jj;
                                    c_new(idx_ab) = c(idx_ab) - bmean(jj);
                                end
                            end
                            
                            exact = 'no';
                            statfun = 'indepAnova2way';
                            stat = FtSimLink_indepANOVA(data,neighbours,c_new,design,statfun,fac,exact);
                            res_A(j) = stat.prob;
                            
                    end
                    
                end
                p_val(r) = length(find(res_A <= 0.05))/Rep;
            paramA(r) = par;
        end
        
        switch method
            case 'ftest'
                plot(paramA,p_val,':d','color',[0.4 0.4 0.4],'linewidth',2)
            case 'raw'
                plot(paramA,p_val,'g:d','linewidth',2)
            case 'exact'
                plot(paramA,p_val,':*','linewidth',2)
            case 'res'
                plot(paramA,p_val,'m:s','linewidth',2)
        end
        title(method)
        cd(resDir)
        save(saveStr,'paramA','p_val*','A','B','sc','method')
        clear paramA p_val* c
    end    
end


