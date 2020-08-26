% Monte Carlo simulations dependent 2-way ANOVA
% Params: fac, Rep, subjInt, subjInt_scale

addpath('<PATH_TO_ADAPTED_FIELDTRIPTOOLBOX>')
ft_defaults

stemFolder = '<INSERT_YOUR_OWN_WORKING_DIR>';
addpath([stemFolder 'stat_util/']) 
resDir = [stemFolder 'SimData/ResultsDepANOVA/']; 

if ~exist(resDir,'dir')
    mkdir(resDir);
end

fac = 'a'; % factor or interaction of interest: 'a', 'b' and 'iaxb' for main factors A and B, and the interaction. 
           % Note: 'b' not tested here yet, as both are within-subject factors. Could be interesting though, as the number of factor levels differs

Int_flag = true; % Interaction between the two within factors?

subjInt_flag = true; % Interaction between within factors and subjects?
subjInt_scale = 100; % strength of interaction

pool_flag = false;
if pool_flag
    statfun = 'depAnova2way_pool_df';
else
    statfun = 'depAnova2way';
end
ASblock_flag = true;

Rep = 10;

Method_List = {'ftest','exact','raw','res','totunres'}; % Note: not all tests are meaningful for each simulation. E.g. it doesn't make sense to calculate an exact test for the interaction. See 
Error_List = {'exp','gauss'};


%% construct save string prefix
savePrefix = sprintf('2way_dep_%s_S',fac);
if subjInt_flag
    savePrefix = [savePrefix sprintf('_subjInt_%d',subjInt_scale)];
end

if ASblock_flag
    savePrefix = [savePrefix '_AS'];
end

if pool_flag
    savePrefix = [savePrefix '_pool'];
end

if Int_flag
    savePrefix = [savePrefix '_AB'];
end

%% build data matrix
a = 3; % number of levels factor A
b = 4; % number of levels factor B
n = 5; % number of UOs (e.g. subjects)

S = [repmat(-50,a,b),repmat(-20,a,b),zeros(a,b),repmat(20,a,b),repmat(50,a,b)];
S = reshape(S,a,b,n);

A = [50; 0; -50]; % see Anderson et al., 2001
B = [50 -50 20 -20]; % 
if Int_flag
    Int = A*B;
end
A = reshape(repmat(A,b,n),a,b,n);
B = reshape(repmat(B,a,n),a,b,n);

if subjInt_flag
    AS = A.*S/subjInt_scale;
    BS = B.*S/subjInt_scale;
end

%% build design matrix
design = [repmat([1,2,3],1,b*n);repmat([1 1 1 2 2 2 3 3 3 4 4 4],1,n); [ones(1,a*b) ones(1,a*b)*2 ones(1,a*b)*3 ones(1,a*b)*4 ones(1,a*b)*5]];
if ASblock_flag
    ASvar = zeros(1,size(design,2));
    for j = 1:n
        for ii = 1:a
            ASvar(a*b*(j-1)+(ii:a:a*b)) = j*a-3+ii;
        end
    end
    design = [design;ASvar];
    
    rearrange = [];
    for rr = 1:n
        for ii = 1:a
            rearrange = [rearrange, (a*b*(rr-1))+ii:a:a*b*rr];
        end
    end
    
    design = design(:,rearrange);    
end

% Fieldtrip wants a neighbourhood matrix for its clustering approach
load([stemFolder filesep 'dummy'],'data');
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg_neighb, data);

sprintf('Working on %s\n',savePrefix)

for ee = 1:length(Error_List)
    fprintf('Error dist.: %s \n',Error_List{ee})
            err_dist = Error_List{ee};
     if strcmp(err_dist,'exp')
            if strcmp(fac, 'a')
                sc = [100000, 230, 140, 90, 65, 48, 34,26, 20]/10;
            elseif strcmp(fac,'iaxb')
                sc = [100000, 400, 200, 130, 80, 55, 40, 30, 24]*3;
            end
        elseif strcmp(err_dist,'gauss')
            if strcmp(fac, 'a')
                sc = [100000, 140, 80, 56, 44, 36, 30, 25, 20]*3;
            elseif strcmp(fac,'iaxb')
                sc = [100000, 200, 100, 70, 50, 40, 34, 28, 24]*100;
            end
     end
        
    figure
    hold on
    for mm = 1:length(Method_List)
        fprintf('Stats: %s \n',Method_List{mm})
        method = Method_List{mm};
  
        saveStr = [savePrefix,'_',err_dist,'_',method];
        p_val = zeros(1,length(sc));
        params = zeros(1,length(sc));
        
        for r = 1:length(sc) % loop through increasing effect sizes
            fprintf('Effect size %d out of %d \n',r,length(sc))
            if strcmp(fac,'a')||(strcmp(fac,'iaxb')&&~Int_flag)
                % gradual increase of main factor
                if r ==1
                    A = [0; 0;0];
                else
                    A = [50; 0;-50]/sc(r);
                end
                m_A = mean(A);
                par = sqrt(sum(((A - m_A).^2)./a));
                A = reshape(repmat(A,b,n),a,b,n);
                if Int_flag
                    AB = A.*B;
                end
                if subjInt_flag
                    AS = A.*S/subjInt_scale;
                    if Int_flag
                        ABS = AB.*S/subjInt_scale;
                    end
                end
                
            elseif strcmp(fac,'iaxb') && Int_flag
                % gradual increase of interaction factor
                if r ==1
                    AB = zeros(size(Int));
                else
                    AB = Int/sc(r);
                end
                m_AB = mean(mean(AB));
                par = sqrt(sum(sum(((AB - m_AB).^2)./(a*b))));
                AB = reshape(repmat(AB,1,n),a,b,n);
                if subjInt_flag
                    ABS = AB.*S/subjInt_scale;
                end
            end
            res = zeros(1,Rep);
            
           for j = 1:Rep % loop through repetitions
                % put the model together & add errors
                if Int_flag
                    y = S + A + B + AB;
                    if subjInt_flag
                        y = y + AS + BS + ABS; % add subject-interaction terms
                    end
                else
                    y = S + A + B;
                     if subjInt_flag
                        y = y + AS + BS;
                    end
                end
                
                % add error terms
                if strcmp(err_dist,'exp')
                    y = y + (exprnd(1,a,b,n)).^3;
                elseif strcmp(err_dist,'gauss')
                    y = y + randn(a,b,n);
                end
                
                switch method
                    case 'ftest'
                        
                        % prepare data input for anova2rm_cell_mod*
                        % functions which are called here directly:
                        % 3 x 4 cell matrix, each containing 5 data points
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
                      
                        if pool_flag
                            [FA, FB, FI, dfa, dfb, dfi] = anova2rm_cell_mod_pool_df(c); % adapted from the resampling statistic toolbox
                        else
                            [FA, FB, FI, dfa, dfb, dfi] = anova2rm_cell_mod(c);
                        end
                       
                        if strcmp(fac,'a')
                            res(j) = 1 - fcdf(FA, dfa(1), dfa(2));
                        elseif strcmp(fac,'b')
                            res(j) = 1 - fcdf(FB, dfb(1), dfb(2));
                        elseif strcmp(fac,'iaxb')
                            res(j) = 1 - fcdf(FI, dfi(1), dfi(2));
                        end
                        
                    case 'raw'
                        c = reshape(y,1,a*b*n); % concatenate data points into a vector. Rearrangment into cell matrix will be done within the statfun (e.g.ft_statfun_depAnova2way)
                        exact = 'no';
                        % permutation ANOVA is called here
                        if ~ASblock_flag
                            stat = FtSimLink_depANOVA(data,neighbours,c,design,statfun,fac,exact);
                        else
                            stat = FtSimLink_depANOVA_AS(data,neighbours,c(rearrange),design,statfun,fac,exact);
                        end
                        res(j) = stat.prob;
                        
                    case 'exact'
                        c = reshape(y,1,a*b*n);
                        exact = 'yes';
                        if ~ASblock_flag
                            stat = FtSimLink_depANOVA(data,neighbours,c,design,statfun,fac,exact);
                        else
                            stat = FtSimLink_depANOVA_AS(data,neighbours,c(rearrange),design,statfun,fac,exact);
                        end
                        res(j) = stat.prob;
                        
                    case 'totunres'
                        c = reshape(y,1,a*b*n);
                        exact = 'totunres';
                        if ~ASblock_flag
                            stat = FtSimLink_depANOVA(data,neighbours,c,design,statfun,fac,exact);
                        else
                            stat = FtSimLink_depANOVA_AS(data,neighbours,c(rearrange),design,statfun,fac,exact);
                        end
                        res(j) = stat.prob;
                        
                    case 'res'
                        
                        c = reshape(y,1,a*b*n);
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
                        if ~ASblock_flag
                            stat = FtSimLink_depANOVA(data,neighbours,c,design,statfun,fac,exact);
                        else
                            stat = FtSimLink_depANOVA_AS(data,neighbours,c(rearrange),design,statfun,fac,exact);
                        end
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
        save(saveStr,'param','p_val','S*','A*','B*','sc','method','res')
        clear param p_val c
        cd(stemFolder)
        
    end
end

