function [s,cfg] = ft_statfun_depAnova2way_pool_df(cfg, dat, design)

% FT_STATFUN_DEPANOVA2WAY calculates the repeated-measurements samples 2-way ANOVA
% F-statistic(s)
% on the biological data in dat (the dependent variable), using the information on
% the independent variables (iv) and the uvar variable indicating the subjects in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'depAnova2way'

% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_depAnova2way(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (iv),  Nfac x Nreplications
%
% Configuration options:
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or
%   cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail = -1, 0, or 1, left, two-sided, or right (default=1)
%              cfg.tail in combination with cfg.computecritval='yes'
%              determines whether the critical value is computed at
%              quantile cfg.alpha (with cfg.tail=-1), at quantiles
%              cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%              quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification:
%   cfg.ivar        = row number(s) of the design that contain(s) the labels of the conditions that must be
%                     compared (default= [1,2]). The labels range from 1 to the number of conditions.
%   cfg.uvar        = row number of design that contains the labels of the UOs (subjects or trials)
%                     (default=3). The labels are assumed to be integers ranging from 1 to 
%                     the number of UOs.
% Optionally: 
%   cfg.cvar        = row number of the design that contains the labels of the conditions of the other factor, 
%                     when testing for the main effect of one factor. 


% Adapted by Saskia Helbling from Eric Maris' statfun_depsamplesT
% Repeated measurement ANOVA is calculated by means of the function
% ANOVA2RM_CELL_MOD based on ANOVA2RM_CELL from the Resampling statistical toolkit by Arnaud Delorme 
% http://www.mathworks.com/matlabcentral/fileexchange/27960-resampling-statistical-toolkit/content/statistics/anova2rm_cell.m
% published under the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


% set the defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end

%% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end


% only 2 way factorial ANOVA for now
ncond_a = length(unique(design(cfg.ivar(1),:)));
ncond_b = length(unique(design(cfg.ivar(2),:)));

ncond = ncond_a*ncond_b;
nrepl=zeros(ncond_a,ncond_b);

for condindx_a = 1:ncond_a
    for condindx_b = 1:ncond_b
        nrepl(condindx_a,condindx_b) = nrepl(condindx_a,condindx_b)+length(find(design(cfg.ivar(1),:)==condindx_a& design(cfg.ivar(2),:)==condindx_b));
    end
end
if min(nrepl)~=max(nrepl)
    error('Use only balanced cell designs for 2-way repeated measurement ANOVA!'); 
end

if sum(sum(nrepl))<size(design,2)
    error('Invalid specification of the independent variable in the design array.');
end
if sum(sum(nrepl))<=ncond
    error('The must be more trials/subjects than levels of the independent variable.');
end

nsmpls = size(dat,1);

%% degrees of freedom
dfA     = ncond_a-1;
dfB     = ncond_b-1; 
dfAB    = (ncond_a-1)*(ncond_b-1); 
dfABS   = (ncond_a-1)*(ncond_b-1)*(nrepl(1)-1);

dfa  = [dfA dfABS];
dfb  = [dfB dfABS];
dfi =  [dfAB dfABS];

%% compute the statistic
if strcmp(cfg.computestat, 'yes')
    % compute the statistic
    s.stat=zeros(nsmpls,1);
    for nfac_a = 1:ncond_a
        for nfac_b = 1:ncond_b
            idx_ab = design(cfg.ivar(1),:) == nfac_a & design(cfg.ivar(2),:) == nfac_b;
            anovaIn{nfac_a,nfac_b} = dat(:,idx_ab);
        end
    end

    [FA, FB, FI] = anova2rm_cell_mod(anovaIn);

    switch cfg.fac
        case 'a'
            s.stat=FA;
        case 'b'
            s.stat=FB;
        case 'iaxb'
            s.stat=FI;
    end
end

%% compute the critical value
if strcmp(cfg.computecritval,'yes')
    % also compute the critical values
  if cfg.tail==-1
      error('For a dependent samples F-statistic, it does not make sense to calculate a left tail critical value.');
  end
  if cfg.tail==0
      error('For a dependent samples F-statistic, it does not make sense to calculate a two-sided critical value.');
  end
    if cfg.tail==1
        switch cfg.fac
            case 'a'
                s.critval = finv(1-cfg.alpha, dfa(1), dfa(2));
            case 'b'
                s.critval = finv(1-cfg.alpha, dfb(1), dfb(2));
            case 'iaxb'
                s.critval = finv(1-cfg.alpha, dfi(1), dfi(2));
        end
    end
end


%% Compute probabilities 
if strcmp(cfg.computeprob,'yes')

   if cfg.tail==-1
      error('For a dependent samples F-statistic, it does not make sense to calculate a left p-value.');
   end
  if cfg.tail==0
      error('For a dependent samples F-statistic, it does not make sense to calculate a two-sided p-value.');
  end
    
    
    if cfg.tail==1      
        switch cfg.fac
            case 'a'
                s.prob = 1 - fcdf(s.stat, dfa(1), dfa(2));
            case 'b'
                s.prob = 1 - fcdf(s.stat, dfb(1), dfb(2));
            case 'iaxb'
                s.prob = 1 - fcdf(s.stat, dfi(1), dfi(2));
        end
    end
end

