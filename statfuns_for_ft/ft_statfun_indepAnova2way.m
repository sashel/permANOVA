function [s,cfg] = ft_statfun_indepAnova2way(cfg, dat, design)

% FT_STATFUN_INDEPANOVA2WAY calculates the independent samples 2-way ANOVA
% F-statistic(s) on the biological data in dat (the dependent variable), using the information on
% the independent variables (iv) in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'indepAnova2way'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_indepAnova2way(cfg, dat, design);
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
% cfg.computeprob='yes'.
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

% Adapted by Saskia Helbling from Eric Maris' statfun_depsamplesT in the Fieldtrip toolbox

% set the defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

%% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    error('cfg.uvar should not exist for an independent samples statistic');
end

nsmpls = size(dat,1);

%% compute the statistic
if strcmp(cfg.computestat, 'yes')
    % compute the statistic
    s.stat=zeros(nsmpls,1);

    [p,tbl,stats,terms] = anovan(dat, {design(1,:),design(2,:)},'model',2,'display','off');

    switch cfg.fac
        case 'a'
            s.stat=tbl{2,6};
        case 'b'
            s.stat=tbl{3,6};
        case 'iaxb'
            s.stat=tbl{4,6};
    end
end

%% compute the critical value
if strcmp(cfg.computecritval,'yes')
    % also compute the critical values
 [p,tbl,stats,terms] = anovan(dat, {design(1,:),design(2,:)},'model',2,'display','off');

    if cfg.tail==-1
        error('For an independent samples F-statistic, it does not make sense to calculate a left tail critical value.');
    end;
    if cfg.tail==0
        error('For an independent samples F-statistic, it does not make sense to calculate a two-sided critical value.');
    end;
    if cfg.tail==1
        switch cfg.fac
            case 'a'
                s.critval = finv(1-cfg.alpha, tbl{2,3}, tbl{5,3});
            case 'b'
                s.critval = finv(1-cfg.alpha, tbl{3,3}, tbl{5,3});
            case 'iaxb'
                s.critval = finv(1-cfg.alpha, tbl{4,3}, tbl{5,3});
        end
    end;
end


%% Compute probabilities 
if strcmp(cfg.computeprob,'yes')

    if cfg.tail==-1
        error('For an independent samples F-statistic, it does not make sense to calculate a left tail p-value.');
    end;
    if cfg.tail==0
        error('For an independent samples F-statistic, it does not make sense to calculate a two-sided p-value.');
    end;
    if cfg.tail==1
        
        switch cfg.fac
            case 'a'
                s.prob = 1 - fcdf(s.stat, tbl{2,3}, tbl{5,3});
            case 'b'
                s.prob = 1 - fcdf(s.stat, tbl{3,3}, tbl{5,3});
            case 'iaxb'
                s.prob = 1 - fcdf(s.stat, tbl{4,3}, tbl{5,3});
        end
    end;
end
