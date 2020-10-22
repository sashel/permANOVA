function [s,cfg] = ft_statfun_mixedAnova_1b_1w(cfg, dat, design)

% FT_STATFUN_MIXEDANOVA_1B_1W calculates the F-statistic for a mixed design
% with one within subjects and one between subjects variable
% on the biological data in dat (the dependent variable), using the information on
% the independent variables (iv) on the subjects indicated by uvar in the design matrix.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'mixedEff_1b_1w'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_mixedEff_1b_1w(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variables (iv),  Nfac x Nreplications &
%   the uvar indicating the subjects
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
%   cfg.ivar        = row number of the design that contains the labels of the conditions that must be
%                        compared (default=1). The labels range from 1 to the number of conditions.

%   cfg.uvar        = row number of design that contains the labels of the UOs (subjects or trials)
%                     (default=3). The labels are assumed to be integers ranging from 1 to 
%                     the number of UOs.
%   or cfg.war      = row number of design that contains the labels of the UOs (subjects or trials)
%                     (default=3), equal labels should be kept together. The labels are assumed to be integers ranging from 1 to 
%                     the number of UOs.

% Optionally: 
%   cfg.cvar        = row number of the design that contains the labels of the conditions of the other factor, 
%                     when testing for the main effect of one factor. 



% Adapted by Saskia Helbling from Eric Maris' statfun_depsamplesT
% Uses MIXED_B1_W1_ANOVA


% set the defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;

if ((~isfield(cfg,'uvar') || isempty(cfg.uvar))) && ((~isfield(cfg,'wvar') || isempty(cfg.wvar)))
    error('uvar or wvar must be specified for a mixed within between factorial design');
end

if (isfield(cfg,'uvar') && ~isempty(cfg.uvar)) && (isfield(cfg,'wvar') && ~isempty(cfg.wvar))
    error('Only uvar OR wvar, not both, must be specified for a mixed within between ANOVA');
end

if (~isfield(cfg,'ivar')) || ~isempty(setdiff(cfg.ivar,[1 2]))
    error('Please make sure that levels of the two factors (ivar) are given in the first two rows of the design matrix!');
end

if (isfield(cfg,'uvar') && ~isempty(setdiff(cfg.uvar,3)))||(isfield(cfg,'wvar') && ~isempty(setdiff(cfg.wvar,3)))
    error('Please make sure that subject indices (uvar or wvar) are given in the third & final row of the design matrix!');
end

if size(design,1)~=3
    error('For a mixed design 2-way ANOVA, the design matrix must have 3 rows: 1. group factor (ivar), 2. within factor (ivar), 3. subject indices (uvar or wvar)');
end


nsmpls = size(dat,1); 
s.stat=zeros(nsmpls,1);

    if strcmp(cfg.computestat, 'yes')
        % compute the statistic
        
        [SSQs, DFs, MSQs, Fs, Ps] = mixed_b1_w1_anova([dat' design'],1);
        
        
        switch cfg.fac
            case 'a' % between
                s.stat=Fs{1};
            case 'b'
                s.stat=Fs{3}; % within
            case 'iaxb'
                s.stat=Fs{4};
        end
    end
    
    if strcmp(cfg.computecritval,'yes')
        % also compute the critical values
        [SSQs, DFs, MSQs, Fs, Ps] = mixed_b1_w1_anova([dat' design'],1);
        if cfg.tail==-1
            error('For an independent samples F-statistic, it does not make sense to calculate a left tail critical value.');
        end;
        if cfg.tail==0
            error('For an independent samples F-statistic, it does not make sense to calculate a two-sided critical value.');
        end;
        if cfg.tail==1
            switch cfg.fac
                case 'a'
                    s.critval = finv(1-cfg.alpha, DFs{1}, DFs{2});
                case 'b'
                    s.critval = finv(1-cfg.alpha, DFs{3}, DFs{5});
                case 'iaxb'
                    s.critval = finv(1-cfg.alpha, DFs{4}, DFs{5});
            end
        end;
    end
    
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
                    s.prob = 1 - fcdf(s.stat, DFs{1}, DFs{2});
                case 'b'
                    s.prob = 1 - fcdf(s.stat, DFs{3}, DFs{5});
                case 'iaxb'
                    s.prob = 1 - fcdf(s.stat, DFs{4}, DFs{5});
            end
        end;
    end

