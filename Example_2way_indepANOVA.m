stemFolder = pwd;
load([stemFolder filesep 'dummy'],'data');

cfg_neighb = [];
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg_neighb, data);
   
%% Gender-Dosage test data set
% Test data set from: http://personality-project.org/r/r.guide/r.anova.html

% Indep 2x2 ANOVA: factors Gender & dosage on alertness level

% Design matrix for the test data
design(1,:) = [ones(1,2*4) ones(1,2*4)*2]; % First factor: Gender
design(2,:) = repmat([ones(1,4) ones(1,4)*2],1,2); % Second Factor: Dosage

data.powspctrm = [8 12 13 12 6 7 23 14 15 12 22 14 15 12 18 22]';

% Results given at the personality-project.org page

%               Df  Sum Sq Mean Sq F value Pr(>F)
% Gender         1  76.562  76.562  2.9518 0.1115
% Dosage         1   5.062   5.062  0.1952 0.6665
% Gender:Dosage  1   0.063   0.063  0.0024 0.9617
% Residuals     12 311.250  25.938     
             


%% Specify configuration structure: the dummy data structure consists of a single channel, time bin and frequency of interest
cfg = [];

cfg.channel          = {'MLC11'};
cfg.latency          = [1 1];
cfg.frequency        = [4 4];
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
cfg.avgoverchan      = 'no';

cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 1; % consider all clusters irrespectively of significance
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 0; % we are looking at one channel only

cfg.statistic        = 'indepAnova2way';
cfg.fac              = 'iaxb'; % main effect Gender: 'a', main effect Dosage: 'b' and interaction: 'iaxb'
cfg.tail             = 1; % F-values can only take positive values
cfg.clustertail      = 1; 
cfg.design           = design;
cfg.neighbours       = neighbours;
cfg.numrandomization = 999; % number of permutations



 %% define permutation strategies
 main_exact_flag = 'yes';
            switch cfg.fac
                case 'a' 
                    switch main_exact_flag
                        case 'no' % unrestricted permutations across levels of both factors
                    cfg.ivar = [1 2]; 
                    cfg.wvar = [];
                    cfg.uvar = [];
                    cfg.cvar = [];
                        case 'yes'  % restricted permutations within levels of the other factor (exact test)
                    cfg.ivar = [1 2];
                    cfg.wvar = [];
                    cfg.uvar = [];
                    cfg.cvar = 2;
                    end
                    
                case 'b' 
                   switch main_exact_flag
                        case 'no'
                    cfg.ivar = [1 2];
                    cfg.wvar = [];
                    cfg.uvar = [];
                    cfg.cvar = [];
                        case 'yes'
                    cfg.ivar = [1 2]; 
                    cfg.wvar = [];
                    cfg.uvar = [];
                    cfg.cvar = 1;
                    end

                case 'iaxb' % unrestricted permutations across levels of both factors
                    cfg.ivar = [1 2]; 
                    cfg.wvar = [];
                    cfg.uvar = [];
                    cfg.cvar = [];
            end

stat = ft_freqstatistics(cfg, data);
