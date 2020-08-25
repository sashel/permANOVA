function stat = FtSimLink_indepANOVA(data,neighbours,y,design,statfun,fac,main_exact_flag)

data.powspctrm =  y'; % data based on a dummy structure

% Specify configuration structure
cfg = [];

cfg.channel          = {'MLC11'};
cfg.latency          = [1 1];
cfg.frequency        = [4 4];
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
cfg.avgoverchan      = 'no';

cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % threshold for a data entry being a cluster candidate
cfg.clusteralpha     = 1; % show all clusters irrespectively of significance
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 0; % the dummy data has only one channel

cfg.statistic        = statfun;
cfg.tail             = 1; % F-values can only take positive values
cfg.clustertail      = 1; 
cfg.design           = design;
cfg.neighbours       = neighbours;
cfg.numrandomization = 999; % number of permutations
cfg.fac              = fac;

 % define permutations to be performed
            switch cfg.fac
                case 'a' % exact main a
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
                    cfg.cvar = 2;
                    end
                    
                case 'b' % exact main b
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

                case 'iaxb' % a x b interaction
                    cfg.ivar = [1 2];
                    cfg.wvar = [];
                    cfg.uvar = [];
                    cfg.cvar = [];
            end


cfg.design = design';
stat = ft_freqstatistics(cfg, data);
