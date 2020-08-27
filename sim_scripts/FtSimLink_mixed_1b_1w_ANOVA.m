function stat = FtSimLink_mixed_1b_1w_ANOVA(data,neighbours,y,design,statfun,fac,main_exact_flag)

origDir = pwd;
data.powspctrm =  y'; % data based on a dummy structure

% specify configuration structure
cfg = [];
cfg.method           = 'montecarlo';
cfg.channel          = {'MLC11'};
cfg.latency          = [1 1];
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 1; % show all clusters regardless of significance 
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 0; % only one channel in dummy data set
cfg.tail             = 1; % F-values can only take positive values
cfg.clustertail      = 1;
cfg.frequency        = [4 4];
cfg.time             = [1 1];
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'yes';
cfg.avgoverchan      = 'no';
cfg.neighbours       = neighbours;
cfg.numrandomization = 999; % number of permutations

cfg.statistic        = statfun;
cfg.fac              = fac;

% define permutations to be done
switch cfg.fac
    case 'a'
        switch main_exact_flag
            case 'yes'
                cfg.ivar = [1 2];
                cfg.wvar = 3;
                cfg.uvar = [];
                cfg.cvar = [];
            case 'totunres'
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = [];
                cfg.cvar = [];
            case 'no'
                [~, ftpath] = ft_version;
                cd([ftpath '/private'])
                clear resample
                % resample group main effect levels
                cfg.ivar = [1 2];
                cfg.wvar = 3;
                cfg.uvar = [];
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                resample_g = resampledesign(cfg, design);
                % resample within main effect levels
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 3;
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                resample_w = resampledesign(cfg, design);
                
                % concatenate permutations
                for j = 1:size(resample_g,1)
                    resample(j,:) = resample_g(j,resample_w(j,:));
                end
                
                cfg.resample = resample;
                cd(origDir)
        end
    case 'b'
        switch main_exact_flag
            case 'yes'
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 3;
                cfg.cvar = [];
            case 'totunres'
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = [];
                cfg.cvar = [];
            case 'no'
                [~, ftpath] = ft_version;
                cd([ftpath '/private'])
                clear resample
                % resample group main effect levels
                cfg.ivar = [1 2];
                cfg.wvar = 3;
                cfg.uvar = [];
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                % resample within main effect levels
                resample_g = resampledesign(cfg, design);
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 3;
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                resample_w = resampledesign(cfg, design);
                
                for j = 1:size(resample_g,1)
                    resample(j,:) = resample_g(j,resample_w(j,:));
                end
                
                cfg.resample = resample;
                cd(origDir)       
        end
               
    case 'iaxb'
        switch main_exact_flag
            case 'no'
                [~, ftpath] = ft_version;
                cd([ftpath '/private'])
                clear resample
                % resample group main effect levels
                cfg.ivar = [1 2];
                cfg.wvar = 3;
                cfg.uvar = [];
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                resample_g = resampledesign(cfg, design);
                % resample within main effect levels
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 3;
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                resample_w = resampledesign(cfg, design);
                
                for j = 1:size(resample_g,1)
                    resample(j,:) = resample_g(j,resample_w(j,:));
                end
                
                cfg.resample = resample;
                cd(origDir)
            case 'totunres'
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = [];
                cfg.cvar = [];
            case 'yes_a'
                cfg.ivar = [1 2];
                cfg.wvar = 3;
                cfg.uvar = [];
                cfg.cvar = [];
            case 'yes_b'
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 3;
                cfg.cvar = [];
        end     
end

cfg.design = design;
stat = ft_freqstatistics(cfg, data);