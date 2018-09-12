function stat = FtSimLink_depANOVA_AS(data,neighbours,y,design,statfun,fac,main_exact_flag)

data.powspctrm =  y'; % data from dummy data set structure

%% Specify configuration structure based on dummy data set
cfg = [];
cfg.channel          = {'MLC11'};
cfg.latency          = [1 1];
cfg.frequency        = [4 4];
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
cfg.avgoverchan      = 'no';

cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % threshold for a data bin constituting a 'cluster' candidate (clusters here are based on a single data entry)
cfg.clusteralpha     = 1; % return all 'clusters' irrespectively of significance
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 0; % we are looking at only one channel

cfg.statistic        = statfun;
cfg.tail             = 1; % F-values can only take positive values
cfg.clustertail      = 1;
cfg.design           = design;
cfg.neighbours       = neighbours;
cfg.numrandomization = 999;
cfg.fac              = fac;

[~, ftpath] = ft_version;
origDir = pwd;

% define permutations to be done
switch cfg.fac
    case 'a' % exact main a
        switch main_exact_flag
            case 'no'
                cd([ftpath '/private'])
                clear resample
                % resampling group main
                cfg.ivar = [1 2];
                cfg.wvar = 4;
                cfg.uvar = 3;
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                [resample_s] = resampledesign(cfg, design);
                
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 4;
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                [resample_b] = resampledesign(cfg, design);
                
                for j = 1:size(resample_s,1)
                    resample(j,:) = resample_s(j,resample_b(j,:));
                end
                cfg.resample = resample;
                cd(origDir)
                
            case 'totunres'
                cd([ftpath '/private'])
                clear resample
                % resampling group main
                cfg.ivar = [1 2];
                cfg.wvar = 4;
                cfg.uvar = [];
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                [resample_s] = resampledesign(cfg, design);
                
                cfg.ivar = [1 2];
                cfg.wvar = [];
                cfg.uvar = 4;
                cfg.cvar = [];
                cfg.resampling = 'permutation';
                [resample_b] = resampledesign(cfg, design);
                
                for j = 1:size(resample_s,1)
                    resample(j,:) = resample_s(j,resample_b(j,:));
                end
                cd(origDir)
                cfg.resample = resample;
            case 'yes'
                cfg.ivar = [1 2];
                cfg.wvar = 4;
                cfg.uvar = 3;
                cfg.cvar = [];
        end
      
end

cfg.design = design;
stat = ft_freqstatistics(cfg, data);