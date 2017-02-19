## Synopsis

permANOVA a set of functions used to perform permutation statistics in the cluster-based permutation test framework. 
Currently, independent 2-way design (both balanced and unbalanced designs) and a mixed-design ANOVA

## Setting up Fieldtrip to run permANOVA

To run permutation ANOVAs within the Fieldtrip toolbox for:

1. Put the ft_statfun* .m files into the statfun folder of your current Fieldtrip version

2. The actual ANOVA calculations are done by the core functions anova2_cell_mod for the independent 2-way ANOVA (adapted from Arnaud Delorme’s resampling statistical toolkit, 
https://de.mathworks.com/matlabcentral/fileexchange/27960-resampling-statistical-toolkit), and between_within for the mixed design ANOVA. Put those functions into the privat subfolder of the statfun folder 

3. Next, you’ll have to modify some Fieldtrip functions. Note that the line numbers given below might be slightly different in your Fieldtrip version. 
In ft_statistics_montecarlo:
   * around line 120 (where the defaults for the main function are set), add the following line to set the default effect of interest to be the interaction: 
   
            cfg.fac = ft_getopt(cfg, 'fac','iaxb');

   * starting at line 213, replace the assignment 
   
            resample = resampledesign(cfg, design); 
     
     with

          if isfield(cfg,'resample') 
             resample = cfg.resample; 
          else 
             resample = resampledesign(cfg, design); 
          end

      This will allow you to use some pre-computed, more complex permutation matrices (not necessary for an independent 2-way  ANOVA, but e.g. for group x condition interactions in a mixed design) 

    * at approximately line 229 add 'tmpcfg.fac = cfg.fac;'
      This configuration struct field will be used to indicate the factor or interaction of interest. This can be 'a' (first factor), 'b' (second factor) or 'iaxb' for the interaction

      In resampledesign (in the private folder of your Fieldtrip version) change line 130 from 

            resample = cat(2, blockres{:}); 
         to 

            resample(:,cat(2, blocksel{:})) = cat(2, blockres{:}); 

      See the following bug note: http://bugzilla.fcdonders.nl/show_bug.cgi?id=1546

      Now your Fieldtrip version is set up and ready to run permutation ANOVAs.


## Example: 2-way balanced independent ANOVA

You can take a look at the example script and the Sim_*.m files to ge a better grasp of how to run permutation ANOVAs. Before running the simulation scripts make sure you have a copy of the core functions (anova2_cell_mod and mixed_b1_w1_anova)and the dummy structure in your working directory (this won’t be necessary if you later run your own analysis in Fieldtrip). The set-up of the design matrix and of the independent and control variables cfg.ivar and cfg.cvar follows the same logic as for one-factorial designs in Fieldtrip (check out the [Fieldtrip Tutorial](http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock) and [Fieldtrip FAQs](http://www.fieldtriptoolbox.org/faq/how_can_i_use_the_ivar_uvar_wvar_and_cvar_options_to_precisely_control_the_permutations?s[]=design&s[]=matrix).


The Sim_indepANOVA_* files exemplify how to perform exact permutation tests with restricted permutations and approximate tests on either the raw data or the residuals under a reduced model as described in ‘Permutation tests for multi-factorial analysis of variance. Marti Anderson and Cajo Ter Braak. Journal of Statistical Computation and Simulation Vol. 73, Iss. 2, 2003.  

You'll have to run the permutation ANOVA for each effect separately by specifying cfg.fac = 'a', 'b' or 'iaxb'. The simulation scripts mainly serve the purpose to demonstrate the validity of the different permutation approaches. You might want to reduce the number of repetitions Rep to a lower number (e.g. 10 instead of 1000) in order to avoid overly long running times.

