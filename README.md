## Synopsis

permANOVA provides a set of functions to perform permutation statistics on factorial designs within the cluster-based permutation test framework of the Fieldtrip toolbox (http://www.fieldtriptoolbox.org/). 
permANOVA implements two-way independent, mixed-design and repeated meausures permutation ANOVAs.

## Set up Fieldtrip

To run permutation ANOVAs within the Fieldtrip toolbox:

1. Put the ft_statfun* .m files in *statfuns_for_ft* into the *statfun* folder of your current Fieldtrip version

2. The actual ANOVA calculations are done by the functions in *statfuns_for_ft/private*. Put those into the private subfolder of the statfun folder 

3. Next, you’ll have to modify some Fieldtrip functions. Note that the line numbers given below might be slightly different in your Fieldtrip version (I used fieldtrip-20200821). 

In ft_statistics_montecarlo:
    * around line 120 (where the defaults for the main function are set), add the following line to set the default effect of interest to be the interaction: 
   
            cfg.fac = ft_getopt(cfg, 'fac','iaxb');

   * starting at line 214, replace the assignment 
   
            resample = resampledesign(cfg, design); 
     
     with

          if isfield(cfg,'resample') 
             resample = cfg.resample; 
          else 
             resample = resampledesign(cfg, design); 
          end

      This will allow you to use some pre-computed, more complex permutation matrices (not necessary for an independent two-way ANOVA, but e.g. for group x condition interactions in a mixed-design ANOVA) 

    * at approximately line 240 add `tmpcfg.fac = cfg.fac;` 
      This configuration struct field will be used to indicate the factor or interaction of interest. This can be 'a' (first factor), 'b' (second factor) or 'iaxb' for the interaction

    In resampledesign (in the private folder of your Fieldtrip version):
    * change line 130 from 

          resample = cat(2, blockres{:}); 

      to

          resample(:,cat(2, blocksel{:})) = cat(2, blockres{:}); 
          
          
    * At about line 160 replace 
    
          error('the design matrix variables should be constant within a block');
          
      with
      
          warning('Fieldtrip:checkinput','the design matrix variables should be constant within a block');
          [lastmsg, wvar_warnid] = lastwarn;
          warning('off', wvar_warnid)
  

  Now your Fieldtrip version is set up and ready to run permutation ANOVAs. Make sure it is in your path and call *ft_defaults*.


## Run simulations

To replicate the simulations first set up Fieldtrip as described above, then:

1. Make a new subfolder *stat_util* in your working directory *<PATH_TO_WORKING_DIR>/stat_util* and copy the files in *statfuns_for_ft/private* into this new folder. This is to allow direct access to the core stats functions.

2. Put the file *dummy.mat* (which contains a data struture in Fieldtrip format) into your working directory.

3. Get the files in *sim_scripts*. Scripts to run the simulations start with the prefix *Sim_*, results can be plotted using the scripts that start with *plotPower_*. *FtSimLink* functions link the simulation scripts to Fieldtrip where the actual computations are performed. Make sure *sim_scripts* is in your path.

## Example: 2-way balanced independent ANOVA

Take a look at the example script and/or the Sim_*  files to ge a better grasp of how to run permutation ANOVAs. The set-up of the design matrix and of the independent and control variables cfg.ivar and cfg.cvar follows the same logic as for one-factorial designs in Fieldtrip (check out the [Fieldtrip Tutorial](http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock) and [Fieldtrip FAQs](http://www.fieldtriptoolbox.org/faq/how_can_i_use_the_ivar_uvar_wvar_and_cvar_options_to_precisely_control_the_permutations?s[]=design&s[]=matrix).

The Sim_indepANOVA_* files demonstrate how to perform exact permutation tests with restricted permutations and approximate tests on either the raw data or the residuals under a reduced model as described in ‘Permutation tests for multi-factorial analysis of variance. Marti Anderson and Cajo Ter Braak. Journal of Statistical Computation and Simulation Vol. 73, Iss. 2, 2003.  

You need to run the permutation ANOVA for each effect separately by specifying cfg.fac = 'a', 'b' or 'iaxb'. The simulation scripts mainly serve the purpose to demonstrate the validity of the different permutation approaches. You may want to reduce the number of repetitions *Rep* to a lower number (e.g. 10 instead of 1000) to avoid very long run times.

