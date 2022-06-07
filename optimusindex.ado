*! optimusindex 1.0.0 June 2022
*! author Grady Killeen, adapted from code by Michael Anderson and Jeremy Magruder
* Implementation of optimus index calculation from Anderson and Magruder (2022), "Highly Powered Analysis Plans"

program define optimusindex, eclass
  syntax varlist [if] [aweight], treatment(varname) ///
  [cluster(varlist) stratify(varlist) covariates(varlist) unweighted(integer 0) bootstrap_cov(integer 0) cov_bootstrapreps(integer 10) ///
  folds(integer 10) fold_seed(integer) fold_iterations(integer 100) prescreen_cutoff(real 1.2) expshare(real 0.5) cov_shrinkage(real 0.5) ///
  concen_weight(real 0.5) ri_iterations(integer 500) onesided(integer 0)]

  tempvar eligible treatment_bs treatment_r temp_resid all_covariates optindex_all b_null less ///
    b_0 fold rank max all_outcomes ///
    treatment_cells treatment_prob sort_order cell_centile

  forvalues f = 1(1)`folds' {
    tempvar optindex_`f' optindex_pos_`f' optindex_neg_`f'
  }

  local varlist: list uniq varlist // Trim any duplicate variables from the list of outcomes

  quietly gen `sort_order' = _n
  if ( "`covariates'" == "" ) {
    quietly gen byte `all_covariates' = 1
  }
  else {
    quietly egen int `all_covariates' = rowmiss( `covariates' )
    quietly replace `all_covariates' = ( `all_covariates' == 0 )
  }
  quietly egen int `all_outcomes' = rowmiss( `varlist' )
  quietly replace `all_outcomes' = ( `all_outcomes' == 0 )
  if ( "`if'" == "" ) {
    local if "if 1"
  }
  if ( !`hc3' ) {
    local std_errors "cluster `cluster'"
  }
  else {
    local std_errors "hc3"
  }

  if ("`fold_seed'" != "") {
    set seed `fold_seed'
    local fold_iterations = 1  // In this case we only need to take a single split of the data
  }

  mata: actual_p = J(`fold_iterations', 1, .)
  mata: ri_p = J(`fold_iterations', `ri_iterations', .)

  forvalues f_iter = 1(1)`fold_iterations' {

    cluster_random_number `cluster', output_var(`random_clust')

    // Stratify folds on treatment assignment and strata used for randomization GSK: Added stratification on randomization strata. Is this correct?
    tempvar strata_grp random_clust_temp
    quietly egen `strata_grp' = group(`treatment' `stratify')
  	replace `fold' = .
  	local step = 100 / `folds'

    local end_cent = `folds' + 1
    quietly summarize `strata_grp', meanonly
    forval st_nm = 1/`r(max)' {
        quietly centile `random_clust' if `strata_grp' == `st_nm', centile( 0(`step')100 )
        forvalues f = 1(1)`folds' {
      		quietly replace `fold' = `f' if ( `random_clust' >= r(c_`f') & `random_clust' < r(c_`++f') ) & `strata_grp' == `st_nm'
      	}
        quietly replace `fold' = `folds' if ( `random_clust' >= r(c_`end_cent') & `random_clust' !=. ) & `strata_grp' == `st_nm'
      	quietly replace `fold' = 1 if `random_clust' <= r(c_1) & `strata_grp' == `st_nm'
    }

    // Calculate the index on the actual data
    // GSK: Is there any reason we can't do this as opposed to passing a treatment_resid variable?
    quietly reg `treatment_bs' `covariates' [`weight' `exp'] `if' // Residualize the treatment for SUR
    quietly predict `treatment_r' if e(sample), resid

    forvalues f = 1(1)`folds' {
      // Estimate coefficients and covariance matrix for optimus
      gen_coeff_and_covar	`varlist' `if' & `fold_var' != `f', treatment_r(`treatment_r') covariates(`covariates') treatment_cluster(`treatment_cluster') ///
        rdraw(`rdraw') all_covariates(`all_covariates') fv(`f') prescreen_cutoff(`prescreen_cutoff') bootstrap_cov(`bootstrap_cov') ///
        cov_bootstrapreps(`cov_bootstrapreps')
      // Get the optimal index
      capture mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///  GSK: Eliminated ebayes shrinkage of beta
        min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `expshare', `cov_shrinkage', ///
        `unweighted_indices', cov_V, `concen_weight' )

      if ( _rc == 430 ) {
        dis "Mata optimizer failed to find a solution on rep `r'. Trying slightly different covariance shrinkage factor for this rep."
        local temp_cov_shrink = `cov_shrinkage' * 0.9
        capture mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
          min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `expshare', `temp_cov_shrink', ///
          `unweighted_indices', cov_V, `concen_weight' )
        if ( _rc == 430 ) {
          dis "Mata optimizer failed to find a solution on rep `r'. Reverting to unweighted index for this rep."
          mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
            min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `expshare', `cov_shrinkage', 1, ///
            cov_V, `concen_weight' )
        }
      }
      mata: index_selector = index_output[2..rows(index_output),1+negative]
      mata: st_view( optindex_`f' = ., ., "`optindex_`f''" )
      mata: optindex_`f'[.] = Y_`f' * index_selector / colsum( index_selector ) // Generate the optimal index for fold f
      quietly replace `optindex_`f'' = 0 `if' & ( `fold_var' != `f' )

      mata: weights_sorted = sort( ( depvar_`f', strofreal( index_selector ) ), -2 )
      mata: st_numscalar( "count_nontrivial_weight", sum( index_selector :> 1 / `vars_in_group' - 10^-5 ) )
      local index_size_sum = count_nontrivial_weight + `index_size_sum'
      mata: st_numscalar( "count_any_weight", sum( index_selector :> 10^-5 ) )
      local count_any_weight = count_any_weight
      forvalues i = 1(1)`count_any_weight' {
        mata: st_local( "depvar_name", weights_sorted[`i', 1] )
        quietly replace `index_var'_`f'_`i' = "`depvar_name'" in `obs_counter'
        mata: st_numscalar( "weight_val", strtoreal( weights_sorted[`i', 2] ) )
        quietly replace `index_var_wgt'_`f'_`i' = weight_val in `obs_counter'
      }
    }
    if ( `onesided' ) {
      mata: st_numscalar( "totalvotes", sum( votes ) )
      local pct_pos = totalvotes / `folds'
      if ( `pct_pos' > 0.66 ) {
        quietly egen `optindex_all' = rowmean( `optindex_pos_1'-`optindex_pos_`folds'' ) `if' // Generate optimus across all folds
        local vote_outcome = 1
        local sign "pos"
      }
      else if ( `pct_pos' < 0.33 ) {
        quietly egen `optindex_all' = rowmean( `optindex_neg_1'-`optindex_neg_`folds'' ) `if' // Generate optimus across all folds
        local vote_outcome = 2
        local sign "neg"
      }
      else {
        local fold_disagree = 1
        local vote_outcome = 0
        quietly gen `optindex_all' = runiform() `if'
      }
      if ( `fold_disagree' == 0) {  // GSK: Find a better way to store average weights
        forvalues f = 1(1)`folds' {
          forvalues i = 1(1)`count_any_weight_`sign'_`f'' {
            mata: st_local( "depvar_name", weights_sorted_`sign'_`f'[`i', 1] )
            quietly replace `index_var'_`f'_`i' = "`depvar_name'" in `obs_counter'
            mata: st_numscalar( "weight_val", strtoreal( weights_sorted_`sign'_`f'[`i', 2] ) )
            quietly replace `index_var_wgt'_`f'_`i' = weight_val in `obs_counter'
          }
        }
      }
    }
    else {
      quietly egen `optindex_all' = rowmean( `optindex_1'-`optindex_`folds'' ) `if' // Generate optimus across all folds
    }
    quietly replace `optindex_all' = `optindex_all' * `folds' `if'

    // Regress optimus on treatment
    quietly reg `optindex_all' `treatment' `covariates' [`weight' `exp'] `if', vce( `std_errors' )
    local b_0 = _b[`treatment']
    if ( `fold_disagree' ) {  // GSK: What's the rationale for this versus raising an error?
      local b_0 = 0
    }
    mata: actual_p[`fold_iter', 1] = `b_0'

    if ( `onesided' ) {  // GSK: Figure out what this is doing and how to adapt it.
      local avg_index_size = ( 1 - `fold_disagree' ) * ///
      ( `index_size_sum_pos' * ( `vote_outcome' == 1 ) + `index_size_sum_neg' * ( `vote_outcome' == 2 ) ) / `folds'
      quietly replace `fold_disagree_var' = `fold_disagree' in `obs_counter'
      if ( `fold_disagree' == 0 ) {
        quietly replace `positive_var' = ( `vote_outcome' == 1 ) in `obs_counter'
      }
    }
    else {
      local avg_index_size = `index_size_sum' / `folds'
    }
  }

  // End calculation on actual data

  quietly gen byte `treatment_bs' = `treatment' if `all_outcomes'
  if ( "`stratify'" != "" ) {
      quietly egen int `treatment_cells' = group( `stratify' `fold_var') `if' & `treatment_bs' != .
  }
  else {
      quietly gen int `treatment_cells' = `fold_var' `if' & `treatment_bs' != .
  }
  sum `treatment_cells', meanonly
  local total_cells = r(max)
  if ( `total_cells' > 500 ) {
    dis "WARNING: You have specified that treatment probability varies across more than 500 cells!"
  }
  quietly bysort `treatment_cells': egen double `treatment_prob' = mean( `treatment_bs' ) `if' & `treatment_bs' != . & `all_covariates'
  sort `sort_order'

  local index_size_sum = 0
  local index_size_sum_pos = 0
  local index_size_sum_neg = 0
    if ( `onesided' ) {
      forvalues f = 1(1)`folds' {
        quietly gen `optindex_pos_`f'' = .
      }
      forvalues f = 1(1)`folds' {
        quietly gen `optindex_neg_`f'' = .
      }
    }
    else {
      forvalues f = 1(1)`folds' {
        quietly gen `optindex_`f'' = .
      }
    }
    cluster_random_number `cluster', output_var(`rdraw')
    if ( `r' > 0 ) {
      if ( "`fwer_type'" == "rw" ) {
        quietly replace `rdraw' = `rw_rdraw_var'`r'
      }
      quietly bysort `treatment_cells': egen long `rank' = rank( `rdraw' ) `if' & `treatment_bs' != ., track
      quietly by `treatment_cells': egen long `max' = max( `rank' ) `if' & `treatment_bs' != .
      sort `sort_order'
      quietly gen double `cell_centile' = `rank' / `max' `if'
      assert `cell_centile' >= 0 & `cell_centile' <= 1 if `cell_centile' != .
      quietly replace `treatment_bs' = `cell_centile' - 10^-9 < `treatment_prob' if `treatment_bs' != .
      quietly replace `treatment_bs' = 1 if abs( `treatment_prob' - 1 ) < 10^-6 & `treatment_bs' != . // Overwrite 100% treated cells
      quietly replace `treatment_bs' = 0 if abs( `treatment_prob' ) < 10^-6 & `treatment_bs' != . // Overwrite 0% treated cells
      quietly gen byte `eligible' = 0
      quietly replace `eligible' = 1 `if'
      quietly replace `treatment_bs' = . if !`eligible'
      quietly reg `treatment_bs' `covariates' [`weight' `exp'] `if' // Residualize the treatment for SUR
      quietly predict `treatment_r' if e(sample), resid
      drop `eligible' `rank' `max' `cell_centile'
    }
    else {
      quietly gen `treatment_r' = `treatment_resid' // Residualized treatment
    }
    // Estimate and store optimus over different folds
    mata: votes = J( `folds', 1, . ) // Will fill in with 1 for positive and 0 for negative
    local fold_disagree = 0
    forvalues f = 1(1)`folds' {
      // Estimate coefficients and covariance matrix for optimus
      gen_coeff_and_covar	`varlist' `if' & `fold_var' != `f', treatment_r(`treatment_r') covariates(`covariates') treatment_cluster(`treatment_cluster') ///
        rdraw(`rdraw') all_covariates(`all_covariates') fv(`f') prescreen_cutoff(`prescreen_cutoff') bootstrap_cov(`bootstrap_cov') ///
        cov_bootstrapreps(`cov_bootstrapreps') depvar_suffix(`depvar_suffix')
      // Get the optimal index
      capture mata: index_output = summary_power_calc( ebayes( beta', 1, `beta_shrinkage' )', omega, group_selector( select_matrix, ///
        min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `expshare', `cov_shrinkage', ///
        `unweighted_indices', cov_V, `concen_weight' )
      if ( _rc == 430 ) {
        dis "Mata optimizer failed to find a solution on rep `r'. Trying slightly different covariance shrinkage factor for this rep."
        local temp_cov_shrink = `cov_shrinkage' * 0.9
        capture mata: index_output = summary_power_calc( ebayes( beta', 1, `beta_shrinkage' )', omega, group_selector( select_matrix, ///
          min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `expshare', `temp_cov_shrink', ///
          `unweighted_indices', cov_V, `concen_weight' )
        if ( _rc == 430 ) {
          dis "Mata optimizer failed to find a solution on rep `r'. Reverting to unweighted index for this rep."
          mata: index_output = summary_power_calc( ebayes( beta', 1, `beta_shrinkage' )', omega, group_selector( select_matrix, ///
            min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `expshare', `cov_shrinkage', 1, ///
            cov_V, `concen_weight' )
        }
      }
      mata: negative = ( index_output[1,1] < index_output[1,2] )
      if ( `onesided ') {
        mata: votes[`f'] = 1 - negative
        mata: index_selector = index_output[2..rows(index_output),.]
        local sign_count = 1
        foreach sign in pos neg {
          mata: st_view( optindex_`sign'_`f' = ., ., "`optindex_`sign'_`f''" )
          mata: optindex_`sign'_`f'[.] = Y_`f' * index_selector[.,`sign_count'] / colsum( index_selector[.,`sign_count++'] ) // Generate the optimal index for fold f
          quietly replace `optindex_`sign'_`f'' = 0 `if' & ( `fold_var' != `f' )
        }
      }
      else {
        mata: index_selector = index_output[2..rows(index_output),1+negative]
        mata: st_view( optindex_`f' = ., ., "`optindex_`f''" )
        mata: optindex_`f'[.] = Y_`f' * index_selector / colsum( index_selector ) // Generate the optimal index for fold f
        quietly replace `optindex_`f'' = 0 `if' & ( `fold_var' != `f' )
      }
      if ( `r' == 0 ) {
        if ( `onesided ') {
          local sign_count = 1
          foreach sign in pos neg {
            mata: weights_sorted_`sign'_`f' = sort( ( depvar_`f', strofreal( index_selector[.,`sign_count'] ) ), -2 )
            mata: st_numscalar( "count_nontrivial_weight", sum( index_selector[.,`sign_count'] :> 1 / `vars_in_group' - 10^-5 ) )
            local index_size_sum_`sign' = count_nontrivial_weight + `index_size_sum_`sign''
            mata: st_numscalar( "count_any_weight", sum( index_selector[.,`sign_count++'] :> 10^-5 ) )
            local count_any_weight_`sign'_`f' = count_any_weight
          }
        }
        else {
          mata: weights_sorted = sort( ( depvar_`f', strofreal( index_selector ) ), -2 )
          mata: st_numscalar( "count_nontrivial_weight", sum( index_selector :> 1 / `vars_in_group' - 10^-5 ) )
          local index_size_sum = count_nontrivial_weight + `index_size_sum'
          mata: st_numscalar( "count_any_weight", sum( index_selector :> 10^-5 ) )
          local count_any_weight = count_any_weight
          forvalues i = 1(1)`count_any_weight' {
            mata: st_local( "depvar_name", weights_sorted[`i', 1] )
            quietly replace `index_var'_`f'_`i' = "`depvar_name'" in `obs_counter'
            mata: st_numscalar( "weight_val", strtoreal( weights_sorted[`i', 2] ) )
            quietly replace `index_var_wgt'_`f'_`i' = weight_val in `obs_counter'
          }
        }
      }
    }
    if ( `onesided' ) {
      mata: st_numscalar( "totalvotes", sum( votes ) )
      local pct_pos = totalvotes / `folds'
      if ( `pct_pos' > 0.66 ) {  // GSK: Where does pct_pos of .66 come from?
        quietly egen `optindex_all' = rowmean( `optindex_pos_1'-`optindex_pos_`folds'' ) `if' // Generate optimus across all folds
        local vote_outcome = 1
        local sign "pos"
      }
      else if ( `pct_pos' < 0.33 ) {
        quietly egen `optindex_all' = rowmean( `optindex_neg_1'-`optindex_neg_`folds'' ) `if' // Generate optimus across all folds
        local vote_outcome = 2
        local sign "neg"
      }
      else {
        local fold_disagree = 1
        local vote_outcome = 0
        quietly gen `optindex_all' = runiform() `if'
      }
      if ( `fold_disagree' == 0 & `r' == 0 ) {
        forvalues f = 1(1)`folds' {
          forvalues i = 1(1)`count_any_weight_`sign'_`f'' {
            mata: st_local( "depvar_name", weights_sorted_`sign'_`f'[`i', 1] )
            quietly replace `index_var'_`f'_`i' = "`depvar_name'" in `obs_counter'
            mata: st_numscalar( "weight_val", strtoreal( weights_sorted_`sign'_`f'[`i', 2] ) )
            quietly replace `index_var_wgt'_`f'_`i' = weight_val in `obs_counter'
          }
        }
      }
    }
    else {
      quietly egen `optindex_all' = rowmean( `optindex_1'-`optindex_`folds'' ) `if' // Generate optimus across all folds
    }
    quietly replace `optindex_all' = `optindex_all' * `folds' `if'
    // Regress optimus on treatment
    quietly reg `optindex_all' `treatment_bs' `covariates' [`weight' `exp'] `if', vce( `std_errors' )
    local b = _b[`treatment_bs']
    if ( `fold_disagree' ) {
      local b = 0
    }
    local se = _se[`treatment_bs']
    local N = e(N)
    if ( `r' == 0 ) {
      quietly gen `b_0' = `b'
      quietly gen `t_0' = `b' / `se'
      if ( `onesided' ) {
        local avg_index_size = ( 1 - `fold_disagree' ) * ///
        ( `index_size_sum_pos' * ( `vote_outcome' == 1 ) + `index_size_sum_neg' * ( `vote_outcome' == 2 ) ) / `folds'
        quietly replace `fold_disagree_var' = `fold_disagree' in `obs_counter'
        if ( `fold_disagree' == 0 ) {
          quietly replace `positive_var' = ( `vote_outcome' == 1 ) in `obs_counter'
        }
      }
      else {
        local avg_index_size = `index_size_sum' / `folds'
      }
      fillin_storage_vars `outcome_name' `kfold_coefficient' `kfold_se' `index_size'  `max_index_size' `iteration' `kfold_N' ///
        `beta_shrinkage_factor' `cov_shrinkage_factor' `concentration_weight' in `obs_counter', outcome_name(optindex_`outcome_name_string') ///
        var2_val(`b') var3_val(`se') var4_val(`avg_index_size') var5_val(`vars_in_group') var6_val(`iteration_val') var7_val(`N') ///
        var8_val(`beta_shrinkage') var9_val(`cov_shrinkage') var10_val(`concen_weight')
      if ( "`fwer_type'" == "rw" ) {
        mata: b_rw_actual[ `rw_hcount' ] = `b' // Store the coefficient for RW (one-sided test already applied above)
        mata: se_rw_actual[ `rw_hcount' ] = `se' // Store the std error for RW
      }
    }
    else {
      quietly replace `t_null' = `b' / `se' in `r'
      quietly replace `b_null' = `b' in `r'
      if ( "`fwer_type'" == "rw" ) {
        mata: b_rw_null[ `r' , `rw_hcount' ] = `b' // Store the coefficient for RW
        mata: se_rw_null[ `r' , `rw_hcount' ] = `se' // Store the std error for RW
      }
    }
    if ( `onesided' ) {
      drop `treatment_r' `rdraw' `optindex_all' `optindex_pos_1'-`optindex_pos_`folds'' `optindex_neg_1'-`optindex_neg_`folds''
    }
    else {
      drop `treatment_r' `rdraw' `optindex_all' `optindex_1'-`optindex_`folds''
    }
  // End of for r loop
  quietly gen `t_null_jitter' = `t_null' + runiform() * 10^-4 // Add some jitter to guarantee no ties
  quietly egen int `rank_big' = rank( `t_null_jitter' ), field
  quietly egen int `rank_small' = rank( `t_null_jitter' ), track
  foreach size in small big {
    local endcount = wordcount("``size'_null_results'")
    tokenize ``size'_null_results'
    forvalues i = 1(1)`endcount' {
      sum `t_null_jitter' if `rank_`size'' == `i', meanonly
      quietly replace ``i'' = r(mean) in `obs_counter'
    }
  }
  foreach stat in t b {
    assert ``stat'_0' != .
    quietly gen byte `less' = abs( ``stat'_0' ) < abs( ``stat'_null' - 10^-7 ) if ``stat'_null' != .
    sum `less', meanonly // Compute fraction of time the actual average k-fold stat is less than the null stats
    local `stat'_p = r(mean)
    drop `less'
  }
  fillin_storage_vars `outcome_name' `kfold_p' in `obs_counter', outcome_name(optindex_`outcome_name_string') var2_val(`t_p')
  fillin_storage_vars `outcome_name' `kfold_beta_p' in `obs_counter', outcome_name(b_p_`outcome_name_string') var2_val(`b_p')
end



// Define helper functions (non-mata)
// Program to generate random numbers at the treatment cluster level
// Takes as input the cluster id variable(s)
cap program drop cluster_random_number
program define cluster_random_number
	syntax varlist, output_var(string)
	tempvar order select random random2
	quietly gen long `order' = _n
	quietly egen byte `select' = tag( `varlist' )
	quietly gen double `random' = runiform() if `select'
	quietly gen double `random2' = runiform()
	quietly gen double `output_var' = `random'
	quietly bysort `varlist' `random2': replace `output_var' = `output_var'[_N]
	sort `order'
end

// Program to do a wild cluster bootstrap draw:
cap program drop wild_resid
program define wild_resid
	syntax varlist [if], new_resid(string) [cluster(varlist)]
	tempvar sign random_var
	cluster_random_number `cluster', output_var(`random_var')
	quietly gen byte `sign' = 2 * ( `random_var' > 0.5 ) - 1 `if'
	quietly gen `new_resid' = `sign' * `1' `if'
end

// Program to residualize variables
// Takes as input variables to be residualized, with the covariates as an option
// If no covariates are specified, it just demeans
// Default suffix is "_r"
// If a variable with the output varname already exists, it gets dropped
cap program drop residualize_variables
program define residualize_variables
	syntax varlist [if] [aweight], [covariates(varlist) suffix_out(string) suffix_in(string)]
	foreach var of varlist `varlist' {
		if ( "`suffix_in'" == "" ) {
			quietly reg `var' `covariates' [`weight' `exp'] `if'
		}
		else {
			quietly reg `var'`suffix_in' `covariates' [`weight' `exp'] `if'
		}
		if ( "`suffix_out'" == "" ) {
			cap drop `var'_r
			quietly predict `var'_r `if', resid
		}
		else {
			cap drop `var'`suffix_out'
			quietly predict `var'`suffix_out' `if', resid
		}
	}
end

// Program to estimate coefficient vector and covariance matrix to use as inputs for optimus index calculator
// Uses sureg to estimate SUR
// Takes as input a varlist of outcome variables for the subgroup in question
// It expects that standardized/residualized versions of the variables exists with suffixes _sdr and _r
cap program drop gen_coeff_and_covar
program define gen_coeff_and_covar
	syntax varlist [if] [aweight], treatment_r(varname) treatment_cluster(varname) rdraw(varname) ///
		all_covariates(varname) fv(integer) [covariates(varlist) prescreen_cutoff(real 1.2) bootstrap_cov(integer 0) cov_bootstrapreps(real 10) ///
		inflation_factor_base(real 0.3) inflation_factor_rho(real 0.3)] // depvar_suffix option is only for compatability with gen_coeff_and_covar_alt
	local count_full = 0
	local count_catch = 0
	if ( "`if'" == "" ) {
		local if "if 1"
	}
	// Ensure K < N and run screening regression for perfect colinearity in outcomes
	local vars_in_group = wordcount( "`varlist'" )
	mata: temp_string_full = J( 1, `vars_in_group', "" ) // Row vector for all outcome names
	mata: temp_string_catch = J( 1, 1, "" ) // Scalar for first outcome name
	if ( `vars_in_group' < 6 ) {
		local prescreen_cutoff = 0
	}
	else if ( `vars_in_group' < 12 ) { // GSK: What exactly is the prescreen cutoff doing?
		local prescreen_cutoff = min( 0.5, `prescreen_cutoff' )
	}
	foreach outcome of varlist `varlist' {
		quietly reg `outcome'_sd `treatment_r' `covariates' [`weight' `exp'] `if', cluster( `treatment_cluster' ) // Screen out small t-stats and drop dep vars with almost no variation to make the SUR system manageable
		if ( `++count_catch' == 1 ) {
			mata: temp_string_catch[1] = "`outcome'" // Record first outcome name
			local outcome1 `outcome'
		}
		if ( abs( _b[`treatment_r'] / _se[`treatment_r'] ) > `prescreen_cutoff' & e(rmse) > 0.05 ) {
			local prescreen_level "`prescreen_level' `outcome'`depvar_suffix'" // List of outcomes (SDR)
			local screened_group_level "`screened_group_level' `outcome'" // List of outcomes (raw)
			mata: temp_string_full[`++count_full'] = "`outcome'" // Fill in vector of outcome names
		}
	}
	if ( `count_full' == 0 ) {
		dis "No outcome cleared the t-threshold of `prescreen_cutoff' for level index!"
		mata: bname = temp_string_catch[1] // Just take the first outcome if nothing looks good
		local prescreen_level "`prescreen_level' `outcome1'`depvar_suffix'" // First outcome (SDR)
		local screened_group_level "`screened_group_level' `outcome1'" // First outcome (raw)
	}
	else {
		mata: bname = temp_string_full[1..`count_full'] // Final vector of outcome names
	}
	quietly reg `rdraw' `prescreen_level' [`weight' `exp'] `if' & `all_covariates'
	mata: b_screen = st_matrix( "e(b)" )
	mata: st_numscalar( "length_b", length( b_screen ) )
	if ( e(rank) < length_b ) {
		mata: end_obs = length( b_screen ) - 1
		mata: b_selected = select( bname, b_screen[1..end_obs] )
		mata: st_local( "indep_outcomes", invtokens( b_selected ) )
		local group_varlist `indep_outcomes'
	}
	else {
		local group_varlist `screened_group_level'
	}
	local group_count = 0
	foreach outcome of varlist `group_varlist' { // GSK: So SUR is used to estimate the coefficients and covariance matrix (within each of the other folds) for index construction?
		local estimation_cmd "`estimation_cmd' (`outcome'`depvar_suffix' `treatment_r')" // Construct the preliminary sureg command
		mata: y`++group_count' = st_data( ., "`outcome'`depvar_suffix'" ) // Get the LHS vars into Mata
		if ( `group_count' == 1 ) {
			mata: Y_`fv' = y1
			mata: depvar_`fv' = "`outcome'"
		}
		else {
			mata: Y_`fv' = ( Y_`fv', y`group_count' )
			mata: depvar_`fv' = ( depvar_`fv' \ "`outcome'" )
		}
	}
	quietly sureg `estimation_cmd' [`weight' `exp'] `if' // Run the preliminary sureg regression for optimus
	// Select the rows and columns of interest from the coefficient and covariance matrices
	local zeros_needed = round( e(k) / `group_count' ) - 1  // GSK: What is this doing?
	mata: select_colsrows = J( 1, `group_count', ( 1 , J( 1, `zeros_needed', 0 ) ) )
	mata: beta = select( st_matrix( "e(b)" ), select_colsrows )
	mata: omega = select( select( st_matrix( "e(V)" ), select_colsrows ), select_colsrows' )
	local N = e(N)
	mata: N = `N'
	local depvars = e(depvar)
	if ( `bootstrap_cov' ) {
		// Bootstrap the covariance matrix std errors if necessary (GSK: Why would these need to be bootstrapped?)
		preserve
		tempvar sureg_sample
		quietly gen byte `sureg_sample' = e(sample)
		foreach depvar of varlist `depvars' {
			tempvar `depvar'_xb `depvar'_r `depvar'_wr
			quietly predict ``depvar'_xb' if `sureg_sample', eq(`depvar') xb
			quietly predict ``depvar'_r' if `sureg_sample', eq(`depvar') resid
			quietly drop `depvar'
		}
		forvalues b = 1(1)`cov_bootstrapreps' {
			foreach depvar of newlist `depvars' {
				wild_resid ``depvar'_r' if `sureg_sample', new_resid(``depvar'_wr') cluster(`treatment_cluster')
				quietly gen `depvar' = ``depvar'_xb' + ``depvar'_wr'
				drop ``depvar'_wr'
			}
			quietly sureg `estimation_cmd' [`weight' `exp']
			assert `zeros_needed' == round( e(k) / `group_count' ) - 1
			mata: select_colsrows = J( 1, `group_count', ( 1 , J( 1, `zeros_needed', 0 ) ) )
			mata: omega_b`b' = select( select( st_matrix( "e(V)" ), select_colsrows ), select_colsrows' )
			if ( `b' == 1 ) {
				mata: bs_V = rowshape( omega_b`b' , 1 )
			}
			else {
				mata: bs_V = ( bs_V \ rowshape( omega_b`b', 1 ) )
			}
			local depvars_new = e(depvar)
			drop `depvars_new'
		}
		mata: cov_V = rowshape( diagonal( variance( bs_V ) )', rows( omega_b1 ) )
		restore
	}
	else {
		quietly corr `depvars' `if'
		mata: rho_mat = st_matrix("r(C)")
		mata: _diag( rho_mat, J( rows( rho_mat ), 1, . ) )
		mata: cov_V = ( 1 :+ `inflation_factor_base' :+ ( `inflation_factor_rho' :* abs( rho_mat ) ) ) :* ( diagonal(omega) * diagonal(omega)' / N )
	}
end

// Define mata programs

mata

// Function for sorting a vector (in descending order) and re-sorting its covariance matrix so that it corresponds to the sorted vector
// Takes a vector and a covariance matrix that has been normalized to have unit variances and correlation coefficients
// Returns the re-sorted covariance matrix (the original vector remains unsorted)
function sort_cov_matrix(real matrix V, real vector beta, real scalar descending)
{
	assert( ( abs( diagonal( V ) :- 1 ) :< 10^-4 ) == J( rows( V ), 1, 1 ) ) // Check diagonals are 1
	assert( ( abs( V ) :< 1.0001 ) == J( rows( V ), cols(V), 1 ) ) // Check off-diagonals are between -1 and 1
	assert( descending == 0 | descending == 1)
	b = beta
	if ( rows( beta ) == 1 ) {
		b = b' // Transpose to column vector if necessary
	}
	b_indexed = ( range( 1, length( b ), 1 ), b )
	b_indexed = sort( b_indexed, 2 - 4 * descending )
	endval = rows( V )
	V_sort = I( endval )
	for (i = 1; i <= endval; i++ ) {
		for( j = 1; j <= i; j++ ) {
			V_sort[i,j] = V[b_indexed[i,1],b_indexed[j,1]]
		}
	}
	return( makesymmetric( V_sort ) )
}

// Function for computing optimal summary index, unweighted or weighted
// Summary_power_calc calculates the highest power summary index for a given a set of betas and covariance matrix
// Takes a length-H row vector and returns a (H + 1) x 2 matrix
// For unweighted index, the row vector should be sorted by t-stat. For weighted index, order shouldn't matter.
// Fast flag assumes the optimal index respects monotonicity of t-stats (in either direction)
// Function searches for the optimal summary index in the both directions (positive and negative)
function summary_power_calc( real rowvector beta, real matrix covariance, real matrix select_variables, real scalar crit_val, real scalar fast, real scalar exp_share, real scalar cov_shrinkage_factor, real scalar unweighted, real matrix cov_covariance, real scalar concen_weight )
{
	//	Standardize everything
	V = covariance
	cov_V = cov_covariance
	beta = beta :/ sqrt( diagonal( V ) )'
	endval = length( beta )
	for (i = 1; i <= endval; i++ ) {
		sigma = sqrt( V[i,i] )
		V[i,.] = V[i,.] / sigma
		cov_V[i,.] = cov_V[i,.] / sigma ^ 2
		V[.,i] = V[.,i] / sigma
		cov_V[.,i] = cov_V[.,i] / sigma ^ 2
	}
	assert( ( abs( diagonal( V ) :- 1 ) :< 10^-4 ) == J( rows( V ), 1, 1 ) ) // Check diagonals of correlation matrix are 1
	assert( ( abs( V ) :< 1.0001 ) == J( rows( V ), cols( V ), 1 ) ) // Check off-diagonals of correlation matrix  are between -1 and 1
	// Shrink the covariance terms
	if ( rows( V ) >= 3 & cov_shrinkage_factor > 0 ) {
		V0 = V
		_diag( V0, J( rows( V0 ), 1, 2 ) ) // Replace diagonal elements with 2 so that we can sort them out below
		V_vec = colshape( V0, 1 ) // Convert matrix to column vector
		V_vec_indexed =  ( range( 1, length( V_vec ), 1 ), V_vec ) // Index the vector
		V_vec_indexed = sort( V_vec_indexed, 2 ) // Sort to place diagonal elements last
		trimlength = rows( V0 ) * cols( V0 ) - rows( V0 ) // Compute number of non-diagonal elements
		V_vec_indexed_trim = V_vec_indexed[1..trimlength,.] // Trim off the diagonal elements
		V_vec_indexed_trim = sort( V_vec_indexed_trim, 1 ) // Get vector in correct order
		cov_V0 = cov_V
		_diag( cov_V0, J( rows( cov_V0 ), 1, . ) ) // Replace diagonal elements with missing so that we can sort them out below
		cov_V_vec = colshape( cov_V0, 1 ) // Convert matrix to column vector
		cov_V_vec_indexed =  ( range( 1, length( cov_V_vec ), 1 ), cov_V_vec ) // Index the vector
		cov_V_vec_indexed = sort( cov_V_vec_indexed, 2 ) // Sort to place diagonal elements last
		cov_V_vec_indexed_trim = cov_V_vec_indexed[1..trimlength,.] // Trim off the diagonal elements
		cov_V_vec_indexed_trim = sort( cov_V_vec_indexed_trim, 1 ) // Get vector in correct order
		assert( sort( cov_V_vec_indexed_trim, 1 )[.,1] == sort( V_vec_indexed_trim, 1 )[.,1] ) // Check that sorting worked as expected
		V_vec_shrunk = ebayes_cov( V_vec_indexed_trim[.,2], ceil( cov_shrinkage_factor ), cov_shrinkage_factor, cov_V_vec_indexed_trim[.,2] )
		V_vec_indexed_shrunk = ( ( V_vec_indexed_trim[., 1] \ V_vec_indexed[trimlength+1..rows(V_vec_indexed), 1]) , ( V_vec_shrunk \ J( rows( V0 ), 1, 1 ) ) ) // Add back in the diagonal elements
		V_vec_shrunk = sort( V_vec_indexed_shrunk, 1 )[.,2] // Get vector in correct order
		V = colshape( V_vec_shrunk , rows( V0 ) ) // Convert back to matrix
	}
	if ( ! issymmetric( V ) ) {
		if ( sum( abs( reldif( V, makesymmetric( V ) ) ) ) > 0.001 ) {
			reldif( V, makesymmetric( V ) )
			"Difference between V and a forced symmetric version of V"
			"You have a covariance matrix that is non-symmetric by more than rounding error (see above). This is not allowed."
			exit(198)
		}
		else {
			V = makesymmetric( V )
		}
	}
	if ( rows( V ) == 1 ) {
		power_pos = ( 1 - normal( crit_val - beta / sqrt( V ) ) )
		power_neg = normal( -crit_val - beta / sqrt( V ) )
		return( ( power_pos, power_neg ) \ ( 1, 1 ) )
	}
	else if ( unweighted ) {
		beta_indexed = ( range( 1, length( beta ), 1 ), beta' )
		if ( fast == 1 ) {
			V = sort_cov_matrix( V, beta_indexed[.,2], 1 )
			beta_indexed = sort( beta_indexed, -2 )
			H = length( beta )
			select_variables_internal = ( lowertriangle( J( H, H, 1 ), . ) \ uppertriangle( J( H, H, 1 ), . )[2..H,.] )
		}
		else {
			select_variables_internal = select_variables
		}
		row_size = rowsum( select_variables_internal )
		beta_avg = rowsum( beta_indexed[.,2]' :* select_variables_internal ) :/ ( row_size + ( row_size :== 0 ) * 0.01 ) // If fast == 0, select_variables is 2 ^ group_size x group_size
		end_val = rows( select_variables_internal )
		se_index = J( end_val, 1, . )
		for ( i = 1; i <= end_val; i++ ) {
			se_index[i] = sqrt( ( exp_share / ( 1 - exp_share ) ) * ( select_variables_internal[i, .] * V * select_variables_internal[i, .]' ) ///
				/ ( row_size[i] ^ 2 + ( row_size[i] == 0 ) * 0.01 ) + ( row_size[i] == 0 ) * 0.01 )
		}
		weights_internal = select_variables_internal :/ ( row_size + ( row_size :== 0 ) * 0.01 )
		power_pos = ( beta_avg :/ se_index ) - concen_weight * rowsum( weights_internal :^2 )
		power_neg = -( beta_avg :/ se_index ) - concen_weight * rowsum( weights_internal :^2 )
		max_power_index_pos = .
		max_power_index_neg = .
		ties = .
		maxindex( power_pos, 1, max_power_index_pos, ties )
		max_power_index_pos = max_power_index_pos[1]
		maxindex( power_neg, 1, max_power_index_neg, ties )
		max_power_index_neg = max_power_index_neg[1]
		return_indexed = ( beta_indexed[., 1], select_variables_internal[max_power_index_pos, .]', select_variables_internal[max_power_index_neg, .]' )
		return( ( power_pos[max_power_index_pos], power_neg[max_power_index_neg] ) \ sort( return_indexed, 1 )[., 2..3] )
	}
	else {
		temp_power = J( 2, 1, . )
		w = J( 2, length( beta ), 0 )
		max_power_index = .
		ties = .
		for (negative = 0; negative <= 1; negative++ ) {
			if ( !negative * !sum( beta :>= 0 ) ) {
				temp_power[1] = 0
				maxindex( beta, 1, max_power_index, ties )
				max_power_index = max_power_index[1]
				w[1,max_power_index] = 1
			}
			else if ( !negative * ( sum( beta :>= 0 ) == 1 ) ) {
				maxindex( beta, 1, max_power_index, ties )
				max_power_index = max_power_index[1]
				w[1,max_power_index] = 1
				temp_power[1] = beta[max_power_index] / sqrt( ( exp_share / ( 1 - exp_share ) ) * V[max_power_index,max_power_index] )
			}
			else if ( negative * !sum( beta :< 0 ) ) {
				temp_power[2] = 0
				minindex( beta, 1, max_power_index, ties )
				w[2,max_power_index] = 1
			}
			else if ( negative * ( sum( beta :< 0 ) == 1 ) ) {
				minindex( beta, 1, max_power_index, ties )
				w[2,max_power_index] = 1
				temp_power[2] = abs( beta[max_power_index] / sqrt( ( exp_share / ( 1 - exp_share ) ) * V[max_power_index,max_power_index] ) )
			}
			else {
				S = optimize_init()
				optimize_init_evaluator( S, &optimus() )
				optimize_init_evaluatortype( S, "d1" )
				b = beta
				optimize_init_argument( S, 1, b )
				optimize_init_argument( S, 2, V )
				optimize_init_argument( S, 3, crit_val )
				optimize_init_argument( S, 4, exp_share )
				optimize_init_argument( S, 5, negative )
				optimize_init_argument( S, 6, concen_weight )
				optimize_init_constraints( S, ( J( 1, length( b ), 1 ), 1 ) ) // Constrain weights to sum to 1
				optimize_init_conv_maxiter( S, 100 )
				optimize_init_params( S, J( 1, length( b ), 1 / length( b ) ) + 0.1 * ( b :== max( b ) ) )
				w[1+negative,.] = optimize(S)
				temp_power[1+negative,1] = abs( optimize_result_value(S) )
			}
		}
		return( ( temp_power, w )' )
	}
}

// Optimus weighted power function, d1 solution (with gradient)
// Code the [0,1] interval constraint a little loosely using a high-order polynomial
// Maximize the expected t-stat, since power is a monotonic transformation of that.
void optimus(todo, w, b, V, crit_val, exp_share, negative, con_wgt, power, g, H)
{
	if ( !negative ) {  // GSK: Can we just take the absolute value of the first term to get rid of the negative accounting?
		power = ( b * w' ) / sqrt( ( exp_share / ( 1 - exp_share ) ) * w * V * w' )  -
			.001 * sum( ( 2 * w :- 1 ) :^ 200 ) - con_wgt * sum( w:^2 )  // GSK: con_wgt is the hyper-parameter on the HHI penalty
			// GSK: The .001 * sum( ( 2 * w :- 1 ) :^ 200 ) term appears to (loosely) force values to be between 0 and 1.
		if ( todo >= 1 ) {
			end_val = length( b )
			Vw = V * w'
			for ( i = 1; i <= end_val; i++ ) {
				g[i] =  ( b[i] / sqrt( w * V * w' ) - b * w' * Vw[i] * ( w * V * w' ) ^ -1.5 ) / sqrt( exp_share / ( 1 - exp_share ) ) -
					0.4 * ( 2 * w[i] - 1 ) ^ 199 - con_wgt * 2 * w[i]
			}
		}
	}
	else {
		power = -( b * w' ) / sqrt( ( exp_share / ( 1 - exp_share ) ) * w * V * w' )  -
			.001 * sum( ( 2 * w :- 1 ) :^ 200 ) - con_wgt * sum( w:^2 )
		if ( todo >= 1 ) {
			end_val = length( b )
			Vw = V * w'
			for ( i = 1; i <= end_val; i++ ) {
				g[i] =  -( b[i] / sqrt( w * V * w' ) - b * w' * Vw[i] * ( w * V * w' ) ^ -1.5 ) / sqrt( exp_share / ( 1 - exp_share ) ) -
					0.4 * ( 2 * w[i] - 1 ) ^ 199 - con_wgt * 2 * w[i]
			}
		}
	}
}

// Function to generate a group selection matrix of arbitrary group size (less than max group size)
function group_selector( real matrix max_selector, real scalar target_size )
{
	max_size = cols( max_selector )
	start_col = max_size - target_size + 1
	r = rows( max_selector ) / 2 ^ ( max_size - target_size )
	selector = max_selector[|1,start_col \ r,max_size|]
	return( selector )
}

// Function to generate Empirical Bayes shrunken correlation coefficients
// There is no conversion to absolute values
// Formula: P(rho | rho-hat) = P(rho-hat | rho) * P(rho) / Sum( P(rho-hat | rho) * P(rho) )
// Note that P(rho) is constant here at 1 / (# stats) since we're using the discrete empirical distribution
// Uses normal approximation for computational efficiency
// Empirical Bayes sometimes seems to shrink too much; shrinkage_factor is an arbitrary factor to reduce the shrinkage
// Takes and returns column vectors
function ebayes_cov( real vector corr, real scalar on, real scalar shrinkage_factor, real vector corr_V )
{
	if ( on == 1 ) {
		s_mat = sqrt( J( 1, rows( corr_V ), corr_V ) )
		post_pre = normalden( ( corr :- J( rows( corr ), 1, corr' ) ) :/ s_mat ) // Entry [i, j] gives P(c-hat[i] | c[j]) * P(c[j]) if c[j] equals the j'th t-stat (keeping in mind P(c) is just a constant). We "integrate" over b by taking the average across the empirical distribution.
		post = post_pre :/ rowsum( post_pre ) // Entry [i, j] gives P(c-hat[i] | c[j]) * P(c[j]) / Sum( P(c-hat | c) * P(c) ) if c[j] equals the j'th t-stat (keeping in mind P(c) is just a constant)
		// post is a group_size x group_size matrix of posterior probabilities, P(c | c-hat)
		// Within each row, each entry gives P(c | c-hat) for a specific value of c from the empirical distribution
		return( ( post * corr ) * shrinkage_factor + corr * ( 1 - shrinkage_factor ) )
	}
	else {
		return( corr )
	}
}

end
