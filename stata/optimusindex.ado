*! optimusindex 1.0.0 June 2022
*! author Grady Killeen, adapted from code by Michael Anderson and Jeremy Magruder
* Implementation of optimus index calculation from Anderson and Magruder (2022), "Highly Powered Analysis Plans"

program define optimusindex, eclass
  syntax varlist [if] [aweight], treatment(varname) ///
  [cluster(varlist) stratify(varlist) covariates(varlist) unweighted(integer 0) bootstrap_cov(integer 0) cov_bootstrapreps(integer 10) ///
  folds(integer 5) fold_seed(integer 0) fold_iterations(integer 100) prescreen_cutoff(real 1.2) cov_shrinkage(real 0.5) ///
  concen_weight(real 0.5) ri_iterations(integer 500) onesided(integer 1) rw(varlist) strat_nulltreat_fold(integer 1)]

  tempfile _torestore
  quietly save `_torestore'

  // Generate group selection matrices for optimal summary index
  clear
  quietly set obs 2

  // Generate binary indicators for all potential index components up to group size 10
  local H = 10
  forvalues v = 1(1)`H' {
  	quietly gen v`v' = 0
  	quietly replace v`v' = 1 in 2
  	local sortlist "`sortlist' v`v'"
  }
  // Generate all combinations of group selection indicators
  fillin v1-v`H'
  // Sort so that we can scale down to smaller max group sizes
  sort `sortlist'
  // Load generated combinations into Mata vectors
  quietly putmata v1-v`H', replace
  //restore
  // use `sbp', clear
  // Merge Mata vectors into single matrix
  mata: select_matrix = v1
  forvalues v = 2(1)`H' {
  	mata: select_matrix = ( select_matrix, v`v')
  }
  use `_torestore', clear

  tempvar eligible treatment_bs treatment_r temp_resid all_covariates optindex_all b_null less ///
    b_0 fold rank max all_outcomes rdraw random_clust strata_grp fold ///
    treatment_cells treatment_prob sort_order cell_centile treatment_r_permuted

  if ("`cluster'" != "") {
    local treatment_cluster `cluster'
  }
  else {
    local treatment_cluster `treatment'
  }

  local expected_crit_val = 2  // Used in power calculations, but not important since we're just maximizing power which is a monotonic function of t-stats

  quietly generate `fold' = .

  forvalues f = 1(1)`folds' {
    tempvar optindex_`f' optindex_pos_`f' optindex_neg_`f'
  }

  local varlist: list uniq varlist // Trim any duplicate variables from the list of outcomes
  local vars_in_group = wordcount( "`varlist'" )

  if ("`rw'" != "") {
    local rw: list uniq rw
    local vars_rw = wordcount("`rw'") + 1  // Plus 1 since the optimus index is also included
  }

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

  quietly replace `treatment' = . if `all_outcomes' == 0
  quietly replace `treatment' = . if `all_covariates' == 0  // Limits the sample to observations with no outcomes or covariates missing

  if ( "`if'" == "" ) {
    local if "if 1"
  }
  if ( "`cluster'" != "" ) {
    local std_errors "cluster `cluster'"
  }
  else {
    local std_errors "hc3"
  }

  // Residualize and standardize outcome variables (only standardize if doesn't exist, assumes treat = 0 gives control)
  foreach var of varlist `varlist' {
    capture confirm variable `var'_sd
    if _rc {
      quietly summarize `var' if `treatment' == 0
      local r_mean = r(mean)
  		local r_sd = r(sd)
  		quietly gen `var'_sd = ( `var' -`r_mean' ) / `r_sd'
    }
  }
  // Residualize, just demeans if no covariates
  if ("`covariates'" != "")  {
    residualize_variables `varlist' `if' [`weight' `exp'], covariates(`covariates')
    foreach var of varlist `varlist' {
      residualize_variables `var'_sd `if' [`weight' `exp'], covariates(`covariates') suffix_out(r)
    }
  }
  else {
    residualize_variables `varlist' `if' [`weight' `exp'], covariates(`covariates')
    foreach var of varlist `varlist' {
      residualize_variables `var'_sd `if' [`weight' `exp'], covariates(`covariates') suffix_out(r)
    }
  }

  if (`fold_seed' != 0) {
    set seed `fold_seed'
    local fold_iterations = 1  // In this case we only need to take a single split of the data
  }

  mata: b_by_iter = J(`fold_iterations', 1, .)  // Store coefficients for each iteration of folds
  if (`onesided') {
        mata: pos_by_iter = J(`fold_iterations', 1, 0)  // Indicate if positive for each iteration
  }
  mata: p_by_iter = J(`fold_iterations', 2, .) // Store unadjusted p-values for each iteration of folds (record row number in second column to preserve index)
  mata: average_weights = J(`vars_in_group', `fold_iterations', 0)  // Will add in average weights across folds
  mata: any_weight = J(1, `fold_iterations', 0)  // Store average number (across folds) of non-zero weights for fold iteration
  mata: nontrivial_weight = J(1, `fold_iterations', 0)  // Store average number of weights > 1/H for fold iteration

  if ("`rw'" != "") {
    mata: p_rw_by_iter = J(`vars_rw', `fold_iterations', .)  // Stores the RW corrected p-value for each variable and iteration if enabled
  }

  forvalues f_iter = 1(1)`fold_iterations' {

    if ( `onesided' ) {
			forvalues f = 1(1)`folds' {
        capture drop `optindex_pos_`f''
				quietly gen `optindex_pos_`f'' = .
			}
			forvalues f = 1(1)`folds' {
        capture drop `optindex_neg_`f''
				quietly gen `optindex_neg_`f'' = .
			}
		}
		else {
			forvalues f = 1(1)`folds' {
        capture drop `optindex_`f''
				quietly gen `optindex_`f'' = .
			}
		}

    // Note: We store p-values so that RI uses a pivotal statistic
    if ("`rw'" != "") {
      mata: p_actual = J(1, `vars_rw', .)  // Store actual (naive) p-values for each variable for RW
      mata: p_student = J(`ri_iterations', `vars_rw', .)  // This stores the permuted p-values for each variable and each RI run
    }
    else {
      mata: p_actual = J(1, 1, .)
      mata: p_student = J(`ri_iterations', 1, .)
    }

    if (`onesided') {
      mata: average_weights_sign = J(`vars_in_group', 2, 0)  // Store average weights for positive/negative folds
    }

    mata: p_by_iter[`f_iter', 2] = `f_iter'  // Record row number
    capture drop `random_clust'
    capture drop `strata_grp'
    cluster_random_number `treatment_cluster', output_var(`random_clust')

    // Stratify folds on treatment assignment and strata used for randomization
    if ("`stratify'" != "") {
      quietly egen `strata_grp' = group(`treatment' `stratify')
    }
    else {
      quietly egen `strata_grp' = group(`treatment')
    }

  	quietly replace `fold' = .
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
    capture drop `treatment_r'
    quietly reg `treatment' `covariates' [`weight' `exp'] `if' // Residualize the treatment for SUR
    capture drop `treatment_r'
    quietly predict `treatment_r' if e(sample), resid

    capture drop `rdraw'
    cluster_random_number `treatment_cluster', output_var(`rdraw')

    mata: votes = J( `folds', 1, . ) // Will fill in with 1 for positive and 0 for negative
		local fold_disagree = 0
    local index_size_sum = 0
  	local index_size_sum_pos = 0
  	local index_size_sum_neg = 0
    local any_weight_sum_pos = 0
    local any_weight_sum_neg = 0
    local any_weight_sum = 0
    forvalues f = 1(1)`folds' {
      // Estimate coefficients and covariance matrix for optimus
      gen_coeff_and_covar	`varlist' `if' & `fold' != `f', treatment_r(`treatment_r') covariates(`covariates') treatment_cluster(`treatment_cluster') ///
        rdraw(`rdraw') all_covariates(`all_covariates') fv(`f') prescreen_cutoff(`prescreen_cutoff') bootstrap_cov(`bootstrap_cov') ///
        cov_bootstrapreps(`cov_bootstrapreps') vce(`std_errors')

      // Get the optimal index
      capture mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
        min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `cov_shrinkage', ///
        `unweighted', cov_V, `concen_weight' )

      if ( _rc == 430 ) {
        dis "Mata optimizer failed to find a solution on rep `r'. Trying slightly different covariance shrinkage factor for this rep."
        local temp_cov_shrink = `cov_shrinkage' * 0.9
        capture mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
          min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `temp_cov_shrink', ///
          `unweighted', cov_V, `concen_weight' )
        if ( _rc == 430 ) {
          dis "Mata optimizer failed to find a solution on rep `r'. Reverting to unweighted index for this rep."
          mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
            min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `cov_shrinkage', 1, ///
            cov_V, `concen_weight' )
        }
      }

      mata: negative = ( index_output[1,1] < index_output[1,2] )  // Indicates that negative index yields better power
      if ( `onesided ') {
				mata: votes[`f'] = 1 - negative
				mata: index_selector = index_output[2..rows(index_output),.]
				local sign_count = 1
				foreach sign in pos neg {
					mata: st_view( optindex_`sign'_`f' = ., ., "`optindex_`sign'_`f''" )
					mata: optindex_`sign'_`f'[.] = Y_`f' * index_selector[.,`sign_count'] / colsum( index_selector[.,`sign_count++'] ) // Generate the optimal index for fold f
					quietly replace `optindex_`sign'_`f'' = 0 `if' & ( `fold' != `f' )
				}
			}
			else {
				mata: index_selector = index_output[2..rows(index_output),1+negative]
				mata: st_view( optindex_`f' = ., ., "`optindex_`f''" )
				mata: optindex_`f'[.] = Y_`f' * index_selector / colsum( index_selector ) // Generate the optimal index for fold f
				quietly replace `optindex_`f'' = 0 `if' & ( `fold' != `f' )
			}

      // Store information about weights (number, number non-trivial, and values)
      if ( `onesided ') {
				local sign_count = 1
				foreach sign in pos neg {
          // Fill weight with 0 if failed prescreen cutoff (since then it is excluded)
          mata: normalized_weights = index_selector[.,`sign_count'] / colsum( index_selector[.,`sign_count'] )
          local j = 1
          forval i = 1(1)`vars_in_group' {
            mata: st_numscalar("var_included", included_vars[`i'])
            local var_included = var_included
            if (`var_included') {
              mata: average_weights_sign[`i',`sign_count'] = average_weights_sign[`i',`sign_count'] + (1/`folds')*normalized_weights[`j']
              local j = `j' + 1
            }
          }
					mata: weights_sorted_`sign'_`f' = sort( ( depvar_`f', strofreal( index_selector[.,`sign_count'] ) ), -2 )
					mata: st_numscalar( "count_nontrivial_weight", sum( index_selector[.,`sign_count'] :> 1 / `vars_in_group' - 10^-5 ) )
					local index_size_sum_`sign' = count_nontrivial_weight + `index_size_sum_`sign''
					mata: st_numscalar( "count_any_weight", sum( index_selector[.,`sign_count++'] :> 10^-5 ) )
					local any_weight_sum_`sign' = count_any_weight + `any_weight_sum_`sign''
				}
			}
			else {
				mata: weights_sorted = sort( ( depvar_`f', strofreal( index_selector ) ), -2 )
				mata: st_numscalar( "count_nontrivial_weight", sum( index_selector :> 1 / `vars_in_group' - 10^-5 ) )
				local index_size_sum = count_nontrivial_weight + `index_size_sum'
				mata: st_numscalar( "count_any_weight", sum( index_selector :> 10^-5 ) )
				local any_weight_sum = count_any_weight + `any_weight_sum'
        // Add in the average weights for this fold
        // Fill weight with 0 if failed prescreen cutoff (since then it is excluded)
        mata: normalized_weights = index_selector / colsum( index_selector )
        local j = 1
        forval i = 1(1)`vars_in_group' {
          mata: st_numscalar("var_included", included_vars[`i'])
          local var_included = var_included
          if (`var_included') {
            mata: average_weights[`i',`f_iter'] = average_weights[`i',`f_iter'] + (1/`folds')*normalized_weights[`j']
            local j = `j' + 1
          }
        }

			}

    }

    capture drop `optindex_all'

    if ( `onesided' ) {
			mata: st_numscalar( "totalvotes", sum( votes ) )
			local pct_pos = totalvotes / `folds'
			if ( `pct_pos' > 0.66 ) {
				quietly egen `optindex_all' = rowmean( `optindex_pos_1'-`optindex_pos_`folds'' ) `if' // Generate optimus across all folds
				local vote_outcome = 1
				local sign "pos"
        mata: average_weights[,`f_iter'] = average_weights_sign[,1]
        mata: pos_by_iter[`f_iter', 1] = 1    // Calculate the index on the actual data
        mata: nontrivial_weight[1, `f_iter'] = (1/`folds')*`index_size_sum_pos'
        mata: any_weight[1, `f_iter'] = (1/`folds')*`any_weight_sum_pos'
			}
			else if ( `pct_pos' < 0.33 ) {
				quietly egen `optindex_all' = rowmean( `optindex_neg_1'-`optindex_neg_`folds'' ) `if' // Generate optimus across all folds
				local vote_outcome = 2
				local sign "neg"
        mata: average_weights[,`f_iter'] = average_weights_sign[,2]
        mata: nontrivial_weight[1, `f_iter'] = (1/`folds')*`index_size_sum_neg'
        mata: any_weight[1, `f_iter'] = (1/`folds')*`any_weight_sum_neg'
			}
			else {
        dis in red "WARNING: Folds do not agree about index sign. Consider re-running with a two-sided test."
        local fold_disagree = 1
        local vote_outcome = 0
        quietly gen `optindex_all' = runiform() `if'
			}
		}
		else {
      mata: nontrivial_weight[1, `f_iter'] = (1/`folds')*`index_size_sum'
      mata: any_weight[1, `f_iter'] = (1/`folds')*`any_weight_sum'
			quietly egen `optindex_all' = rowmean( `optindex_1'-`optindex_`folds'' ) `if' // Generate optimus across all folds
		}
		quietly replace `optindex_all' = `optindex_all' * `folds' `if'

    // Regress optimus on treatment
    quietly reg `optindex_all' `treatment' `covariates' [`weight' `exp'] `if', vce( `std_errors' )
    local b = _b[`treatment']
    if ( `fold_disagree' ) {
			local b = 0
		}
    local se = _se[`treatment']
    mata: b_by_iter[`f_iter', 1] = `b'  // Store the optimus coefficient
    local t_0 = `b' / `se'

    mata: N_0 = `e(N)'  // Store the sample size

    if ("`rw'" != "") {
      local rw_index = 2  // index 1 is the optimus which is handled separately
      foreach y of varlist `rw' {
        quietly reg `y' `treatment' `covariates' [`weight' `exp'] `if', vce( `std_errors' )
        local b = _b[`treatment']
        local se = _se[`treatment']
        local t_`rw_index' = `b' / `se'
        local rw_index = `rw_index' + 1
      }
    }
    if (`onesided') {
      if ("`sign'" == "pos") {
        mata: p_actual[1, 1] = 1-normal(`t_0')
        if ("`rw'" != "") {
          forval rw_index = 2(1)`vars_rw' {
              mata: p_actual[1, `rw_index'] = 1-normal(`t_`rw_index'')
          }
        }
      }
      else {
        mata: p_actual[1, 1] = 1-normal(-1*`t_0')
        if ("`rw'" != "") {
          forval rw_index = 2(1)`vars_rw' {
              mata: p_actual[1, `rw_index'] = 1-normal(-1*`t_`rw_index'')
          }
        }
      }
    }
    else {
      mata: p_actual[1, 1] = 2*(1-normal(abs(`t_0')))
      if ("`rw'" != "") {
        forval rw_index = 2(1)`vars_rw' {
            mata: p_actual[1, `rw_index'] = 2*(1-normal(abs(`t_`rw_index'')))
        }
      }
    }

    *********************************************************************************
    * Now calculate the RI p-value for this split
    *********************************************************************************

    forval ri = 1(1)`ri_iterations' {

      foreach x in `cell_centile' `rdraw' `treatment_cells' `treatment_bs' `rank' `max' `eligible' `treatment_r_permuted' `treatment_prob' {
        capture drop `x'
      }

      if ( `onesided' ) {
  			forvalues f = 1(1)`folds' {
  				quietly replace `optindex_pos_`f'' = .
  			}
  			forvalues f = 1(1)`folds' {
  				quietly replace `optindex_neg_`f'' = .
  			}
  		}
  		else {
  			forvalues f = 1(1)`folds' {
  				quietly replace `optindex_`f'' = .
  			}
  		}

      // Generate a permuted treatment assignment
      quietly gen byte `treatment_bs' = `treatment' if `all_outcomes'
      if (`strat_nulltreat_fold') {
        if ( "`stratify'" != "" ) {
            quietly egen int `treatment_cells' = group( `stratify' `fold') `if' & `treatment_bs' != .
        }
        else {
            quietly gen int `treatment_cells' = `fold' `if' & `treatment_bs' != .
        }
      }
      else {
        if ( "`stratify'" != "" ) {
            quietly egen int `treatment_cells' = group( `stratify') `if' & `treatment_bs' != .
        }
        else {
            quietly gen int `treatment_cells' = 1 `if' & `treatment_bs' != .
        }
      }

      sum `treatment_cells', meanonly
      local total_cells = r(max)
      if ( `total_cells' > 500 ) {
        if (`strat_nulltreat_fold') {
          dis "WARNING: You have specified that treatment probability varies across more than 500 cells! Consider setting strat_nulltreat_fold(0)."
        }
        else {
          dis "WARNING: You have specified that treatment probability varies across more than 500 cells!"
        }
      }
      quietly bysort `treatment_cells': egen double `treatment_prob' = mean( `treatment_bs' ) `if' & `treatment_bs' != . & `all_covariates'
      sort `sort_order'

      cluster_random_number `treatment_cluster', output_var(`rdraw')  // Random assignment to each cluster (ensures same treatment status of clusters)
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
      quietly predict `treatment_r_permuted' if e(sample), resid

      // Now calculate the optimus matrix coefficient on the permuted treatment assignment
      mata: votes = J( `folds', 1, . ) // Will fill in with 1 for positive and 0 for negative
  		local fold_disagree = 0
      forvalues f = 1(1)`folds' {
        // Estimate coefficients and covariance matrix for optimus
        gen_coeff_and_covar	`varlist' `if' & `fold' != `f', treatment_r(`treatment_r_permuted') covariates(`covariates') treatment_cluster(`cluster') ///
          rdraw(`rdraw') all_covariates(`all_covariates') fv(`f') prescreen_cutoff(`prescreen_cutoff') bootstrap_cov(`bootstrap_cov') ///
          cov_bootstrapreps(`cov_bootstrapreps') vce(`std_errors') permuted(1)

        // Get the optimal index
        capture mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
          min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `cov_shrinkage', ///
          `unweighted', cov_V, `concen_weight' )

        if ( _rc == 430 ) {
          dis "Mata optimizer failed to find a solution on rep `r'. Trying slightly different covariance shrinkage factor for this rep."
          local temp_cov_shrink = `cov_shrinkage' * 0.9
          capture mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
            min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `temp_cov_shrink', ///
            `unweighted', cov_V, `concen_weight' )
          if ( _rc == 430 ) {
            dis "Mata optimizer failed to find a solution on rep `r'. Reverting to unweighted index for this rep."
            mata: index_output = summary_power_calc( beta, omega, group_selector( select_matrix, ///
              min( ( length( beta ), 10 ) ) ), `expected_crit_val', ( length( beta ) > 10 ), `cov_shrinkage', 1, ///
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
  					quietly replace `optindex_`sign'_`f'' = 0 `if' & ( `fold' != `f' )
  				}
  			}
  			else {
  				mata: index_selector = index_output[2..rows(index_output),1+negative]
  				mata: st_view( optindex_`f' = ., ., "`optindex_`f''" )
  				mata: optindex_`f'[.] = Y_`f' * index_selector / colsum( index_selector ) // Generate the optimal index for fold f
  				quietly replace `optindex_`f'' = 0 `if' & ( `fold' != `f' )
  			}

      }

      capture drop `optindex_all'
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
      local t_0 = `b' / `se'

      if ("`rw'" != "") {
        local rw_index = 2  // index 1 is the optimus which is handled separately
        foreach y of varlist `rw' {
          quietly reg `y' `treatment_bs' `covariates' [`weight' `exp'] `if', vce( `std_errors' )
          local b = _b[`treatment_bs']
          local se = _se[`treatment_bs']
          local t_`rw_index' = `b' / `se'
          local rw_index = `rw_index' + 1
        }
      }

      if (`onesided') {
        if ("`sign'" == "pos") {
          mata: p_student[`ri', 1] = 1-normal(`t_0')
          if ("`rw'" != "") {
            forval rw_index = 2(1)`vars_rw' {
                mata: p_student[`ri', `rw_index'] = 1-normal(`t_`rw_index'')
            }
          }
        }
        else {
          mata: p_student[`ri', 1] = 1-normal(-1*`t_0')
          if ("`rw'" != "") {
            forval rw_index = 2(1)`vars_rw' {
                mata: p_student[`ri', `rw_index'] = 1-normal(-1*`t_`rw_index'')
            }
          }
        }
      }
      else {
        mata: p_student[`ri', 1] = 2*(1-normal(abs(`t_0')))
        if ("`rw'" != "") {
          forval rw_index = 2(1)`vars_rw' {
              mata: p_student[`ri', `rw_index'] = 2*(1-normal(abs(`t_`rw_index'')))
          }
        }
      }

    }
    // The unadjusted (for multiple hypothesis testing) p-value of the optimus index is just the share of the time p_actual > p_student (permuted)
    // This already accounts for 1 or 2-sided tests based on how p_actual and p_student are calculated
    mata: p_greater = p_actual[1, 1] :> p_student[ , 1]
    mata: p_by_iter[`f_iter', 1] = mean(p_greater)  // Fraction of time actual p-value of optimus index exceeds null

    // Calculate the Romano-Wolf vector of p-values for this split, if relevant
    if ("`rw'" != "") {
      mata: temp = stepdown(p_actual, p_student)  // The second row has the RW p-values
      mata: p_rw_by_iter[ , `f_iter'] = temp[2, ]'
    }
  }

  // Now generate statistics to return
  mata: p_w_index = matrix_median(p_by_iter)
  mata: p = p_w_index[1, 1]

  // Record the index (or indices) of the CV draw associated with the median p value
  if (mod(`fold_iterations', 2) == 0) {  // Then the median is associated with 2 index values
    local med_indices = 2
    mata: med_index_1 = p_w_index[1, 2]
    mata: med_index_2 = p_w_index[1, 3]
  }
  else {
    local med_indices = 1
    mata: med_index = p_w_index[1, 2]
  }

  // Now extract the RW p-values associated with the median unadjusted p-value split
  if ("`rw'" != "") {
    if (`med_indices' == 1) {
      mata: p_rw = p_rw_by_iter[, med_index]'
    }
    else {
      mata: p_rw = 0.5*p_rw_by_iter[, med_index_1]' + 0.5*p_rw_by_iter[, med_index_2]'
    }
  }

  // Chernozhukov et al. (2018) and Romano and DiCiccio (2019) bound on p-value across multiple fold iterations if no seed specified
  if (`fold_seed' == 0) {
    mata: p = 2*p
    if ("`rw'" != "") {
      mata: p_rw = 2*p_rw
    }
  }

  // Diplay version of RW p-values
  if ("`rw'" != "") {
    mata: p_rw_rounded = round(p_rw, 0.01)
    mata: p_rw_display = J(2, `vars_rw', "")
    local full_rw_names = "optimus `rw'"
    tokenize "`full_rw_names'"
    forval rw_index = 1(1)`vars_rw' {
      mata: p_rw_display[1, `rw_index'] = "``rw_index''"
      mata: st_numscalar("tmp", p_rw_rounded[1, `rw_index'])
      local tmp = tmp
      mata: p_rw_display[2, `rw_index'] = "`tmp'"
    }
  }

  if (`med_indices' == 1) {
    mata: median_b = b_by_iter[med_index, 1]
    mata: median_weights = average_weights[,med_index]
    mata: median_any_weight = any_weight[med_index]
    mata: median_nontrivial_weight = nontrivial_weight[med_index]
    if (`onesided') {
      mata: median_pos = pos_by_iter[med_index]
    }
  }
  else {
    mata: median_b = 0.5*b_by_iter[med_index_1, 1] + 0.5*b_by_iter[med_index_2, 1]
    mata: median_weights = 0.5*average_weights[,med_index_1] + 0.5*average_weights[,med_index_2]
    mata: median_any_weight = 0.5*any_weight[med_index_1] + 0.5*any_weight[med_index_2]
    mata: median_nontrivial_weight = 0.5*nontrivial_weight[med_index_1] + 0.5*nontrivial_weight[med_index_2]
    if (`onesided') {
      mata: median_pos = 0.5*pos_by_iter[med_index_1] + 0.5*pos_by_iter[med_index_2]
    }
  }

  mata: mean_b = mean(b_by_iter)
  if (`onesided') {
    mata: mean_pos = mean(pos_by_iter)
  }
  mata: mean_weights = mean(average_weights')
  mata: mean_any_weight = mean(any_weight')
  mata: mean_nontrivial_weight = mean(nontrivial_weight')

  // Diplay version of median weights
  mata: weights_rounded = round(median_weights, 0.01)
  mata: weights_display = J(2, `vars_in_group', "")
  tokenize "`varlist'"
  forval ind = 1(1)`vars_in_group' {
    mata: weights_display[1, `ind'] = "``ind''"
    mata: st_numscalar("tmp", weights_rounded[`ind'])
    local tmp = tmp
    mata: weights_display[2, `ind'] = "`tmp'"
  }

  use `_torestore', clear

  mata: st_numscalar("e(_b)", median_b)  // Returning directly to e(b) fails since protected
  tempname b
  matrix `b' = e(_b)
  ereturn post `b'

  mata: st_numscalar("e(p)", p)
  mata: st_numscalar("e(N)", N_0)
  mata: st_numscalar("e(_b)", median_b)  // Re-create to display as a scalar
  if (`onesided') {
      mata: st_numscalar("e(pos)", median_pos)
      mata: st_numscalar("e(mean_pos)", mean_pos)
      mata: st_matrix("e(pos_by_split)", pos_by_iter)
  }
  mata: st_numscalar("e(any_weight)", median_any_weight)
  mata: st_numscalar("e(nontrivial_weight)", median_nontrivial_weight)
  mata: st_numscalar("e(mean_b)", mean_b)
  mata: st_numscalar("e(mean_any_weight)", mean_any_weight)
  mata: st_numscalar("e(mean_nontrivial_weight)", mean_nontrivial_weight)

  mata: st_matrix("e(weights)", median_weights)
  mata: st_matrix("e(mean_weights)", mean_weights)

  mata: st_matrix("e(p_by_split)", p_by_iter[,1])
  mata: st_matrix("e(b_by_split)", b_by_iter)
  mata: st_matrix("e(weights_by_split)", average_weights)
  mata: st_matrix("e(any_weight_by_split)", any_weight)
  mata: st_matrix("e(nontrivial_weight_by_split)", nontrivial_weight)

  if ("`rw'" != "") {
    mata: st_matrix("e(p_rw)", p_rw)
    dis "Index coefficient: `e(_b)'"
    if (`onesided') {
      dis "Index positive: `e(pos)'"
    }
    dis "N: `e(N)'"
    dis "Index unadjusted p-value: `e(p)'"
    dis "Romano-Wolf adjusted p-values"
    mata: p_rw_display
    dis "Average number of weights > 1/H across folds: `e(nontrivial_weight)'"
    dis "Average weights across folds:"
    mata: weights_display
  }
  else {
    dis "Index coefficient: `e(_b)'"
    if (`onesided') {
      dis "Index positive: `e(pos)'"
    }
    dis "N: `e(N)'"
    dis "Index p-value: `e(p)'"
    dis "Average number of weights > 1/H across folds: `e(nontrivial_weight)'"
    dis "Average weights across folds:"
    mata: weights_display
  }

end



// Define helper functions (non-mata)
// Program to generate random numbers at the treatment cluster level
// Takes as input the cluster id variable(s)
cap program drop cluster_random_number
program define cluster_random_number
	syntax varlist, output_var(string)
	tempvar order select random
	quietly gen long `order' = _n
	quietly egen byte `select' = tag( `varlist' )
	quietly gen double `random' = runiform() if `select'
	quietly bysort `varlist': egen double `output_var' = mean( `random' )
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
		all_covariates(varname) fv(integer) vce(string) [covariates(varlist) prescreen_cutoff(real 1.2) bootstrap_cov(integer 0) cov_bootstrapreps(real 10) ///
		inflation_factor_base(real 0.3) inflation_factor_rho(real 0.3) permuted(integer 0)]
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
	else if ( `vars_in_group' < 12 ) {
		local prescreen_cutoff = min( 0.5, `prescreen_cutoff' )
	}
	foreach outcome of varlist `varlist' {
		quietly reg `outcome'_sd `treatment_r' `covariates' [`weight' `exp'] `if', vce(`vce') // Screen out small t-stats and drop dep vars with almost no variation to make the SUR system manageable
		if ( `++count_catch' == 1 ) {
			mata: temp_string_catch[1] = "`outcome'" // Record first outcome name
			local outcome1 `outcome'
		}
		if ( abs( _b[`treatment_r'] / _se[`treatment_r'] ) > `prescreen_cutoff' & e(rmse) > 0.05 ) {
			local prescreen_level "`prescreen_level' `outcome'_sdr" // List of outcomes (SDR)
			local screened_group_level "`screened_group_level' `outcome'" // List of outcomes (raw)
			mata: temp_string_full[`++count_full'] = "`outcome'" // Fill in vector of outcome names
		}
	}
	if ( `count_full' == 0 ) {
    if (`permuted' == 0) {  // Only display warning if on actual data
  		dis "No outcome cleared the t-threshold of `prescreen_cutoff' for level index!"
    }
		mata: bname = temp_string_catch[1] // Just take the first outcome if nothing looks good
		local prescreen_level "`prescreen_level' `outcome1'_sdr" // First outcome (SDR)
		local screened_group_level "`screened_group_level' `outcome1'" // First outcome (raw)

    mata: included_vars = J(1, `vars_in_group', 0)  // For assigning 0 weights to all excluded variables in recording average weights
    tokenize "`varlist'"
    forval i = 1(1)`vars_in_group' {
      if "``i''" == "`temp_name'" {
        mata: included_vars[`i'] = 1
      }
    }
	}
	else {
		mata: bname = temp_string_full[1..`count_full'] // Final vector of outcome names
    mata: included_vars = J(1, `vars_in_group', 0)  // For assigning 0 weights to all excluded variables in recording average weights
    tokenize "`varlist'"
    forval j = 1(1)`count_full' {
      mata: st_local("temp_name", bname[`j'])
      forval i = 1(1)`vars_in_group' {
        if "``i''" == "`temp_name'" {
          mata: included_vars[`i'] = 1
        }
      }
    }
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
	foreach outcome of varlist `group_varlist' {
		local estimation_cmd "`estimation_cmd' (`outcome'_sdr `treatment_r')" // Construct the preliminary sureg command
		mata: y`++group_count' = st_data( ., "`outcome'_sdr" ) // Get the LHS vars into Mata
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
	local zeros_needed = round( e(k) / `group_count' ) - 1
	mata: select_colsrows = J( 1, `group_count', ( 1 , J( 1, `zeros_needed', 0 ) ) )
	mata: beta = select( st_matrix( "e(b)" ), select_colsrows )
	mata: omega = select( select( st_matrix( "e(V)" ), select_colsrows ), select_colsrows' )
	local N = e(N)
	mata: N = `N'
	local depvars = e(depvar)
	if ( `bootstrap_cov' ) {
		// Bootstrap the covariance matrix std errors if necessary
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
function summary_power_calc( real rowvector beta, real matrix covariance, real matrix select_variables, real scalar crit_val, real scalar fast, real scalar cov_shrinkage_factor, real scalar unweighted, real matrix cov_covariance, real scalar concen_weight )
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
			se_index[i] = sqrt( ( select_variables_internal[i, .] * V * select_variables_internal[i, .]' ) ///
				/ ( row_size[i] ^ 2 + ( row_size[i] == 0 ) * 0.01 ) + ( row_size[i] == 0 ) * 0.01 )
		}
		weights_internal = select_variables_internal :/ ( row_size + ( row_size :== 0 ) * 0.01 )
		power_pos = ( beta_avg :/ se_index ) - concen_weight * rowsum( weights_internal :^2 )
		power_neg = -( beta_avg :/ se_index ) - concen_weight * rowsum( weights_internal :^2 )
		max_power_index_pos = .
		max_power_index_neg = .
		ties = .
		max_power_index_pos = max_power_index_pos[1]
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
				// Positive direction and no positive betas
				temp_power[1] = 0 // Zero out power in that direction
				maxindex( beta, 1, max_power_index, ties ) // Chooose the least negative beta (shouldn't really matter)
				max_power_index = max_power_index[1]
				w[1,max_power_index] = 1 // Assign all weight to that beta
			}
			else if ( !negative * ( sum( beta :>= 0 ) == 1 ) ) {
				// Positive direction and one positive beta
				maxindex( beta, 1, max_power_index, ties ) // Chooose the positive beta
				max_power_index = max_power_index[1]
				w[1,max_power_index] = 1 // Assign all weight to that beta
				temp_power[1] = beta[max_power_index] / sqrt( V[max_power_index,max_power_index] )
			}
			else if ( negative * !sum( beta :< 0 ) ) {
				// Negative direction and no negative betas
				temp_power[2] = 0 // Zero out power in that direction
				minindex( beta, 1, max_power_index, ties ) // Chooose the least positive beta (shouldn't really matter)
				max_power_index = max_power_index[1]
				w[2,max_power_index] = 1 // Assign all weight to that beta
			}
			else if ( negative * ( sum( beta :< 0 ) == 1 ) ) {
				// Negative direction and one negative beta
				minindex( beta, 1, max_power_index, ties )
				max_power_index = max_power_index[1]
				w[2,max_power_index] = 1  // Assign all weight to that beta
				temp_power[2] = abs( beta[max_power_index] / sqrt( V[max_power_index,max_power_index] ) )
			}
			else {
				S = optimize_init()
				optimize_init_evaluator( S, &optimus() )
				optimize_init_evaluatortype( S, "d1" )
				b = beta
				optimize_init_argument( S, 1, b )
				optimize_init_argument( S, 2, V )
				optimize_init_argument( S, 3, crit_val )
				optimize_init_argument( S, 4, negative )
				optimize_init_argument( S, 5, concen_weight )
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
void optimus(todo, w, b, V, crit_val, negative, con_wgt, power, g, H)
{
	if ( !negative ) {
		power = ( b * w' ) / sqrt( w * V * w' )  -
			.001 * sum( ( 2 * w :- 1 ) :^ 200 ) - con_wgt * sum( w:^2 )
		if ( todo >= 1 ) {
			end_val = length( b )
			Vw = V * w'
			for ( i = 1; i <= end_val; i++ ) {
				g[i] =  ( b[i] / sqrt( w * V * w' ) - b * w' * Vw[i] * ( w * V * w' ) ^ -1.5 ) -
					0.4 * ( 2 * w[i] - 1 ) ^ 199 - con_wgt * 2 * w[i]
			}
		}
	}
	else {
		power = -( b * w' ) / sqrt( w * V * w' )  -
			.001 * sum( ( 2 * w :- 1 ) :^ 200 ) - con_wgt * sum( w:^2 )
		if ( todo >= 1 ) {
			end_val = length( b )
			Vw = V * w'
			for ( i = 1; i <= end_val; i++ ) {
				g[i] =  -( b[i] / sqrt( w * V * w' ) - b * w' * Vw[i] * ( w * V * w' ) ^ -1.5 ) -
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

// Define a function to calculate the median p-value and index
function matrix_median(real matrix mat)
{
  sorted_mat = sort(mat, 1)
  tmp = rows(sorted_mat)/2
  if (mod(tmp, 1) == 0) {
    med = 0.5*(sorted_mat[tmp, 1] + sorted_mat[tmp + 1, 1])
    med_index = (sorted_mat[tmp, 2], sorted_mat[tmp + 1, 2])
  }
  else {
    med = sorted_mat[ceil(tmp), 1]
    med_index = (ceil(tmp))
  }
  return (med, med_index)
}

// Function to calculate Romano-Wolf p-values from actual and permuted p-values calculated via RI
// Output is a 2 x H matrix sorted in original order with RW p-vals in second row
function stepdown(real matrix p_actual, real matrix p_student)
{
	// Organize the permuted statistics
	p_total = ( range( 1, cols( p_actual ), 1 )' \ p_actual \ p_student )
	p_sort = sort( p_total', 2 )' // (2 + rw_reps) x H sorted matrix with original order as first row, actual p-val as second
	endrow = rows( p_sort )
	endcol = cols( p_sort )
	p_actual = p_sort[2, .] // Resorted vector of actual p-vals
	// Get the min p-val for each rep, with stepdown below
	p_min = rowmin( p_sort[| 3, 1 \ endrow, endcol |] ) // Syntax note: can't have space between a bracket and a pipe!
	for ( h = 2; h <= endcol; h++ ) {
		p_min = ( p_min, rowmin( p_sort[| 3, h \ endrow, endcol |] ) )
	}
	p_greater = p_actual :> p_min // Mark cases in which actual p-val exceeds null permuted p-val
	rw_pval = colsum( p_greater ) / rows( p_greater ) // Compute fraction of times actual p-val exceeds null permuted p-val
	results = ( p_sort[| 1, 1 \ 1, endcol |] \ rw_pval ) // 2 x H matrix with original order in first row and RW p-vals in second row
	// Enforce monotonicity
	for ( h = 2; h <= endcol; h++ ) {
		results[2, h] = rowmax( results[2, h-1..h] )
	}
	results = sort( results', 1 )' // 2 x H matrix sorted in original order with RW p-vals in second row
	return( results )
}

end
