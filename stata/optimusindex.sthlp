{smcl}
{viewerjumpto "Syntax" "optimusindex##syntax"}{...}
{viewerjumpto "Description" "optimusindex##description"}{...}
{viewerjumpto "Options" "optimusindex##options"}{...}
{viewerjumpto "Stored results" "optimusindex##results"}{...}

{* June 15th, 2022}{...}
{hline}
help for {hi:optimusindex}
{hline}

{title:Optimus index calculation for use in highly powered analysis plans (Anderson and Magruder 2022)}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:optimusindex} {it:outcomes} [{it:aweight}]
[{cmd:if} {it:exp}]{cmd:, {cmd:treatment(}{it:varname}{cmd:)} }
{bind:[{cmd:cluster(}{it:varlist}{cmd:)}}
{cmd:stratify(}{it:varlist}{cmd:)}
{cmd:covariates(}{it:varlist}{cmd:)}
{cmd:unweighted(}{it:binary integer}{cmd:)}
{cmd:onesided(}{it:binary integer}{cmd:)}
{cmd:bootstrap_cov(}{it:binary integer}{cmd:)}
{cmd:cov_bootstrapreps(}{it:integer}{cmd:)}
{cmd:folds(}{it:integer}{cmd:)}
{cmd:fold_seed(}{it:integer}{cmd:)}
{cmd:fold_iterations(}{it:integer}{cmd:)}
{cmd:ri_iterations(}{it:integer}{cmd:)}
{cmd:prescreen_cutoff(}{it:real}{cmd:)}
{cmd:cov_shrinkage(}{it:real}{cmd:)}
{cmd:concen_weight(}{it:real}{cmd:)}
{cmd:rw(}{it:varlist}{cmd:)}
{bind:{cmd:strat_nulltreat_fold(}{it:binary integer}{cmd:)}]}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Primary}
{synopt: {cmdab:treatment(varname)}} Variable identifying random treatment assignment {p_end}
{synopt: {cmdab:cluster(varlist)}} Cluster variables used for calculating p-values. {p_end}
{synopt: {cmdab:stratify(varlist)}} Denotes any variables used to stratify treatment assignment. It is important to
specify these if they were used in randomization to get accurate p-values. {p_end}
{synopt: {cmdab:covariates(varlist)}} List of covariates to control for. {p_end}
{synopt: {cmdab:unweighted(binary integer)}} Specify 1 to calculate an unweighted optimus index. {p_end}
{synopt: {cmdab:onesided(binary integer)}} Default is 1. Specify 0 to calculate a two-sided hypothesis test. {p_end}

{syntab:Covariance matrix estimation}
{synopt: {cmdab:bootstrap_cov(binary integer)}} Specify 1 to bootstrap covariance matrix estimation. {p_end}
{synopt: {cmdab:cov_bootstrapreps(integer)}} Only relevant if bootstrap_cov set to 1. Specifies the number of bootstrap reps used. {p_end}

{syntab:Sample splitting}
{synopt: {cmdab:folds(integer)}} Number of folds. Default is 5. {p_end}
{synopt: {cmdab:fold_seed(integer)}} Specifies the seed used for sample splitting. This should only be used if a seed was pre-specified.
Otherwise set a seed using the set seed command before running this program. {p_end}
{synopt: {cmdab:fold_iterations(integer)}} Number of sample splits. The default is 100. If fold_seed is specified,
this is overridden and set to 1. {p_end}
{synopt: {cmdab:ri_iterations(integer)}} Number of iterations for randomization inference. The default is 500.
The total number of iterations will be ri_iterations*fold_iterations.

{syntab:Hyperparameters}
{synopt: {cmdab:prescreen_cutoff(real)}} Important parameter for ensuring nonsingular covariance matrices with sureg.  Positive real number with a default of 1.2.
Variables with a t-stat below this threshold are assigned a weight of 0 in the optimus index. {p_end}
{synopt: {cmdab:cov_shrinkage(real)}} Real number between 0 and 1. A value of 0 does not shrink the
covariance matrix, while 1 applies the full empirical Bayesian shrinkage estimate. Default is 0.5. {p_end}
{synopt: {cmdab:concen_weight(real)}} Positive real number. Penalty weight on index concentration. Default is 0.5. {p_end}

{syntab:Additional options}
{synopt: {cmdab:rw(varlist)}} The input is a subset of the variables used in the index. If included,
a Romano-Wolf correction for multiple hypothesis testing is calculated across these variables and the optimus index. {p_end}
{synopt: {cmdab:strat_nulltreat_fold(binary integer)}} Default is 1. If 0, permuted treatment assignments
are not stratified by fold. It is recommended to keep this option on unless the number of stratification variables is high. {p_end}

{marker description}{...}
{p} {title:Description}

{p}{cmd:optimusindex} calculates the optimus index as defined in Anderson and Magruder (2022) {it:Highly Powered Analysis Plans}.
The algorithm calculates a power-maximizing index over a family of related outcomes using k-fold cross-validation.
This program is intended for use with data from a randomized controlled trial (RCT).

{p} If a seed for sample splitting is pre-specified and passed via the fold_seed option, then a single sample split
is generated using this seed. This is the most computationally efficient option and generally
produces lower p-values. If no seed is pre-specified, then multiple independent sample splits
are taken. Following Chernozhukov et al. (2018) and Romano and DiCiccio (2019), we bound the
p-value as 2 times the median p-value across sample splits.

{p} The command can return a p-value, calculated using randomization inference, that is
unadjusted for multiple hypothesis testing or a vector of p-values that uses the Romano-Wolf
stepdown procedure to control the familywise error rate (FWER) (Romano and Wolf 2015a,b; 2016).
The p-values are constructed by first taking either a one or two-sided p-value from
a regression of the optimusindex on treatment, controlling for covariates. We then
calculate the same p-values for each permuted draw of the data and calculate the
RI p-value as the share of the time the permuted p-values are smaller than the
p-value estimated using actual treatment assignment.
The multiple testing adjusted p-values may be calculated over an arbitrary subset of outcomes entered into the
optimus index calculation using the {cmd:rw(}{it:varlist}{cmd:)} option. The program
will return the adjusted p-values for these variables as well as the optimus index.

{p} Those wishing to use the optimus index for gatekeeping purposes will generally
wish to use unadjusted p-values. Researchers that wish to use other methods to control
the FWER or false discovery rate may generally apply these methods ex-post
based on the unadjusted p-value since most common methods other than the Romano-Wolf correction only require unadjusted
p-values as an input.

{p} A limited set of results are displayed to the Stata output. These outputs plus
a larger set of results are also returned to e-class macros. These may be easily used
in a program or written to an output such as a LaTex or Excel document for display. Details
about the returned values are presented in the section {help optimusindex##results:Stored results}.

{p} The command uses Mata for calculations. To avoid overwriting any Mata results stored by the user,
it does not clear Mata at the end of the command. However, you may wish to clear Mata at the
end of the run to reduce memory usage.

{marker options}{...}
{p} {title:Options}

{dlgtab:Primary}

{phang}
{opth treatment(varname)} defines the variable denoting treatment assignment. If there are multiple treatment variables,
they should be combined into a single categorical variable before running the command. This must be a numeric
variable.

{phang}
{opth cluster(varlist)} is an optional parameter capturing dependence in the data. If a list of categorical variables is passed,
p-values will be clustered by the groups. Specifically, permuted treatment assignments are constant across clusters.
Clustered standard errors are also used for the pre-screen cutoff if this option is passed. If not cluster variables
are included, HC3 standard errors are used for the pre-screen cutoff and other internal steps.

{phang}
{opth stratify(varlist)} specifies a list of categorical variables that treatment assignment
was stratified on. This parameter is optional to account for designs that do not stratify
treatment assignment, but {it: must} be specified if randomization was stratified.

{pmore}
Fold assignments are stratified by treatment and by each of the strata variables.
Furthermore, the permuted treatment assignments for randomization inference are
stratified on these variables. This is important to obtain correct p-values since
the RI procedure must mirror that used to assign initial treatment assignment.

{phang}
{opth covariates(varlist)} specifies control variables to include in all regressions.
If randomization was stratified, the strata variables should by included in this list.
Covariates can include continuous variables or categorical variables (e.g. fixed effects),
but categorical variables should be passed as a series of dummy variables. High-dimensional
fixed-effects are not currently supported.

{phang}
{opth unweighted(binary integer)} controls whether an unweighted optimus index is
calculated. The default is 0, meaning a weighted index is calculated. If an unweighted
index is specified, all (standardized) outcomes receive a weight of 0 or 1. This option
can be computationally easier if the number of outcomes in the index is large. This
differs from simply using the standardized mean of the index outcomes since the
algorithm selects the set of variables that maximizes power.

{phang}
{opth onesided(binary integer)} determines whether a one or two-sided test is conducted.
The default is 0, meaning a one-sided test is conducted.
{dlgtab:Primary}

{dlgtab:Covariance estimation}
{phang}
These options determine how the covariance matrix of index variables is calculated. This covariance
matrix is used for constructing the index.

{phang}
{opth bootstrap_cov(binary integer)} controls whether the covariance matrix of index
variables is bootstrapped. The default is 0, which results in the covariance matrix
outputted from sureg being used in the index construction. If enabled, the Wild bootstrap
is used to estimate the covariance matrix. The bootstrap can be desirable because
it better accounts for dependence (clusters) in the data. However, in practice we
found this is generally not important, and bootstrapping significantly raises computational
costs (total bootstrap calculations is bootstrap_reps*folds*(1_ri_iterations + 1)*fold_iterations).

{phang}
{opth cov_bootstrapreps(integer)} specifies the number of bootstrap draws used for estimating the covariance matrix,
which is only relevant if the previous option is enabled. The default is 10. The total number of
bootstrap draws is given by bootstrap_reps*folds*(1_ri_iterations + 1)*fold_iterations.

{dlgtab:Sample splitting}

{phang}
{opth folds(integer)} determines the number of folds. The default is 5. With a large
sample, more folds can help improve the precision of the results, but we've generally
found little gain from a large number of folds.

{phang}
{opth fold_seed(integer)} specifies the seed used to split the sample. In general,
this should only be used if a seed was specified in a pre-analysis plan. If this option
is enabled, only a single split of the data is constructed using this seed, and p-values are
not adjusted for multiple splits.

{pmore}
This option should {it: not} be used if one simply wishes to have replicable code. Rather,
the seed should be set before running this command using the set seed command.

{phang}
{opth fold_iterations(integer)} controls the number of sample splits. The default is 100.
If fold_seed is specified, this is set to 1 by the program. For each fold iteration,
an independent sample split is constructed and the index, index coefficient, and p-values are
calculated. Following Chernozhukov et al. (2018) and Romano and DiCiccio (2019), we then return
2 times the median p-value across folds. The treatment coefficient and other parameters associated
with this median draw are also returned, with details provided below. Users may also view the mean
coefficient and weights across all splits.

{phang}
{opth ri_iterations(integer)} determined the number of permuted treatment draws used for
randomization inference. The default is 500.

{dlgtab:Hyperparameters}

{phang}
{opth prescreen_cutoff(real)} should be a strictly positive number, with a default of 1.2.
For each fold, we regress each outcome in the index on treatment controlling for covariates
with either clustered or HC3 standard errors using observations from all other folds. If
the two-sided t-stat falls below this value, the variable is excluded from the index for that
fold (weighted to 0), and any variable clearing the cutoff are used for index construction.

{pmore}
This parameter is important to ensure that covariance matrices are non-singular with sureg.
The number of variables used in index calculation is limited to ensure that the code can execute.
If the number of covariates is below 6, this is ignored and all variable are included.
If the number of index variables is below 12, we keep all variables with a t-stat above half the
prescreen cutoff.

{phang}
{opth cov_shrinkage(real)} is a number between 0 and 1. The default
value is 0.5. This controls how much
covariance terms in the covariance matrix terms are shrunk towards 0. This
prevents outcomes with small negative covariances due to noise from being
weighted too strongly. Note that variance terms are not adjusted.

{pmore}
If 0 is passed, no regularization is applied to the covariance terms. If 1 is specified,
empirical Bayes shrinkage is used to regularize the covariance terms. For values between
0 and 1, a weighted average of these two covariance matrices is calculated.

{phang}
{opth concen_weight(real)} is a positive number with a default of 0.5. This is the
value multiplying the HHI which penalizes a lack of complexity in the optimus index.

{dlgtab:Additional options}

{phang}
{opth rw(varlist)} calculates Romano-Wolf adjusted p-values across the optimus
index and each of the variables passed. The variables passed should generally be
a susbet of the index variables. The same test used for the optimus index is
used for each of the variables included (e.g. a one or two-sided test controlling
for the same set of covariates).

{phang}
{opth strat_nulltreat_fold(binary integer)} is set to 1 by default. This
stratifies the permuted treatment assignment by fold, which mirrors the fact that
fold assignment is stratified by treatment. This generally produces better results,
but if a large number of strata are present it may be set to 0.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:optimusindex} stores the following in {cmd:e()}:

{pstd}

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations{p_end}

{syntab:Scalars - Median split}
{synopt:{cmd:e(p)}}optimus index p-value (unadjusted for multiple testing if RW used){p_end}
{synopt:{cmd:e(_b)}}optimus index coefficient{p_end}
{synopt:{cmd:e(pos)}}Index sign positive, if one-sided test only{p_end}
{synopt:{cmd:e(any_weight)}}Average variables receiving any weight{p_end}
{synopt:{cmd:e(nontrivial_weight)}}Average variables with weight > 1/(# outcomes){p_end}

{syntab:Scalars - Mean across split}
{synopt:{cmd:e(mean_b)}}optimus index coefficient{p_end}
{synopt:{cmd:e(mean_pos)}}Index sign positive, if one-sided test only{p_end}
{synopt:{cmd:e(mean_any_weight)}}Average variables receiving any weight{p_end}
{synopt:{cmd:e(mean_nontrivial_weight)}}Average variables with weight > 1/(# outcomes){p_end}

{syntab:Matrices - Median split}
{synopt:{cmd:e(p_rw)}}Romano-Wolf adjusted p-values{p_end}
{synopt:{cmd:e(b)}}optimus index coefficient{p_end}
{synopt:{cmd:e(weights)}}average variable weights across folds{p_end}

{syntab:Matrices - Average across splits}
{synopt:{cmd:e(mean_weights)}}average variable weights across folds{p_end}

{syntab:Matrices - Results by split}
{synopt:{cmd:e(p_by_split)}}optimus index p-value by split{p_end}
{synopt:{cmd:e(b_by_split)}}optimus index coefficient by split{p_end}
{synopt:{cmd:e(weights_by_split)}}optimus index weights by split{p_end}
{synopt:{cmd:e(any_weight_by_split)}}optimus index number of variables with any weight by split{p_end}
{synopt:{cmd:e(nontrivial_weight_by_split)}}optimus index weights > 1/(# outcomes) by split{p_end}


{marker contact}{...}
{title:Author}

{pstd}Michael Anderson, Jeremy Magruder, and Grady Killeen{break}
Email: {browse "mailto:gkilleen@berkeley.edu":gkilleen@berkeley.edu}
{p_end}


{title:References}

{p 0 4}Anderson, Michael L and Jeremy Magruder (2022). Highly Powered Analysis Plans. {it: NBER Working Paper}.

{p 0 4}Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen, Whitney Newey, James Robins (2018).
Double/debiased machine learning for treatment and structural parameters. {it: The Econometrics Journal} 21:1, https://doi.org/10.1111/ectj.12097.

{p 0 4}DiCiccio, Cyrus and Joseph Romano (2019). Multiple Data Splitting for Testing. Technical Report No. 2019-03,
{it: Stanford University Department of Statistics}.

{p 0 4}Romano, Joseph P and Michael Wolf (2005a). Exact and Approximate Stepdown Methods for Multiple Hypothesis Testing,
{it: Journal of the American Statistical Association}, 100:469, 94-108, DOI: 10.1198/016214504000000539.

{p 0 4}Romano, Joseph P. and Michael Wolf (2005b). Stepwise Multiple Testing as Formalized Data Snooping. {it: Econometrica} 73-4,
https://doi.org/10.1111/j.1468-0262.2005.00615.x.

{p 0 4}Romano, Joseph P. and Michael Wolf (2016). Efficient computation of adjusted p-values for resampling-based stepdown
multiple testing. {it: Statistics and Probability Letters} 113, https://doi.org/10.1016/j.spl.2016.02.012.
