# Analysis of Conservative Effects (ACE)
## Research repository for the most and least differnece in means

Note: run "set_wd_here.R" in base directoy first to set working directory to the base folder for the repository. All code assumes that is the location for the workign directory.

Folder Structure
1. __R/__: general r code for calculating mdm, ldm, agreement contests, row effect size metrics etc.
   1. __aces.R__: functions to calculate most difference in means and least difference in means statistics
   2. __agreement_contests.R__: simulation toolbox to test candidate measures in identifying experiments with higher agreement
   3. __coverage_error_toolbox.R__: simulation toolbox to test coverage error of mdm and rmdm (coverage error is the complement of coverage probability used for confidence intervals, code tests how often the mdm and rmdm is wrongly less than the pop. difference in means and rel. differnece in means).
   4. __credibility_rate_toolbox.R__: simulation toolbox to test coverage error of mdm and rmdm (coverage error is the complement of coverage probability used for confidence intervals, code tests how often the mdm and rmdm is wrongly less than the pop. difference in means and rel. differnece in means).
   5. __parallel_utils.R__: helper functions for parallel processing.
   6. __row_stats_toolbox.R__: helper functions for parallel processing.
3. __mdm_t/__: scripts used to generate figures for the most difference 
   1. __figure1_intro_example.R__: A series of examples to illustrate the three parameters of unscaled agreement and relative agreement. Each example compares the difference in means or relative difference in means between a hypothetical control and experiment group.
   2. __figure2_covary_with_disagreement_SF2-5.R__: Sample a series of population parameter configurations of a control and experiment group where one parameter of agreement is swept towards increasing disagreement. Candidate statistics means are calculated based on repeated samples. Correlation between mean value of each statistic and the agreement parameters are calculated and visualized in a heat map table. 
   3. __figure2_illustrative_disagreement_examples.R__: A series of examples to illustrate the three parameters of unscaled agreement and relative agreement. Each example compares the difference in means or relative difference in means between a hypothetical control and experiment group.
   4. __figure3_risk_assessment_summary.R__: Runs a collection of risk assessments for test how often each candidate statistic incorrectly predicts which of two results have greater noteworthiness (higher evidence of practical significance)
   5. __figure4-5_calculate_stats_from_sxls.R__: creates excel table of results used as applied examples in paper
   7. __sfigure1_mdm_credibility_rate.R__: calculates credibility rate of mdm with Monte Carlo method
   8. __sfigure1_rmdm_credibility_rate.R__: calculates credibility rate of mdm with Monte Carlo method.
   9. __sfigure6-7_unscaled_risk_assessment_null.r__: Risk assessment for mu_dm, sigma_dm, df_dm, alpha_dm with null results.
   10. __sfigure8-9_unscaled_risk_assessment_critical.r__: Risk assessment for mu_dm, sigma_dm, df_dm, alpha_dm with positive results.
   11. __sfigure10-11_relative_risk_assessment_null.R__: Risk assessment for rmu_dm, rsigma_dm, df_dm, alpha_dm with null results.
   12. __sfigure12-13_risk_assessment_contest_crit.r__: Risk assessment for rmu_dm, rsigma_dm, df_dm, alpha_dm with positive results.

