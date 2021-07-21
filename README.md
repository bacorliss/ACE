# Analysis of Conservative Effects (ACE)
## Research repository for the most and least differnece in means

Note: run "set_wd_here.R" in base directoy first to set working directory to the base folder for the repository. All code assumes that is the location for the workign directory.

Folder Structure
1. __R/__: general r code for calculating mdm, ldm, agreement contests, row effect size metrics etc.
   1. __mdm.R__: toolbox for calculating most difference in means and relative most difference in means.
   2. __ldm.R__: toolbox for calculating least difference in means and relative least difference in means.
   3. __rationormal_toolbox.R__: toolbox for calculating confidence intervals of the ratio of two normal variables.
   4. __agreement_contests.R__: simulation toolbox to test candidate measures in identifying experiments with higher agreement
   5. __coverage_error_toolbox.R__: simulation toolbox to test coverage error of mdm and rmdm (coverage error is the complement of coverage probability used for confidence intervals, code tests how often the mdm and rmdm is wrongly less than the pop. difference in means and rel. differnece in means).
   6. __row_stats_toolbox.R__: toolbox for simulation testing that adds supports for candidate effect size measures to evaluate a large number of samples (row based operations).
   7. __parallel_utils.R__: helper functions for parallel processing.
3. __mdm_z/__: scripts used to generate figures for the most difference in means z distribution paper (note that figure numbers are not currently in sync with paper)
   1. __figure1_consistency_with_similarity.r__:
   2. __figure1_illustrative_examples.R__:
   3. __figure2_agreement_contests_summary.R__:
   4. __figure3_applied_examples_tot_cholesterol.R__:
   5. __figure3_calculate_stats_from_sxls.R__:
   6. __figure4_applied_examples_plaque_area.R__:
   7. __sfigure1_fold_normal_dependence_of_mu_sigma.r__:
   8. __sfigure1_null_es_follows_folded_normal.R__:
   9. __sfigure2_integration_over_folded_normal.r__:
   10. __sfigure3_mdm_compared_to_rTOST.r__:
   11. __sfigure3_mdm_relation_confidence_intervals.R__:
   12. __sfigure4_mdm_error_rate_2d.r__:
   13. __sfigure5_rmdm_error_rate_2d.r__:
   14. __sfigure8_unscaled_agreement_contest_null.R__:
   15. __sfigure9_unscaled_agreement_contest_critical.R__:
   16. __sfigure10_relative_agreement_contest_null.R__:
   17. __sfigure11_relative_agreement_contest_crit.R__:
   18. __sfigure12_power_contest.R__:
   19. __stable1_area_folded_normals.R__:
   20. __stable2_mmd95_coeff_scaled_by_CL95_CL90.R__:

