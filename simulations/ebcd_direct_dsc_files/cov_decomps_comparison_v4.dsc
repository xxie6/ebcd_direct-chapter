# simulate samples from 4 population balanced tree
baltree_4pop: 4pop_balancedtree_sim.R
    # module input and variables
    pop_sizes: rep(40, 4), c(20,50,30,60)
    n_genes: 1000
    branch_sds: rep(2, 7)
    indiv_sd: 1
    seed: R{1:10}
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

# simulate samples from 10 group overlapping binary structure
group_overlap: group_overlapping_sim.R
    # module input and variables
    n: 100
    p: 1000
    k: 10
    group_sd: 1 
    indiv_sd: 1
    seed: R{1:10}
    pi0: 0.1, 0.3
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

# simulate samples from non-overlapping binary structure
group_nonoverlap: group_nonoverlapping_sim.R
    # module input and variables
    pop_sizes: rep(40, 4), c(20,50,30,60)
    n_genes: 1000
    branch_sds: rep(2, 4)
    indiv_sd: 1
    seed: R{1:10}
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

# simulate samples from 4 population unbalanced tree
unbaltree_4pop: 4pop_unbalancedtree_sim.R
    # module input and variables
    pop_sizes: rep(40, 4)
    n_genes: 1000
    branch_sds: rep(2, 7)
    indiv_sd: 1
    seed: R{1:10}
    # module decoration
    @ALIAS: args = List()
    # module output
    $data: data
    $true_L: data$LL
    $dat: data$YYt

pca: runpca.R
    # module input and variables
    Kmax_factor: 1, 2
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: pca_data$pca_fit
    $est_L: pca_data$est_L
    $est_LLt: pca_data$est_LLt

ebcd: runebcd.R
    # module input and variables
    seed: 1
    ebnm_fn: ebnm::ebnm_generalized_binary, ebnm::ebnm_point_exponential
    Kmax_factor: 1, 2
    # module decoration
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: ebcd_data$ebcd_fit
    $est_L: ebcd_data$ebcd_fit$EL
    $est_LLt: ebcd_data$est_LLt

ebmfcov_diag: runebmfcov_diag.R
    # module input and variables
    ebnm_fn: ebnm::ebnm_generalized_binary, ebnm::ebnm_point_exponential
    Kmax_factor: 1, 2
    backfit: TRUE, FALSE
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: ebmfcov_diag_data$ebmfcov_diag_fit
    $est_L: ebmfcov_diag_data$ebmfcov_diag_fit$fl$L_pm
    $est_LLt: ebmfcov_diag_data$est_LLt

gbcd: rungbcd_L_only.R
    # module input and variables
    ebnm_fn: ebnm::ebnm_generalized_binary, ebnm::ebnm_point_exponential
    Kmax_factor: 1, 2
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: gbcd_data$fit_obj
    $est_L: gbcd_data$fit_obj$L
    $est_LLt: gbcd_data$est_LLt

sindclus: runsindclus.R
    # module input and variables
    seed: 1
    init_method: "default"
    additive_term: TRUE, FALSE
    off_diagonal: FALSE
    Kmax_factor: 1, 2
    # module decoration
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: sindclus_data$sindclus_fit
    $est_L: sindclus_data$sindclus_fit$P
    $est_LLt: sindclus_data$est_LLt

sympres: runsympres.R
    # module input and variables
    seed: 1
    init_method: "default"
    additive_term: TRUE, FALSE
    off_diagonal: FALSE
    Kmax_factor: 1, 2
    # module decoration
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: sympres_data$sympres_fit
    $est_L: sympres_data$sympres_fit$P
    $est_LLt: sympres_data$est_LLt

codesymnmf: runcodesymnmf.R
    # module input and variables
    Kmax_factor: 1, 2
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: codesymnmf_data$codesymnmf_fit
    $est_L: codesymnmf_data$codesymnmf_fit$H
    $est_LLt: codesymnmf_data$est_LLt

flash_seminmf: runflash_seminmf.R
    # module input and variables
    ebnm_fn: ebnm::ebnm_generalized_binary
    Kmax_factor: 1, 2
    backfit: TRUE, FALSE
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: flash_seminmf_data$flash_fit
    $est_L: flash_seminmf_data$flash_fit$L_pm
    $est_LLt: flash_seminmf_data$est_LLt

laplace_split_ebmfcov_diag: runlaplace_split_ebmfcov_diag.R
    # module input and variables
    ebnm_fn: ebnm::ebnm_generalized_binary, ebnm::ebnm_point_exponential
    Kmax_factor: 1, 2
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: ls_ebmfcov_diag_data$ls_ebmfcov_diag_fit
    $est_L: ls_ebmfcov_diag_data$ls_ebmfcov_diag_fit$fl$L_pm
    $est_LLt: ls_ebmfcov_diag_data$est_LLt

symebcovmf: runsymebcovmf.R
    # module input and variables
    ebnm_fn: ebnm::ebnm_generalized_binary, ebnm::ebnm_point_exponential
    Kmax_factor: 1, 2
    backfit: TRUE, FALSE
    @ALIAS: args = List()
    input: $data
    # module output
    $fit_obj: symebcovmf_data$symebcovmf_fit
    $est_L: symebcovmf_data$symebcovmf_fit$L_pm
    $est_LLt: symebcovmf_data$est_LLt

#odsymnmf: runODsymnmf.R
    # module input and variables
    #method_type: 2
    # module decoration
    #@ALIAS: args = List()
    #input: $data
    # module output
    #$fit_obj: odsymnmf_data$odsymnmf_fit
    #$est_L: odsymnmf_data$odsymnmf_fit
    #$est_LLt: odsymnmf_data$est_LLt

#symnmf_anls: runsymnmfanls.R
    # module input and variables
    #input: $data
    # module output
    #$filename: symnmfanls_data$filename
    #$fit_obj: symnmfanls_data$symnmfanls_fit
    #$est_L: symnmfanls_data$symnmfanls_fit
    #$est_LLt: symnmfanls_data$est_LLt
    
#symnmf_newton: runsymnmfnewton.R
    # module input and variables
    #input: $data
    # module output
    #$filename: symnmfnewton_data$filename
    #$fit_obj: symnmfnewton_data$symnmfnewton_fit
    #$est_L: symnmfnewton_data$symnmfnewton_fit
    #$est_LLt: symnmfnewton_data$est_LLt

crossprod_similarity: score_crossprod_sim2.R
    # module input and variables
    est: $est_L
    truth: $true_L
    # module output
    $result: result

cov_L2_fit: score_L2_fit.R
    # module input and variables
    est: $est_LLt
    dat: $dat
    # module output
    $result: result

prop_est_high_cos_sim: prop_est_high_cos_sim.R
    # module input and variables
    threshold: 0.9, 0.99
    @ALIAS: args = List()
    est: $est_L
    truth: $true_L
    # module output
    $result: result

prop_true_high_cos_sim: prop_true_high_cos_sim.R
   # module input and variables
   threshold: 0.9, 0.99
   @ALIAS: args = List()
   est: $est_L
   truth: $true_L
   # module output
   $result: result 

diff_num_factors: compute_diff_num_factors.R
   # module input and variables
   est: $est_L
   truth: $true_L
   # module output
   $result: result

symebcovmf_fail_check: symebcovmf_fail_check.R
   # module input and variables
   symebcovmf_fit: $fit_obj
   # module output
   $result: result

DSC:
    # module ensembles
    define:
      simulate: baltree_4pop, group_overlap, group_nonoverlap, unbaltree_4pop
      analyze: pca, ebmfcov_diag, gbcd, ebcd, sindclus, sympres, codesymnmf, laplace_split_ebmfcov_diag, flash_seminmf, symebcovmf
      #analyze2:
      score: crossprod_similarity, cov_L2_fit, prop_est_high_cos_sim, prop_true_high_cos_sim, diff_num_factors
    # pipelines
    run:
      standard: simulate * analyze * score
      symebcovmf_fail_yn: simulate * symebcovmf * symebcovmf_fail_check
    output: cov_decomps_comparison_v4

