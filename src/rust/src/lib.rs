#![allow(non_snake_case)]

// Note: extendr only export the documentation above impl blocks
// It's then necessary to repeat the struct documentation above impl blocks
// to have them appear in the R documentation.

///////////////////////////////////////////////////////////////
/// Main dependencies
///////////////////////////////////////////////////////////////

use extendr_api::prelude::*;
use extendr_api::wrapper::symbol::{names_symbol, row_names_symbol};
use gpredomics::individual::AdditionalMetrics;
use gpredomics::individual::{RAW_TYPE, PREVALENCE_TYPE, LOG_TYPE};
use gpredomics::individual::{BINARY_LANG, TERNARY_LANG, POW2_LANG, RATIO_LANG, MCMC_GENERIC_LANG};
use gpredomics::param::Param as GParam;
use gpredomics::data::Data as GData;
use gpredomics::data::{FeatureAnnotations, SampleAnnotations};
use gpredomics::param::get as GParam_get;
use gpredomics::population::Population as GPopulation;
use gpredomics::individual::Individual as GIndividual;
use gpredomics::experiment::Experiment as GExperiment;
use gpredomics::experiment::ExperimentMetadata;
use gpredomics::experiment::{ImportanceAggregation, ImportanceCollection, ImportanceScope, ImportanceType};
use gpredomics::voting::Jury as GJury;
use gpredomics::{ run, run_on_data };
use log::warn;

use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use flexi_logger::{Logger, LoggerHandle, WriteMode, FileSpec, LogSpecification};
use chrono::Local;
use std::collections::HashMap;
use std::collections::HashSet;
use once_cell::sync::OnceCell;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::{RngCore, SeedableRng};
use rayon::iter::IntoParallelRefIterator;   
use rayon::ThreadPoolBuilder;
use rayon::iter::ParallelIterator;

#[inline]
fn r_message(msg: impl AsRef<str>) {
    let safe_msg = msg.as_ref().replace("\"", "\\\"");
    let _ = extendr_api::eval_string(&format!("message(\"gpredomicsR: {}\")", safe_msg));
}

#[inline]
fn r_warning(msg: impl AsRef<str>) {
    let safe_msg = msg.as_ref().replace("\"", "\\\"");
    let _ = extendr_api::eval_string(&format!("warning(\"gpredomicsR: {}\")", safe_msg));
}

#[inline]
fn r_error(msg: impl AsRef<str>) -> Error {
    let safe_msg = msg.as_ref().replace("\"", "\\\"");
    let _ = extendr_api::eval_string(&format!("stop(\"gpredomicsR: {}\")", safe_msg));
    Error::from(format!("gpredomicsR: {}", msg.as_ref()))
}

#[inline]
fn r_print_error(msg: impl AsRef<str>) {
    let safe_msg = msg.as_ref().replace("\"", "\\\"");
    let _ = extendr_api::eval_string(&format!("stop(\"gpredomicsR: {}\")", safe_msg));
}

/// Reconstruct a Data subset from a full dataset using sample IDs.
/// 
/// This function is used to extract cross-validation fold data by matching
/// sample IDs from the fold specification to the full training dataset.
/// 
/// # Arguments
/// * `full_data` - The complete dataset containing all samples
/// * `sample_ids` - Vector of sample IDs to extract
/// 
/// # Returns
/// A new GData object containing only the specified samples
fn reconstruct_fold_data(full_data: &GData, sample_ids: &[String]) -> GData {
    let indices: Vec<usize> = sample_ids
        .iter()
        .filter_map(|id| {
            full_data.samples.iter().position(|s| s == id)
        })
        .collect();
    
    full_data.subset(indices)
}

///////////////////////////////////////////////////////////////
/// Running Flag
///////////////////////////////////////////////////////////////

#[derive(Debug, Clone)]
#[extendr]
pub struct RunningFlag {
    flag: Arc<AtomicBool>,
}

/// @title RunningFlag
/// @name RunningFlag
/// @description A struct to manage the running flag for controlling algorithm execution.
/// 
/// @details
/// The RunningFlag object allows you to control the execution of long-running algorithms.
/// It can be used to stop algorithms gracefully from another thread or after user interruption.
/// 
/// @section Methods:
/// \describe{
///   \item{\code{new()}}{Create a new RunningFlag (initially set to TRUE).}
///   \item{\code{stop()}}{Set the running flag to FALSE to signal algorithm termination.}
///   \item{\code{is_running()}}{Check the current value of the running flag (TRUE or FALSE).}
///   \item{\code{reset()}}{Reset the running flag to TRUE.}
/// }
/// 
/// @export
#[extendr]
impl RunningFlag {
    pub fn new() -> Self {
        Self {
            flag: Arc::new(AtomicBool::new(true)),
        }
    }

    pub fn stop(&self) {
        self.flag.store(false, Ordering::Relaxed);
    }

    pub fn is_running(&self) -> bool {
        self.flag.load(Ordering::Relaxed)
    }

    pub fn reset(&self) {
        self.flag.store(true, Ordering::Relaxed);
    }

}

impl RunningFlag {
    /// Get a clone of the inner Arc (useful for passing to long-running functions).
    fn get_arc(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.flag)
    }
}

///////////////////////////////////////////////////////////////
/// Param object
///////////////////////////////////////////////////////////////


#[extendr]
#[derive(Clone)]
pub struct Param {
    intern: GParam
}

/// @title Param
/// @name Param
/// @description Gpredomics parameter object that stores all algorithm settings.
/// 
/// @details
/// The Param object contains all configuration settings for running gpredomics algorithms
/// including genetic algorithm parameters, data parameters, cross-validation settings, etc.
/// 
/// @section Methods:
/// \describe{
///   \item{\code{new()}}{Create a new empty Param object.}
///   \item{\code{load(file_path)}}{Load a param.yaml file to create a new Param. 
///     \itemize{
///       \item \code{file_path}: Path to param.yaml file
///     }
///   }
///   \item{\code{get()}}{Returns an R list representing the current state of Param with all settings.}
///   \item{\code{set(variable, value)}}{Set a numeric parameter by name.
///     \itemize{
///       \item \code{variable}: Name of the parameter to set
///       \item \code{value}: New numeric value for the parameter
///     }
///   }
///   \item{\code{set_string(variable, string)}}{Set a string parameter by name.
///     \itemize{
///       \item \code{variable}: Name of the parameter to set
///       \item \code{string}: New string value for the parameter
///     }
///   }
///   \item{\code{set_bool(variable, value)}}{Set a boolean parameter by name.
///     \itemize{
///       \item \code{variable}: Name of the parameter to set
///       \item \code{value}: New boolean value for the parameter
///     }
///   }
///   \item{\code{address()}}{Get memory address of this Param object as a string.}
/// }
/// 
/// @export
#[extendr]
impl Param {
    pub fn new() -> Self {
        Self {
            intern: GParam::new()
        }
    }

    /// @description Load a param.yaml file to create a new Param.
    /// @param file_path Path to param.yaml file
    /// @return A new Param object
    pub fn load(file_path: String) -> Result<Self> {
        if let Ok(this_param) = GParam_get(file_path.clone()) {
            Ok(Self {
                intern: this_param
            })
        }
        else {
            r_error(&format!("Unsuitable param file: {}",file_path));
            Err(format!("Unsuitable param file: {}",file_path))?
        }
    } 


    /// @title Get Parameter description
    /// @name Param$get
    /// @description Returns an R list representing the current state of Param with all settings.
    /// @return R object containing parameter details
    pub fn get(&self) -> Robj {
        // Convert General fields
        let general = List::from_pairs(vec![
            ("seed", Robj::from(self.intern.general.seed)),
            ("algo", Robj::from(self.intern.general.algo.clone())),
            ("language", Robj::from(self.intern.general.language.clone())),
            ("data_type", Robj::from(self.intern.general.data_type.clone())),
            ("data_type_epsilon", Robj::from(self.intern.general.data_type_epsilon)),
            ("thread_number", Robj::from(self.intern.general.thread_number)),
            ("log_base", Robj::from(self.intern.general.log_base.clone())),
            ("log_suffix", Robj::from(self.intern.general.log_suffix.clone())),
            ("log_level", Robj::from(self.intern.general.log_level.clone())),
            ("fit", Robj::from(format!("{:?}",self.intern.general.fit))),
            ("k_penalty", Robj::from(self.intern.general.k_penalty)),
            ("fr_penalty", Robj::from(self.intern.general.fr_penalty)),
            ("bias_penalty", Robj::from(self.intern.general.bias_penalty)),
            ("threshold_ci_penalty", Robj::from(self.intern.general.threshold_ci_penalty)),
            ("threshold_ci_alpha", Robj::from(self.intern.general.threshold_ci_alpha)),
            ("threshold_ci_n_bootstrap", Robj::from(self.intern.general.threshold_ci_n_bootstrap)),
            ("threshold_ci_frac_bootstrap", Robj::from(self.intern.general.threshold_ci_frac_bootstrap)),
            ("user_penalties_weight", Robj::from(self.intern.general.user_penalties_weight)),
            ("n_model_to_display", Robj::from(self.intern.general.n_model_to_display)),
            ("gpu", Robj::from(self.intern.general.gpu)),
            ("cv", Robj::from(self.intern.general.cv)),
            ("display_colorful", Robj::from(self.intern.general.display_colorful)),
            ("keep_trace", Robj::from(self.intern.general.keep_trace)),
            ("save_exp", Robj::from(self.intern.general.save_exp.clone()))
        ]);

        // Convert Data fields
        let data = List::from_pairs(vec![
            ("X", Robj::from(self.intern.data.X.clone())),
            ("y", Robj::from(self.intern.data.y.clone())),
            ("Xtest", Robj::from(self.intern.data.Xtest.clone())),
            ("ytest", Robj::from(self.intern.data.ytest.clone())),
            ("holdout_ratio", Robj::from(self.intern.data.holdout_ratio)),
            ("feature_annotations", Robj::from(self.intern.data.feature_annotations.clone())),
            ("sample_annotations", Robj::from(self.intern.data.sample_annotations.clone())),
            ("features_in_rows", Robj::from(self.intern.data.features_in_rows)),
            ("max_features_per_class", Robj::from(self.intern.data.max_features_per_class)),
            ("feature_selection_method", Robj::from(format!("{:?}", self.intern.data.feature_selection_method))),
            ("feature_minimal_prevalence_pct", Robj::from(self.intern.data.feature_minimal_prevalence_pct)),
            ("feature_maximal_adj_pvalue", Robj::from(self.intern.data.feature_maximal_adj_pvalue)),
            ("feature_minimal_feature_value", Robj::from(self.intern.data.feature_minimal_feature_value)),
            ("feature_minimal_log_abs_bayes_factor", Robj::from(self.intern.data.feature_minimal_log_abs_bayes_factor)),
            ("inverse_classes", Robj::from(self.intern.data.inverse_classes)),
            ("n_validation_samples", Robj::from(self.intern.data.n_validation_samples)),
            ("classes", Robj::from(self.intern.data.classes.clone()))
        ]);

        // Convert GA fields
        let ga = List::from_pairs(vec![
            ("population_size", Robj::from(self.intern.ga.population_size)),
            ("max_epochs", Robj::from(self.intern.ga.max_epochs)),
            ("min_epochs", Robj::from(self.intern.ga.min_epochs)),
            ("max_age_best_model", Robj::from(self.intern.ga.max_age_best_model)),
            ("kmin", Robj::from(self.intern.ga.kmin)),
            ("kmax", Robj::from(self.intern.ga.kmax)),
            ("select_elite_pct", Robj::from(self.intern.ga.select_elite_pct)),
            ("select_niche_pct", Robj::from(self.intern.ga.select_niche_pct)),
            ("select_random_pct", Robj::from(self.intern.ga.select_random_pct)),
            ("mutated_children_pct", Robj::from(self.intern.ga.mutated_children_pct)),
            ("mutated_features_pct", Robj::from(self.intern.ga.mutated_features_pct)),
            ("mutation_non_null_chance_pct", Robj::from(self.intern.ga.mutation_non_null_chance_pct)),
            ("forced_diversity_pct", Robj::from(self.intern.ga.forced_diversity_pct)),
            ("forced_diversity_epochs", Robj::from(self.intern.ga.forced_diversity_epochs)),
            ("random_sampling_pct", Robj::from(self.intern.ga.random_sampling_pct)),
            ("random_sampling_epochs", Robj::from(self.intern.ga.random_sampling_epochs))
        ]);

        // Convert GA fields
        let beam = List::from_pairs(vec![
            ("method", Robj::from(format!("{:?}", self.intern.beam.method))),
            ("kmin", Robj::from(self.intern.beam.kmin)),
            ("kmax", Robj::from(self.intern.beam.kmax)),
            ("best_models_ci_alpha", Robj::from(self.intern.beam.best_models_ci_alpha)),            
            ("max_nb_of_models", Robj::from(self.intern.beam.max_nb_of_models)),
        ]);

        // Convert MCMC fields
        let mcmc = List::from_pairs(vec![
            ("n_iter", Robj::from(self.intern.mcmc.n_iter)),
            ("n_burn", Robj::from(self.intern.mcmc.n_burn)),
            ("lambda", Robj::from(self.intern.mcmc.lambda)),
            ("nmin", Robj::from(self.intern.mcmc.nmin)),
            ("save_trace_outdir", Robj::from(self.intern.mcmc.save_trace_outdir.clone()))
        ]);


        // Convert Importance fields
        let importance = List::from_pairs(vec![
            ("compute_importance", Robj::from(self.intern.importance.compute_importance)),
            ("n_permutations_oob",Robj::from(self.intern.importance.n_permutations_oob)),
            ("scaled_importance", Robj::from(self.intern.importance.scaled_importance)),
            ("importance_aggregation", Robj::from(format!("{:?}",self.intern.importance.importance_aggregation)))
        ]);

        // Convert CV fields
        let cv = List::from_pairs(vec![
            ("inner_folds", Robj::from(self.intern.cv.inner_folds)),
            ("overfit_penalty", Robj::from(self.intern.cv.overfit_penalty)),
            ("resampling_inner_folds_epochs", Robj::from(self.intern.cv.resampling_inner_folds_epochs)),
            ("outer_folds", Robj::from(self.intern.cv.outer_folds)),
            ("fit_on_valid", Robj::from(self.intern.cv.fit_on_valid)),
            ("cv_best_models_ci_alpha", Robj::from(self.intern.cv.cv_best_models_ci_alpha)),
            ("stratify_by", Robj::from(self.intern.cv.stratify_by.clone()))
        ]);
        
        // Convert Voting fields
        let voting = List::from_pairs(vec![
            ("vote", Robj::from(self.intern.voting.vote)),
            ("min_perf", Robj::from(self.intern.voting.min_perf)),
            ("min_diversity", Robj::from(self.intern.voting.min_diversity)),
            ("fbm_ci_alpha", Robj::from(self.intern.voting.fbm_ci_alpha)),
            ("method", Robj::from(format!("{:?}",self.intern.voting.method))),
            ("method_threshold", Robj::from(self.intern.voting.method_threshold)),
            ("threshold_windows_pct", Robj::from(self.intern.voting.threshold_windows_pct)),
            ("complete_display", Robj::from(self.intern.voting.complete_display)),
            ("prune_before_voting", Robj::from(self.intern.voting.prune_before_voting)),
        ]);

        // Combine all sections into a single R list object
        let param_list = List::from_pairs(vec![
            ("general", Robj::from(general)),
            ("data", Robj::from(data)),
            ("ga", Robj::from(ga)),
            ("beam", Robj::from(beam)),
            ("mcmc", Robj::from(mcmc)),
            ("importance", Robj::from(importance)),
            ("cv", Robj::from(cv)),
            ("voting", Robj::from(voting)),
        ]);

        Robj::from(param_list)
    }   

    /// @description Set a numeric parameter by name.
    /// @param variable Name of the parameter to set
    /// @param value New numeric value for the parameter
    /// @export 
    pub fn set(&mut self, variable: &str, value: f64) {
        match variable {
            // General parameters
            "seed" => self.intern.general.seed = value as u64,
            "epsilon" | "data_type_epsilon" => self.intern.general.data_type_epsilon = value,
            "thread_number" => self.intern.general.thread_number = value as usize,
            "k_penalty" => self.intern.general.k_penalty = value,
            "fr_penalty" => self.intern.general.fr_penalty = value,
            "bias_penalty" => self.intern.general.bias_penalty = value,
            "threshold_ci_penalty" => self.intern.general.threshold_ci_penalty = value,
            "threshold_ci_alpha" => self.intern.general.threshold_ci_alpha = value,
            "threshold_ci_n_bootstrap" => self.intern.general.threshold_ci_n_bootstrap = value as usize,
            "threshold_ci_frac_bootstrap" => self.intern.general.threshold_ci_frac_bootstrap = value,
            "user_penalties_weight" => self.intern.general.user_penalties_weight = value,
            "n_model_to_display" => self.intern.general.n_model_to_display = value as u32,
            
            // Data parameters
            "holdout_ratio" => self.intern.data.holdout_ratio = value,
            "feature_minimal_prevalence_pct" => self.intern.data.feature_minimal_prevalence_pct = value,
            "feature_maximal_adj_pvalue" => self.intern.data.feature_maximal_adj_pvalue = value,
            "feature_minimal_feature_value" => self.intern.data.feature_minimal_feature_value = value,
            "feature_minimal_log_abs_bayes_factor" => self.intern.data.feature_minimal_log_abs_bayes_factor = value,
            "max_features_per_class" => self.intern.data.max_features_per_class = value as usize,
            "n_validation_samples" => self.intern.data.n_validation_samples = value as usize,
            
            // GA parameters
            "population_size" => self.intern.ga.population_size = value as u32,
            "max_epochs" => self.intern.ga.max_epochs = value as usize,
            "min_epochs" => self.intern.ga.min_epochs = value as usize,
            "max_age_best_model" => self.intern.ga.max_age_best_model = value as usize,
            "k_min" => self.intern.ga.kmin = value as usize,
            "k_max" => self.intern.ga.kmax = value as usize,
            "select_elite_pct" => self.intern.ga.select_elite_pct = value,
            "select_niche_pct" => self.intern.ga.select_niche_pct = value,
            "select_random_pct" => self.intern.ga.select_random_pct = value,
            "mutated_children_pct" => self.intern.ga.mutated_children_pct = value,
            "mutated_features_pct" => self.intern.ga.mutated_features_pct = value,
            "mutation_non_null_chance_pct" => self.intern.ga.mutation_non_null_chance_pct = value,
            "forced_diversity_pct" => self.intern.ga.forced_diversity_pct = value,
            "forced_diversity_epochs" => self.intern.ga.forced_diversity_epochs = value as usize,
            "random_sampling_pct" => self.intern.ga.random_sampling_pct = value,
            "random_sampling_epochs" => self.intern.ga.random_sampling_epochs = value as usize,
            "n_epochs_before_global" => self.intern.ga.n_epochs_before_global = value as usize,
            
            // BEAM parameters
            "best_models_ci_alpha" => self.intern.beam.best_models_ci_alpha = value,
            "max_nb_of_models" => self.intern.beam.max_nb_of_models = value as usize,
            
            // MCMC parameters
            "n_iter" => self.intern.mcmc.n_iter = value as usize,
            "n_burn" => self.intern.mcmc.n_burn = value as usize,
            "lambda" => self.intern.mcmc.lambda = value,
            "n_min" => self.intern.mcmc.nmin = value as u32,
            
            // CV parameters
            "overfit_penalty" => self.intern.cv.overfit_penalty = value,
            "inner_folds" => self.intern.cv.inner_folds = value as usize,
            "resampling_inner_folds_epochs" => self.intern.cv.resampling_inner_folds_epochs = value as usize,
            "outer_folds" => self.intern.cv.outer_folds = value as usize,
            "cv_best_models_ci_alpha" => self.intern.cv.cv_best_models_ci_alpha = value,
            
            // Voting parameters
            "min_perf" => self.intern.voting.min_perf = value,
            "min_diversity" => self.intern.voting.min_diversity = value,
            "fbm_ci_alpha" => self.intern.voting.fbm_ci_alpha = value,
            "method_threshold" => self.intern.voting.method_threshold = value,
            "threshold_windows_pct" => self.intern.voting.threshold_windows_pct = value,
            
            // Importance parameters
            "n_permutations_oob" => self.intern.importance.n_permutations_oob = value as usize,
            
            // GPU parameters
            "max_total_memory_mb" => self.intern.gpu.max_total_memory_mb = value as u64,
            "max_buffer_size_mb" => self.intern.gpu.max_buffer_size_mb = value as u32,

            "X"|"y"|"Xtest"|"ytest"|"feature_selection_method"|"algo"|"language"|"data_type"|"fit" => r_warning(format!("Use param$set_string() for {}",variable)),
            "keep_trace"|"gpu"|"cv"|"display_colorful"|"features_in_rows"|"inverse_classes"|"fit_on_valid"|"vote"|"complete_display"|"prune_before_voting"|"compute_importance"|"scaled_importance"|"fallback_to_cpu" => r_warning(format!("Use param$set_bool() for {}",variable)),
            "log_level"|"log_base"|"log_suffix" => r_print_error("Cannot set logs this way, create or get back your GLogger object"),
            _ => r_print_error(format!("Unknown variable: {} ", variable))
        }
    }

    /// @title Get memory address
    /// @name Param$address
    /// @description Get memory address of this Param object
    /// @return String representing the memory address
    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern)
    }

    /// @description Set a string parameter by name.
    /// @param variable Name of the parameter to set
    /// @param string New string value for the parameter
    /// @export
    pub fn set_string(&mut self, variable: &str, string: String) {
        match variable {
            "stratify_by" => self.intern.cv.stratify_by = string,

            // Data parameters
            "X" => self.intern.data.X = string,
            "y" => self.intern.data.y = string,
            "Xtest" => self.intern.data.Xtest = string,
            "ytest" => self.intern.data.ytest = string,
             "feature_annotations" => self.intern.data.feature_annotations = string,
            "sample_annotations" => self.intern.data.sample_annotations = string,
            "feature_selection_method" => self.intern.data.feature_selection_method = match string.to_lowercase().as_str() {
                "wilcoxon" => gpredomics::data::PreselectionMethod::wilcoxon,
                "studentt" => gpredomics::data::PreselectionMethod::studentt,
                "bayesian_fisher" => gpredomics::data::PreselectionMethod::bayesian_fisher,
                _ => { r_warning(format!("Unknown feature selection method: {}. Currently only 'wilcoxon', 'studentt', and 'bayesian_fisher' are supported.", string)); 
                self.intern.data.feature_selection_method.clone()}
            },
            
            // General parameters
            "algo" => self.intern.general.algo = string,
            "language" => self.intern.general.language = string,
            "data_type" => self.intern.general.data_type = string,
            
            // BEAM parameters
            "method" => self.intern.beam.method = match string.to_lowercase().as_str() { 
                "LimitedExhaustive" => gpredomics::beam::BeamMethod::LimitedExhaustive,
                "ParallelForward" => gpredomics::beam::BeamMethod::ParallelForward,
                _ => {r_print_error(format!("Unknown BEAM method: {}. Currently only 'LimitedExhaustive' and 'ParallelForward' are supported.", string)); 
                self.intern.beam.method.clone()}
                },
            
            // MCMC parameters
            "fit" => self.intern.general.fit = match string.to_lowercase().as_str() { 
                "auc" => gpredomics::param::FitFunction::auc,
                "specificity" => gpredomics::param::FitFunction::specificity,
                "sensitivity" => gpredomics::param::FitFunction::sensitivity,
                "mcc" =>  gpredomics::param::FitFunction::mcc,
                "npv" => gpredomics::param::FitFunction::npv,
                "ppv" => gpredomics::param::FitFunction::ppv,
                "g_mean" => gpredomics::param::FitFunction::g_mean,
                "f1_score" => gpredomics::param::FitFunction::f1_score,
                _ => { r_print_error(format!("Unknown fit function: {}. Currently supported: auc, specificity, sensitivity, mcc, npv, ppv, g_mean, f1_score.", string)); 
                self.intern.general.fit.clone() }
             }, //self.intern.general.fit = string,
            "log_level"|"log_base"|"log_suffix" => r_print_error("Cannot set logs this way, create or get back your GLogger object"),
            _ => r_print_error(format!("Variable unknown or not settable byt set_string: {} ", variable))
        }
    }

    /// @description Set a boolean parameter by name.
    /// @param variable Name of the parameter to set
    /// @param value New boolean value for the parameter
    /// @export
    pub fn set_bool(&mut self, variable: &str, value: bool) {
        match variable {
            // General parameters
            "gpu" => self.intern.general.gpu = value,
            "cv" => self.intern.general.cv = value,
            "display_colorful" => self.intern.general.display_colorful = value,
            "keep_trace" => { warn!("keep_trace modified. Please note that data from previous generations is only accessible when keep_trace is set to true. \
            Some future analyses may be unavailable if keep_trace is set to false.");
             self.intern.general.keep_trace = value }
            
            // Data parameters
            "features_in_rows" => self.intern.data.features_in_rows = value,
            "inverse_classes" => self.intern.data.inverse_classes = value,
           
            
            // CV parameters
            "fit_on_valid" => self.intern.cv.fit_on_valid = value,
            
            // Voting parameters
            "vote" => self.intern.voting.vote = value,
            "complete_display" => self.intern.voting.complete_display = value,
            "prune_before_voting" => self.intern.voting.prune_before_voting = value,
            
            // Importance parameters
            "compute_importance" => self.intern.importance.compute_importance = value,
            "scaled_importance" => self.intern.importance.scaled_importance = value,
            
            // GPU parameters
            "fallback_to_cpu" => self.intern.gpu.fallback_to_cpu = value,
            
            "log_level"|"log_base"|"log_suffix" => r_print_error("Cannot set logs this way, create or get back your GLogger object"),
            _ => r_print_error(format!("Unknown boolean variable: {}", variable)),
            }
    }

}

///////////////////////////////////////////////////////////////
/// Data object
///////////////////////////////////////////////////////////////

/// Convert a sparse matrix (HashMap) to an R data.frame.
/// 
/// This internal helper function transforms a sparse matrix representation into an R data.frame.
/// Missing entries in the sparse matrix are filled with 0.0.
/// 
/// # Arguments
/// * `matrix` - Sparse matrix as HashMap with keys (col_index, row_index) and f64 values
/// * `column_names` - Names for each column (length must equal number of columns)
/// * `row_names` - Names for each row (length must equal number of rows)
/// 
/// # Returns
/// An R data.frame object with the specified dimensions and names
fn sparse_matrix_to_dataframe(
    matrix: &HashMap<(usize, usize), f64>,
    column_names: &Vec<String>,
    row_names: &Vec<String>,
) -> Robj {

    let mut columns: Vec<Vec<f64>> = vec![ vec![0.0; row_names.len()]; column_names.len() ];

    for ((col,row), value) in matrix {
        columns[*col][*row] = *value;
    }

    let mut df_list = List::from_values(columns.iter().map(|c| {c.into_robj()}));
    let _ = df_list.set_names(column_names);

    df_list.set_attrib(row_names_symbol(), Robj::from(row_names)).ok();

    df_list.set_class(&["data.frame"]).ok();

    df_list.into()

}

fn feature_tags_to_dataframe_named(
    tags: &HashMap<usize, Vec<String>>,
    features: &[String],
    tag_column_names: &[String],
) -> Robj {
    if tags.is_empty() {
        return Robj::from(());
    }
    
    let max_tag_cols = tags.values()
        .map(|v| v.len())
        .max()
        .unwrap_or(0);
    
    if max_tag_cols == 0 {
        return Robj::from(());
    }
    
    let mut feature_names = Vec::new();
    let mut tag_columns: Vec<Vec<String>> = vec![Vec::new(); max_tag_cols];
    
    let mut sorted_indices: Vec<&usize> = tags.keys().collect();
    sorted_indices.sort();
    
    for &fidx in sorted_indices.clone() {
        let fname = if fidx < features.len() {
            features[fidx].clone()
        } else {
            format!("Feature_{}", fidx)
        };
        feature_names.push(fname);
        
        let tag_vec = &tags[&fidx];
        for (col_idx, tag_col_vec) in tag_columns.iter_mut().enumerate() {
            if col_idx < tag_vec.len() {
                tag_col_vec.push(tag_vec[col_idx].clone());
            } else {
                tag_col_vec.push(String::new());
            }
        }
    }

    let n_cols = 1 + max_tag_cols;
    let mut df = List::new(n_cols);
    let mut col_names = Vec::with_capacity(n_cols);
    
    df.set_elt(0, Robj::from(feature_names)).ok();
    col_names.push("feature".to_string());
    
    for (i, col) in tag_columns.into_iter().enumerate() {
        let col_name = if i < tag_column_names.len() {
            tag_column_names[i].clone()
        } else {
            format!("tag_{}", i + 1)
        };
        
        df.set_elt(i + 1, Robj::from(col)).ok();
        col_names.push(col_name);
    }
    
    df.set_names(col_names).ok();
    
    let rownames: Vec<String> = (1..=sorted_indices.len())
        .map(|i| i.to_string())
        .collect();
    df.set_attrib(row_names_symbol(), rownames).ok();
    df.set_class(&["data.frame"]).ok();
    
    Robj::from(df)
}


/// Convert R named vectors and data.frame to FeatureAnnotations
fn process_feature_annotations(
    features: &[String],
    prior_weight_robj: Option<Robj>,
    feature_penalty_robj: Option<Robj>,
    feature_tags_robj: Option<Robj>,
) -> Result<FeatureAnnotations> {
    
    let mut prior_weight = HashMap::<usize, f64>::new();
    let mut feature_penalty = HashMap::<usize, f64>::new();
    let mut feature_tags = HashMap::<usize, Vec<String>>::new();
    let mut tag_column_names = Vec::<String>::new();
    
    // Process prior_weight: named numeric vector
    if let Some(pw_robj) = prior_weight_robj {
        let pw_vec = pw_robj.as_real_vector()
            .ok_or_else(|| r_error("prior_weight must be a numeric vector"))?;
        
        let pw_names = pw_robj.get_attrib(names_symbol())
            .and_then(|n| n.as_string_vector())
            .ok_or_else(|| r_error("prior_weight must be a named vector"))?;
        
        if pw_vec.len() != pw_names.len() {
            return Err(r_error("prior_weight: length mismatch between values and names"));
        }
        
        for (fname, &weight) in pw_names.iter().zip(pw_vec.iter()) {
            if let Some(fidx) = features.iter().position(|f| f == fname) {
                prior_weight.insert(fidx, weight);
            } else {
                r_warning(format!("prior_weight: feature '{}' not found in data, ignored", fname));
            }
        }
    }
    
    // Process feature_penalty: named numeric vector
    if let Some(fp_robj) = feature_penalty_robj {
        let fp_vec = fp_robj.as_real_vector()
            .ok_or_else(|| r_error("feature_penalty must be a numeric vector"))?;
        
        let fp_names = fp_robj.get_attrib(names_symbol())
            .and_then(|n| n.as_string_vector())
            .ok_or_else(|| r_error("feature_penalty must be a named vector"))?;
        
        if fp_vec.len() != fp_names.len() {
            return Err(r_error("feature_penalty: length mismatch between values and names"));
        }
        
        for (fname, &penalty) in fp_names.iter().zip(fp_vec.iter()) {
            if let Some(fidx) = features.iter().position(|f| f == fname) {
                feature_penalty.insert(fidx, penalty);
            } else {
                r_warning(format!("feature_penalty: feature '{}' not found in data, ignored", fname));
            }
        }
    }
    
    // Process feature_tags: data.frame with first column = feature name
    if let Some(ft_robj) = feature_tags_robj {
        let ft_list = ft_robj.as_list()
            .ok_or_else(|| r_error("feature_tags must be a data.frame"))?;
        
        if ft_list.len() == 0 {
            return Err(r_error("feature_tags is empty"));
        }
        
        // Extract column names
        let col_names = ft_robj.get_attrib(names_symbol())
            .and_then(|n| n.as_string_vector())
            .ok_or_else(|| r_error("feature_tags must have column names"))?;
        
        // First column should be feature identifier
        let first_col = ft_list.values().next().unwrap();
        let feature_names = first_col.as_string_vector()
            .ok_or_else(|| r_error("First column of feature_tags must contain feature names"))?;
        
        let nrows = feature_names.len();
        
        // Store other columns as tags
        tag_column_names = col_names.iter().skip(1).map(|s| s.to_string()).collect();
        
        for row_idx in 0..nrows {
            let fname = &feature_names[row_idx];
            
            if let Some(fidx) = features.iter().position(|f| f == fname) {
                let mut row_tags = Vec::new();
                
                // Extract values from each tag column
                for (col_idx, _col_name) in tag_column_names.iter().enumerate() {
                    let col_robj = ft_list.values().nth(col_idx + 1).unwrap();
                    
                    let tag_value = if let Some(str_vec) = col_robj.as_string_vector() {
                        str_vec.get(row_idx).unwrap_or(&"".to_string()).to_string()
                    } else if let Some(int_vec) = col_robj.as_integer_vector() {
                        int_vec.get(row_idx).map(|i| i.to_string()).unwrap_or_default()
                    } else {
                        String::new()
                    };
                    
                    row_tags.push(tag_value);
                }
                
                if !row_tags.is_empty() {
                    feature_tags.insert(fidx, row_tags);
                }
            } else {
                r_warning(format!("feature_tags: feature '{}' not found in data, ignored", fname));
            }
        }
    }
    
    Ok(FeatureAnnotations {
        tag_column_names,
        feature_tags,
        prior_weight,
        feature_penalty,
    })
}

/// Convert R stratification vector and data.frame to SampleAnnotations
fn process_sample_annotations(
    samples: &[String],
    sample_tags_robj: Option<Robj>,
) -> Result<SampleAnnotations> {
    
    let mut sample_tags = HashMap::<usize, Vec<String>>::new();
    let mut tag_column_names = Vec::<String>::new();
    
    // Process sample_tags: data.frame with first column = sample name
    if let Some(st_robj) = sample_tags_robj {
        let st_list = st_robj.as_list()
            .ok_or_else(|| r_error("sample_tags must be a data.frame"))?;
        
        if st_list.len() == 0 {
            return Err(r_error("sample_tags is empty"));
        }
        
        let col_names = st_robj.get_attrib(names_symbol())
            .and_then(|n| n.as_string_vector())
            .ok_or_else(|| r_error("sample_tags must have column names"))?;
        
        let first_col = st_list.values().next().unwrap();
        let sample_names = first_col.as_string_vector()
            .ok_or_else(|| r_error("First column of sample_tags must contain sample names"))?;
        
        let nrows = sample_names.len();
        
        // Merge with existing tag_column_names or replace
        let new_tag_cols: Vec<String> = col_names.iter().skip(1).map(|s| s.to_string()).collect();
        
        if tag_column_names.is_empty() {
            tag_column_names = new_tag_cols;
        } else {
            // Merge strata with tags: extend each sample's tag vector
            tag_column_names.extend(new_tag_cols);
        }
        
        for row_idx in 0..nrows {
            let sname = &sample_names[row_idx];
            
            if let Some(sidx) = samples.iter().position(|s| s == sname) {
                let mut row_tags = sample_tags.get(&sidx).cloned().unwrap_or_default();
                
                // Extract values from each tag column
                for col_idx in 1..st_list.len() {
                    let col_robj = st_list.values().nth(col_idx).unwrap();
                    
                    let tag_value = if let Some(str_vec) = col_robj.as_string_vector() {
                        str_vec.get(row_idx).unwrap_or(&"".to_string()).to_string()
                    } else if let Some(int_vec) = col_robj.as_integer_vector() {
                        int_vec.get(row_idx).map(|i| i.to_string()).unwrap_or_default()
                    } else {
                        String::new()
                    };
                    
                    row_tags.push(tag_value);
                }
                
                sample_tags.insert(sidx, row_tags);
            } else {
                r_warning(format!("sample_tags: sample '{}' not found in data, ignored", sname));
            }
        }
    }
    
    Ok(SampleAnnotations {
        tag_column_names,
        sample_tags,
    })
}


fn sample_tags_to_dataframe_named(
    tags: &HashMap<usize, Vec<String>>,
    samples: &[String],
    tag_column_names: &[String],
) -> Robj {
    if tags.is_empty() {
        return Robj::from(());
    }
    
    let max_tag_cols = tags.values()
        .map(|v| v.len())
        .max()
        .unwrap_or(0);
    
    if max_tag_cols == 0 {
        return Robj::from(());
    }
    
    let mut sample_names = Vec::new();
    let mut tag_columns: Vec<Vec<String>> = vec![Vec::new(); max_tag_cols];
    
    let mut sorted_indices: Vec<&usize> = tags.keys().collect();
    sorted_indices.sort();
    
    for &sidx in sorted_indices {
        let sname = if sidx < samples.len() {
            samples[sidx].clone()
        } else {
            format!("Sample_{}", sidx)
        };
        
        sample_names.push(sname);
        
        let tag_vec = &tags[&sidx];
        
        for (col_idx, tag_col_vec) in tag_columns.iter_mut().enumerate() {
            if col_idx < tag_vec.len() {
                tag_col_vec.push(tag_vec[col_idx].clone());
            } else {
                tag_col_vec.push(String::new());
            }
        }
    }
    

    let n_cols = 1 + max_tag_cols;
    let mut df = List::new(n_cols);
    let mut col_names = Vec::with_capacity(n_cols);

    df.set_elt(0, Robj::from(sample_names.clone())).ok();
    col_names.push("sample".to_string());
    
    for (i, col) in tag_columns.into_iter().enumerate() {
        let col_name = if i < tag_column_names.len() {
            tag_column_names[i].clone()
        } else {
            format!("tag_{}", i + 1)
        };
        
        df.set_elt(i + 1, Robj::from(col)).ok();
        col_names.push(col_name);
    }
    
    df.set_attrib(names_symbol(), col_names).ok();
    
    
    // Row names
    let rownames: Vec<String> = (1..=sample_names.len())
        .map(|i| i.to_string())
        .collect();
    
    df.set_attrib(row_names_symbol(), rownames).ok();
    // OU: df.set_attrib("row.names", rownames).ok();
    
    df.set_class(&["data.frame"]).ok();
    
    Robj::from(df)
}



/// Convert a sparse vector (HashMap usize f64) to a dense integer vector for R.
/// 
/// This internal helper converts a sparse representation into a dense vector,
/// filling missing indices with zeros.
/// 
/// # Arguments
/// * `vector` - Sparse vector as HashMap with usize keys (indices) and u8 values
/// * `length` - Length of the resulting dense vector
/// 
/// # Returns
/// An R integer vector of the specified length
fn sparse_vector_to_vector_usize_f64(vector: &HashMap<usize, f64>, length: usize) -> Robj {
    // Create a dense vector initialized with zeros
    let mut dense_vector: Vec<i32> = vec![0; length];

    // Populate the dense vector with non-zero values from the sparse vector
    for (index, value) in vector {
        dense_vector[*index] = *value as i32;
    }

    // Convert the dense vector to an R object
    Robj::from(dense_vector)
}

/// Convert a sparse vector (HashMap usize u8) to a dense integer vector for R.
/// 
/// This internal helper converts a sparse representation into a dense vector,
/// filling missing indices with zeros.
/// 
/// # Arguments
/// * `vector` - Sparse vector as HashMap with usize keys (indices) and u8 values
/// * `length` - Length of the resulting dense vector
/// 
/// # Returns
/// An R integer vector of the specified length
fn sparse_vector_to_vector_usize_u8(vector: &HashMap<usize, u8>, length: usize) -> Robj {
    // Create a dense vector initialized with zeros
    let mut dense_vector: Vec<i32> = vec![0; length];

    // Populate the dense vector with non-zero values from the sparse vector
    for (index, value) in vector {
        dense_vector[*index] = *value as i32;
    }

    // Convert the dense vector to an R object
    Robj::from(dense_vector)
}

/// Create a Gpredomics Data object from an R DataFrame `df` (X) and a named label vector `y`.
///
/// Orientation and alignment:
/// - If `features_in_columns = true`: rows are samples, columns are features; `names(y)` must match `rownames(df)`.
/// - If `features_in_columns = false`: columns are samples, rows are features; `names(y)` must match `colnames(df)`.
/// - Samples are strictly aligned by exact string matching; any mismatch of names triggers an error.
/// - Internally, X is stored as a sparse map with keys `(sample_idx, feature_idx)`.
///
/// Label handling (`y`):
/// - Accepts integer, character, or 2-level factor.
/// - If `y` has exactly two classes:
///   - Character or factor: classes are extracted from the actual labels/levels and stored in `Data.classes`.
///   - Integer:
///     - If values are {0,1}, they are kept as-is (classes = ["0","1"]).
///     - If two distinct integers {a,b} are detected, a warning is emitted and a→0, b→1 mapping is applied (classes = [str(a), str(b)]).
/// - If more than two non-missing classes are present, associated samples are classified as unknown.
/// - The binary vector `y` is reordered to match the chosen sample order derived from `df`.
///
/// Requirements:
/// - `df` must be an R data.frame-like object where feature values are numeric/integer/logical.
/// - All sample names must be provided via `rownames(df)` (if `features_in_columns`) or `colnames(df)` (otherwise).
/// - `y` must be named with the exact same sample identifiers used on the chosen axis of `df`.
///
/// Behavior on missing or invalid inputs:
/// - Missing `names(y)`: a warning is emitted; a fallback alignment to the chosen sample order is attempted, but exact names are recommended.
/// - Non-numeric feature columns (character): error advising to encode to numeric first.
/// - Column/row length mismatches: error with the offending index and dimensions.
/// - Any sample without a corresponding label in `y`: error with the sample name.
/// Errors:
/// - Detailed, user-facing messages are emitted via R console and returned as `Err(Error::Other(...))`.
///
/// Notes:
/// - Name normalization is not performed; ensure consistent punctuation (e.g., "-" vs ".") on both `df` and `y`.
/// - For factors, R’s 1-based codes are mapped to their `levels` to retrieve human-readable class names.
/// 
/// @param df R data.frame-like object containing feature values.
/// @param y_vec R vector containing binary labels, named with sample identifiers.
/// @param features_in_columns Boolean indicating if features are in columns (`true`) or rows (`false`).
/// @param prior_weight Named numeric vector: feature name -> prior weight (optional)
/// @param feature_penalty Named numeric vector: feature name -> penalty (optional)
/// @param feature_tags R data.frame with feature annotations (optional)
/// @param sample_tags R data.frame with sample metadata (optional)
/// @return Gpredomics Data object with aligned and processed X and y.
/// @export
#[extendr]
pub fn as_gpredomics_data(
    df: Robj,
    y_vec: Robj,
    features_in_columns: bool,
    prior_weight: Option<Robj>,
    feature_penalty: Option<Robj>,
    feature_tags: Option<Robj>,
    sample_tags: Option<Robj>,
) -> Result<Data> {
    let prior_weight = match prior_weight {
        Some(robj) if !robj.is_null() => Some(robj),
        _ => None,
    };
    let feature_penalty = match feature_penalty {
        Some(robj) if !robj.is_null() => Some(robj),
        _ => None,
    };
    let feature_tags = match feature_tags {
        Some(robj) if !robj.is_null() => Some(robj),
        _ => None,
    };
    let sample_tags = match sample_tags {
        Some(robj) if !robj.is_null() => Some(robj),
        _ => None,
    };

    // DataFrame as list
    let df_list = df.as_list().ok_or_else(|| r_error("df must be a data.frame (list-like)."))?;
    let ncol = df_list.len();

    // First column to infer nrow
    let first_col = df_list.values().next().ok_or_else(|| r_error("data.frame has zero columns."))?;
    let nrow = if let Some(v) = first_col.as_real_vector() {
        v.len()
    } else if let Some(v) = first_col.as_integer_vector() {
        v.len()
    } else if let Some(v) = first_col.as_logical_vector() {
        v.len()
    } else if let Some(v) = first_col.as_string_vector() {
        v.len()
    } else {
        return Err(r_error("Unsupported first column type; expected numeric/integer/logical/character."));
    };

    // Column names (accept string or integer-encoded rownames)
    let colnames: Vec<String> = if let Some(names_obj) = df.get_attrib::<Robj>(names_symbol().into()) {
        if let Some(v) = names_obj.as_string_vector() {
            v.to_vec()
        } else if let Some(v) = names_obj.as_integer_vector() {
            v.iter().map(|i| i.to_string()).collect()
        } else {
            rprintln!("Column names attribute is not string/integer; synthesizing.");
            (0..ncol).map(|i| format!("V{}", i + 1)).collect()
        }
    } else {
        (0..ncol).map(|i| format!("V{}", i + 1)).collect()
    };

    // Row names (accept string or integer-encoded rownames)
    let rownames: Vec<String> = if let Some(rownames_obj) = df.get_attrib::<Robj>(row_names_symbol().into()) {
        if let Some(v) = rownames_obj.as_string_vector() {
            v.to_vec()
        } else if let Some(v) = rownames_obj.as_integer_vector() {
            v.iter().map(|i| i.to_string()).collect()
        } else {
            rprintln!("rownames attribute is not string/integer; synthesizing.");
            (0..nrow).map(|i| format!("row_{}", i + 1)).collect()
        }
    } else {
        (0..nrow).map(|i| format!("row_{}", i + 1)).collect()
    };

    // Orientation
    let (features, samples, feature_len, sample_len) = if features_in_columns {
        (colnames.clone(), rownames.clone(), ncol, nrow)
    } else {
        (rownames.clone(), colnames.clone(), nrow, ncol)
    };

    // y vector + names
    let y_names_raw: Vec<String> = if let Some(names_obj) = y_vec.get_attrib::<Robj>(names_symbol().into()) {
        if let Some(v) = names_obj.as_string_vector()      { v.to_vec() }
        else if let Some(v) = names_obj.as_integer_vector(){ v.iter().map(|i| i.to_string()).collect() }
        else { samples.clone() }
    } else {
        r_warning("y has no names(); set names(y) to chosen sample names before calling.");
        samples.clone()
    };

    // Test for factor FIRST (integer codes + levels attribute)
    // because factors can also be read as string vectors by extendr!
    let (y_labels_u8, classes): (Vec<u8>, Vec<String>) = if let Some(vi) = y_vec.as_integer_vector() {
        // Try factor: integer codes + levels
        let levels: Option<Vec<String>> = y_vec
            .get_attrib::<Robj>(sym!(levels).into())
            .and_then(|lv| lv.as_string_vector().map(|v| v.to_vec()));

        let vals: Vec<i32> = vi.to_vec();
        let mut uniq = vals.clone();
        uniq.sort_unstable();
        uniq.dedup();

        if let Some(lv) = levels {
            // Factor: preserve level order by code (1-based)
            // Allow unknowns (values outside the two main levels) → classified as 2
            let valid_codes: Vec<i32> = uniq.iter().filter(|&&code| code >= 1 && (code as usize) <= lv.len()).copied().collect();
            
            if valid_codes.len() < 2 {
                return Err(r_error(&format!("y factor must have at least two valid levels; found {}.", valid_codes.len())));
            }
            
            // Use first two valid codes only
            let code0 = valid_codes[0];
            let code1 = valid_codes[1];
            let class0 = lv[(code0 - 1) as usize].clone();
            let class1 = lv[(code1 - 1) as usize].clone();
            
            // Map: code0→0, code1→1, others→2 (unknown)
            let mapped: Vec<u8> = vals.iter().map(|&z| {
                if z == code0 { 0 }
                else if z == code1 { 1 }
                else { 2 } // Unknown class
            }).collect();
            
            let unknown_count = mapped.iter().filter(|&&x| x == 2).count();
            if unknown_count > 0 {
                r_warning(&format!("{} samples with values outside the two main factor levels will be classified as 'unknown' (class 2).", unknown_count));
            }
            
            (mapped, vec![class0, class1, "unknown".to_string()])
        } else {
            // Pure integer: {0,1} or two other ints {a,b}, others→2 (unknown)
            match uniq.as_slice() {
                [0] | [1] | [0,1] => {
                    let mapped: Vec<u8> = vals.iter().map(|&z| {
                        if z == 0 { 0 }
                        else if z == 1 { 1 }
                        else { 2 } // Unknown
                    }).collect();
                    let unknown_count = mapped.iter().filter(|&&x| x == 2).count();
                    if unknown_count > 0 {
                        r_warning(&format!("{} samples with values outside {{0,1}} will be classified as 'unknown' (class 2).", unknown_count));
                    }
                    (mapped, vec!["0".into(), "1".into(), "unknown".into()])
                },
                [a, b] => {
                    let (lo, hi) = (*a, *b);
                    r_message(&format!("y has two distinct integer values {{{},{}}}; mapping {}→0, {}→1, others→2.", lo, hi, lo, hi));
                    let mapped: Vec<u8> = vals.iter().map(|&z| {
                        if z == lo { 0 }
                        else if z == hi { 1 }
                        else { 2 } // Unknown
                    }).collect();
                    let unknown_count = mapped.iter().filter(|&&x| x == 2).count();
                    if unknown_count > 0 {
                        r_warning(&format!("{} samples with values outside {{{},{}}} will be classified as 'unknown' (class 2).", unknown_count, lo, hi));
                    }
                    (mapped, vec![lo.to_string(), hi.to_string(), "unknown".into()])
                },
                _ => {
                    // More than 2 distinct values: use first two, rest→unknown
                    r_message(&format!("y has {} distinct integer values; using first two ({}, {}), others→unknown.", uniq.len(), uniq[0], uniq[1]));
                    let (lo, hi) = (uniq[0], uniq[1]);
                    let mapped: Vec<u8> = vals.iter().map(|&z| {
                        if z == lo { 0 }
                        else if z == hi { 1 }
                        else { 2 }
                    }).collect();
                    (mapped, vec![lo.to_string(), hi.to_string(), "unknown".into()])
                }
            }
        }
    } else if let Some(vs) = y_vec.as_string_vector() {
        // Character: 2 distinct labels → lexicographic order, others→unknown
        let vals = vs.to_vec();
        let mut uniq = vals.clone();
        uniq.sort();
        uniq.dedup();
        
        if uniq.len() < 2 {
            return Err(r_error(&format!("y must have at least two string classes; found {}.", uniq.len())));
        }
        
        // Use first two classes alphabetically
        let class0 = uniq[0].clone();
        let class1 = uniq[1].clone();
        
        let mapped: Vec<u8> = vals.iter().map(|s| {
            if *s == class0 { 0 }
            else if *s == class1 { 1 }
            else { 2 } // Unknown
        }).collect();
        
        let unknown_count = mapped.iter().filter(|&&x| x == 2).count();
        if unknown_count > 0 {
            r_warning(&format!("{} samples with values outside {{{}, {}}} will be classified as 'unknown' (class 2).", unknown_count, class0, class1));
        }
        
        (mapped, vec![class0, class1, "unknown".into()])
    } else {
        return Err(r_error("Unsupported y type; provide integer, character, or 2-level factor."));
    };

    if y_names_raw.len() != y_labels_u8.len() {
        return Err(r_error("Length mismatch between names(y) and y values."));
    }

    let y_map: std::collections::HashMap<String, u8> = y_names_raw
        .into_iter()
        .zip(y_labels_u8.into_iter())
        .collect();

    let y_final: Vec<u8> = samples
        .iter()
        .map(|sname| {
            y_map
                .get(sname)
                .copied()
                .ok_or_else(|| r_error(&format!(
                    "No label found for sample '{}'; ensure names(y) match the chosen sample names.",
                    sname
                )))
        })
        .collect::<extendr_api::Result<_>>()?;

    // Build sparse X
    let mut X: HashMap<(usize, usize), f64> = HashMap::new();
    for (col_idx, col_obj) in df_list.values().enumerate() {
        if let Some(v) = col_obj.as_real_vector() {
            if v.len() != nrow {
                return Err(r_error(&format!(
                    "Numeric column length mismatch at {}: got {}, expected {}.",
                    col_idx, v.len(), nrow
                )));
            }
            for (row_idx, &val) in v.iter().enumerate() {
                if val != 0.0 && !val.is_nan() {
                    let (sample_idx, feature_idx) =
                        if features_in_columns { (row_idx, col_idx) } else { (col_idx, row_idx) };
                    X.insert((sample_idx, feature_idx), val);
                }
            }
        } else if let Some(v) = col_obj.as_integer_vector() {
            if v.len() != nrow {
                return Err(r_error(&format!(
                    "Integer column length mismatch at {}: got {}, expected {}.",
                    col_idx, v.len(), nrow
                )));
            }
            for (row_idx, &ival) in v.iter().enumerate() {
                if ival != 0 {
                    let (sample_idx, feature_idx) =
                        if features_in_columns { (row_idx, col_idx) } else { (col_idx, row_idx) };
                    X.insert((sample_idx, feature_idx), ival as f64);
                }
            }
        } else if let Some(v) = col_obj.as_logical_vector() {
            if v.len() != nrow {
                return Err(r_error(&format!(
                    "Logical column length mismatch at {}: got {}, expected {}.",
                    col_idx, v.len(), nrow
                )));
            }
            for (row_idx, lv) in v.iter().enumerate() {
                if lv.is_true() {
                    let (sample_idx, feature_idx) =
                        if features_in_columns { (row_idx, col_idx) } else { (col_idx, row_idx) };
                    X.insert((sample_idx, feature_idx), 1.0);
                }
            }
        } else if col_obj.as_string_vector().is_some() {
            return Err(r_error(&format!(
                "Character column at {} not supported; encode to numeric first.",
                col_idx
            )));
        } else {
            return Err(r_error(&format!(
                "Unsupported column type at index {}.",
                col_idx
            )));
        }
    }

    // Build Data
    let mut gdata = GData::new();
    gdata.X = X;
    gdata.y = y_final;
    gdata.features = features;
    gdata.samples = samples;
    gdata.feature_len = feature_len;
    gdata.sample_len = sample_len;
    gdata.classes = classes;
    gdata.feature_class = HashMap::new();
    gdata.feature_selection = Vec::new();

    if prior_weight.is_some() || feature_penalty.is_some() || feature_tags.is_some() {
        gdata.feature_annotations = Some(process_feature_annotations(
            &gdata.features,
            prior_weight,
            feature_penalty,
            feature_tags,
        )?);
    }
    
    if sample_tags.is_some() {
        gdata.sample_annotations = Some(process_sample_annotations(
            &gdata.samples,
            sample_tags,
        )?);
    }

    Ok(Data { intern: gdata })
}

#[extendr]
#[derive(Clone)]
pub struct Data {
    intern: GData
}

/// @title Data
/// @name Data
/// @description Gpredomics Data object containing feature matrix and labels.
/// 
/// @details
/// The Data object stores the feature matrix (X), labels (y), sample names, feature names,
/// and other metadata needed for machine learning algorithms.
/// 
/// @section Methods:
/// \describe{
///   \item{\code{new()}}{Create a new empty Data object.}
///   \item{\code{get()}}{Returns an R list with all Data fields (X, y, features, samples, etc.).}
///   \item{\code{address()}}{Get memory address of this Data object as a string.}
///   \item{\code{print()}}{Get a formatted string summary of the Data dimensions and content.}
///   \item{\code{train_test_split(test_ratio, stratify_by, seed)}}{Stratified train/test split. Parameters: test_ratio (fraction for test set, 0-1), stratify_by (optional column name in sample annotations for double stratification), seed (optional random seed). Returns a list with train and test Data objects.}
/// }
/// 
/// @export
#[extendr]
impl Data {
    pub fn new() -> Self {
        Self {
            intern: GData::new()
        }
    }

    pub fn get(&self) -> Robj {
        // Convert Data fields to R objects
        let mut pairs = vec![
            ("X", sparse_matrix_to_dataframe(&self.intern.X, &self.intern.samples, &self.intern.features)),
            ("y", Robj::from(self.intern.y.iter().map(|x|{*x as i32}).collect::<Vec<i32>>())),
            ("features", Robj::from(&self.intern.features)),
            ("samples", Robj::from(&self.intern.samples)),
            ("feature_class", sparse_vector_to_vector_usize_u8(&self.intern.feature_class, self.intern.feature_len)),
            // indexes start by 1 in R and need to be adapted
            ("feature_selection", Robj::from(&self.intern.feature_selection
                .iter()
                .map(|x| {x+1})
                .collect::<Vec<usize>>()
            )),
            ("feature_len", Robj::from(&self.intern.feature_len)),
            ("sample_len", Robj::from(&self.intern.sample_len)),
            ("classes", Robj::from(&self.intern.classes))
            
        ];

        // Feature annotations
        if let Some(ref fa) = self.intern.feature_annotations {
            let fa_list = List::from_pairs(vec![
                ("tag_column_names", Robj::from(&fa.tag_column_names)),
                ("feature_tags", feature_tags_to_dataframe_named(&fa.feature_tags, &self.intern.features, &fa.tag_column_names)),
                ("prior_weight", sparse_vector_to_vector_usize_f64(&fa.prior_weight, self.intern.feature_len)),
                ("feature_penalty", sparse_vector_to_vector_usize_f64(&fa.feature_penalty, self.intern.feature_len)),
            ]);
            pairs.push(("feature_annotations", Robj::from(fa_list)));
        } else {
            pairs.push(("feature_annotations", Robj::from(())));
        }
        
        // Sample annotations 
        if let Some(ref sa) = self.intern.sample_annotations {
            let sa_df = if sa.sample_tags.is_empty() {
                Robj::from(())
            } else {
                sample_tags_to_dataframe_named(
                    &sa.sample_tags,
                    &self.intern.samples,
                    &sa.tag_column_names
                )
            };
            
            pairs.push(("sample_annotations", sa_df));
        } else {
            pairs.push(("sample_annotations", Robj::from(())));
        }

        // Return the data as an R list object
        Robj::from(List::from_pairs(pairs))
    }

    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern as *const _)
    }
    
    pub fn print(&self) -> String {
        let addr = self.address();
        
        // === DATA DIMENSIONS ===
        let total_features = self.intern.feature_len;
        let selected_features = self.intern.feature_selection.len();
        let total_samples = self.intern.sample_len;
        
        // === CLASS DISTRIBUTION ===
        let mut class_counts = std::collections::HashMap::new();
        for &class in &self.intern.y {
            *class_counts.entry(class).or_insert(0) += 1;
        }
        
        let class0_count = class_counts.get(&0).copied().unwrap_or(0);
        let class1_count = class_counts.get(&1).copied().unwrap_or(0);
        let other_count = total_samples - class0_count - class1_count;
        
        // === SPARSITY ===
        let total_elements = total_samples * total_features;
        let sparse_elements = self.intern.X.len();
        let sparsity_pct = if total_elements > 0 {
            (1.0 - sparse_elements as f64 / total_elements as f64) * 100.0
        } else { 0.0 };
        
        // === FEATURE SELECTION INFO ===
        let selection_status = if selected_features == total_features {
            "All features".to_string()
        } else if selected_features == 0 {
            "No selection".to_string()
        } else {
            format!("{:.1}% selected", (selected_features as f64 / total_features as f64) * 100.0)
        };

        // === FEATURE ANNOTATIONS ===
        let feature_annot_info = if let Some(ref fa) = self.intern.feature_annotations {
            let mut parts = vec!["Feature annotations present:".to_string()];
            
            if !fa.prior_weight.is_empty() {
                parts.push(format!("  - Prior weights: {} features", fa.prior_weight.len()));
            }
            if !fa.feature_penalty.is_empty() {
                parts.push(format!("  - Feature penalties: {} features", fa.feature_penalty.len()));
            }
            if !fa.feature_tags.is_empty() {
                parts.push(format!("  - Feature tags: {} features, {} columns", fa.feature_tags.len(), fa.tag_column_names.len()));
            }
            
            parts.join("\n")
        } else {
            "Feature annotations: none".to_string()
        };

        // === SAMPLE ANNOTATIONS ===
        let sample_annot_info = if let Some(ref sa) = self.intern.sample_annotations {
            let count = sa.sample_tags.len();
            let cols = sa.tag_column_names.len();
            let example_cols = sa.tag_column_names.iter().take(3).cloned().collect::<Vec<_>>().join(", ");
            
            format!(
                "Sample annotations present:\n  - Samples with tags: {}\n  - Tag columns: {}\n  - Tag columns names: {}",
                count, cols, example_cols
            )
        } else {
            "Sample annotations: none".to_string()
        };
        
        format!(
            "@{}\n\
            Data: {} samples × {} features\n\
            Dimensions\n\
            \tSamples   : {}\n\
            \tFeatures  : {} total\n\
            \tSelected  : {} ({})\n\
            \tSparsity  : {:.1}%\n\
            Class Distribution\n\
            \tClass 0   : {} samples\n\
            \tClass 1   : {} samples\n\
            \tOther     : {} samples\n\
            Memory\n\
            \tSparse Elements: {}\n\
            \n{}\n\n{}",
            addr, total_samples, total_features,
            total_samples,
            total_features,
            selected_features, selection_status,
            sparsity_pct,
            class0_count,
            class1_count,
            other_count,
            sparse_elements,
            feature_annot_info,
            sample_annot_info
        )
    }

    // @title Train/test split on Gpredomics Data
    /// @description Stratified split by class, with optional stratification by sample annotation.
    /// @param test_ratio Fraction of samples in test set (0 < test_ratio < 1).
    /// @param stratify_by Optional column name in sample annotations for double stratification.
    /// @param seed Optional integer seed for reproducibility.
    /// @return A list with elements `train` and `test` (both Data objects).
    pub fn train_test_split(
        &self,
        test_ratio: f64,
        stratify_by: Option<String>,
        seed: Option<u64>,
    ) -> Result<Robj> {
        if test_ratio <= 0.0 || test_ratio >= 1.0 {
            return Err(r_error("test_ratio must be in (0,1)."));
        }

        let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(42));
        let strat = stratify_by.as_ref().map(|s| s.as_str());

        let (train, test) = self.intern.train_test_split(
            test_ratio,
            &mut rng,
            strat,
        );

        let train_data = Data { intern: train };
        let test_data = Data { intern: test };

        Ok(list!(
            train = train_data,
            test  = test_data
        ).into_robj())
    }
}

///////////////////////////////////////////////////////////////
/// Individual object
///////////////////////////////////////////////////////////////

#[derive(Debug, Clone)]
#[extendr]
pub struct Individual {
    intern: GIndividual,
    data: Arc<GData>,
    param: Arc<GParam>
}


impl Individual {
    
    /// Create a new Individual object from a GIndividual.
    /// This internal constructor wraps a Rust GIndividual with shared references\n    
    /// to Data and Param for use in the R interface.
    /// # Arguments
    /// * `ind` - Reference to the internal GIndividual
    /// * `data` - Shared reference to the training data
    /// * `param` - Shared reference to algorithm parameters
    pub fn new(ind: &GIndividual, data: Arc<GData>, param: Arc<GParam>) -> Self {
        Self {
            intern: ind.clone(),
            data,
            param
        }
    }

}

/// @title Individual
/// @name Individual
/// @description A single predictive model from the Gpredomics algorithm
/// @details 
/// Individual represents a single model (individual) in the population.
/// It contains:
/// - Feature set with coefficients (sparse representation)
/// - Performance metrics (AUC, accuracy, sensitivity, specificity, fit value)
/// - Optional threshold confidence interval and rejection rate
/// - Optional beta coefficients for certain model types
/// - Optional additional metrics (MCC, NPV, PPV, F1-score, G-mean)
/// - Genealogical information (parents, generation/epoch)
/// 
/// @section Methods:
/// \describe{
///   \item{\code{get()}}{Get complete individual description including all fields and metrics}
///   \item{\code{get_metrics()}}{Get base metrics (AUC, fit, accuracy, sensitivity, specificity, threshold, rejection_rate)}
///   \item{\code{compute_metrics(data)}}{Compute all metrics including additional ones (MCC, NPV, PPV, F1-score, G-mean) on new data}
///   \item{\code{predict(data)}}{Predict classes and scores for samples in the provided Data object}
///   \item{\code{fit(data, param)}}{Refit individual on new data with complete metric computation including penalties}
///   \item{\code{refit()}}{Refit individual on its associated data using complete Gpredomics logic}
///   \item{\code{evaluate()}}{Compute prediction scores on associated data}
///   \item{\code{to_string()}}{Get debug string representation}
///   \item{\code{address()}}{Get memory address as hexadecimal string}
///   \item{\code{print()}}{Print formatted individual summary to console}
///   \item{\code{set_threshold(threshold)}}{Set prediction threshold for binary classification}
///   \item{\code{compute_importance(data, n_perm, seed, used_only)}}{Compute OOB permutation feature importance}
///   \item{\code{prune_by_threshold(threshold, n_perm, seed, min_k)}}{Remove low-importance features using absolute threshold}
///   \item{\code{prune_by_quantile(quantile, eps, n_perm, seed, min_k)}}{Remove low-importance features using quantile threshold}
///   \item{\code{get_genealogy(experiment, max_depth)}}{Retrieve ancestry tree across generations for visualization}
///   \item{\code{explain_sample(data, sample_index)}}{Explain individual prediction for one sample with feature contributions. Parameters: data (Data object), sample_index (1-based R index). Returns data.frame with Feature, Value, Coefficient, Contribution, CumScore.}
/// }
/// @export
#[extendr]
impl Individual {

    /// Create a new Individual object
    //pub fn new() -> Self {
    //    Self {
    //        intern: GIndividual::new(),
    //        features: Vec::new()
    //    }
    //}

    pub fn get(&self) -> Robj {
        let mut coeff = Vec::new();
        let mut indexes = Vec::new();
        for (index, coefficient) in self.intern.features.iter() {
            coeff.push(*coefficient as i32);     // R integers are 32-bit
            indexes.push(*index as i32 + 1);     // indexes start by 1 in R and need to be adapted
        }
    
        // Build all fields including optionals before constructing the List
        let mut individual_fields = vec![
            ("features", Robj::from(
                self.intern.features.keys()
                    .filter_map(|&idx| self.data.features.get(idx).cloned())
                    .collect::<Vec<String>>()
            )),
            ("auc", Robj::from(self.intern.auc)),
            ("fit", Robj::from(self.intern.fit)),
            ("specificity", Robj::from(self.intern.specificity)),
            ("sensitivity", Robj::from(self.intern.sensitivity)),
            ("accuracy", Robj::from(self.intern.accuracy)),
            ("threshold", Robj::from(self.intern.threshold)),
            ("k", Robj::from(self.intern.k as i32)),
            ("epoch", Robj::from(self.intern.epoch as i32)),
            ("language", Robj::from(self.intern.get_language())),
            ("data_type", Robj::from(self.intern.get_data_type())),
            ("hash", Robj::from(self.intern.hash)),
            ("epsilon", Robj::from(self.intern.epsilon)),
            ("coefficients", Robj::from(coeff)),
            ("rust", self.clone().into_robj())
        ];
    
        // Add "parents"
        let parents_robj = if let Some(parents) = &self.intern.parents {
            Robj::from(parents.clone())
        } else {
            Robj::from(())
        };
        individual_fields.push(("parents", parents_robj));
    
        // Add "betas"
        let betas_robj = if let Some(betas) = &self.intern.betas {
            Robj::from(vec![betas.a, betas.b, betas.c])
        } else {
            Robj::from(())
        };
        
        individual_fields.push(("betas", betas_robj));
    
        let individual_robj = List::from_pairs(individual_fields);
    
        individual_robj.into_robj()

    }

    /// @title Get already computed base metrics
    /// @name Individual$get_metrics
    /// @description Get the already computed base metrics stored in this individual (no computation).
    /// Returns only the core metrics that are always calculated during training: AUC, fit, accuracy, 
    /// sensitivity, specificity, threshold, and rejection_rate (if threshold_ci was computed, otherwise 0.0).
    /// For additional metrics (MCC, NPV, PPV, F1-score, G-mean), use compute_metrics(data).
    /// @return A list with base metrics (auc, fit, accuracy, sensitivity, specificity, threshold, rejection_rate)
    pub fn get_metrics(&self) -> Robj {
        let rejection_rate = self.intern.threshold_ci
            .as_ref()
            .map(|ci| ci.rejection_rate)
            .unwrap_or(0.0);
            
        list!(
            auc = self.intern.auc.into_robj(),
            fit = self.intern.fit.into_robj(),
            accuracy = self.intern.accuracy.into_robj(),
            sensitivity = self.intern.sensitivity.into_robj(),
            specificity = self.intern.specificity.into_robj(),
            threshold = self.intern.threshold.into_robj(),
            rejection_rate = rejection_rate.into_robj(),
        ).into()
    }

    /// @title Compute all metrics on new data
    /// @name Individual$compute_metrics
    /// @description Compute all metrics including base metrics (AUC, threshold, accuracy, sensitivity, specificity)
    /// and additional metrics (MCC, NPV, PPV, F1-score, G-mean, rejection rate) on the provided Data object.
    /// @param data The Data object to compute metrics on
    /// @return A list with all computed metrics (base + additional)
    pub fn compute_metrics(&self, data: &Data) -> Robj {
        let additional_metrics = AdditionalMetrics { 
            mcc: Some(0.0),
            npv: Some(0.0),
            ppv: Some(0.0),
            f1_score: Some(0.0),
            g_mean: Some(0.0),
        };

        let mut individual = self.intern.clone();
        individual.metrics = additional_metrics;
        individual.compute_auc(&data.intern);

        let (accuracy, sensitivity, specificity, rejection_rate, additional) = individual.compute_metrics(&data.intern);
        
        list!(
            auc = individual.auc.into_robj(),
            threshold = individual.threshold.into_robj(),
            accuracy = accuracy.into_robj(),
            sensitivity = sensitivity.into_robj(),
            specificity = specificity.into_robj(),
            rejection_rate = rejection_rate.into_robj(),
            mcc = additional.mcc.into_robj(),
            npv = additional.npv.into_robj(),
            ppv = additional.ppv.into_robj(),
            f1_score = additional.f1_score.into_robj(),
            g_mean = additional.g_mean.into_robj(),
        ).into()
    }

    /// @title Fit individual on new data
    /// @name Individual$fit
    /// @description Fit the individual on new data by recomputing all metrics using the complete Gpredomics fitting logic.
    /// This includes computing AUC, threshold (with optional confidence interval and rejection rate), accuracy, sensitivity, 
    /// specificity, fit value (with penalties), and additional metrics based on the param configuration.
    /// @param data The Data object to fit on
    /// @param param The Param object containing fitting configuration (fit function, penalties, bootstrap settings, etc.)
    /// @return A new Individual with all metrics fully recomputed
    /// @export
    pub fn fit(&self, data: &Data, param: &Param) -> Individual {
        let mut gi = self.clone();
        
        // Create a temporary population with just this individual to use Population::fit
        let mut temp_pop = GPopulation::new();
        temp_pop.individuals.push(gi.intern.clone());
        
        // Use the full Population::fit logic (without GPU)
        temp_pop.fit(&data.intern, &mut None, &None, &None, &param.intern);
        
        // Extract the fitted individual
        gi.intern = temp_pop.individuals[0].clone();
        
        Individual {
            intern: gi.intern,
            data: Arc::new(data.intern.clone()),
            param: Arc::new(param.intern.clone())
        }
    }

    /// @title Refit individual on associated data
    /// @name Individual$refit
    /// @description Refit the individual on its associated data by recomputing all metrics using the complete Gpredomics fitting logic.
    /// This includes computing AUC, threshold (with optional confidence interval and rejection rate), accuracy, sensitivity, 
    /// specificity, fit value (with penalties), and additional metrics based on the param configuration.
    /// @return A new Individual with all metrics fully recomputed
    /// @export
    pub fn refit(&self) -> Individual {
        let mut gi = self.clone();
        
        // Create a temporary population with just this individual to use Population::fit
        let mut temp_pop = GPopulation::new();
        temp_pop.individuals.push(gi.intern.clone());
        
        // Use the full Population::fit logic (without GPU)
        temp_pop.fit(&self.data, &mut None, &None, &None, &self.param);
        
        // Extract the fitted individual
        gi.intern = temp_pop.individuals[0].clone();
        gi
    }

    /// @title Evaluate individual
    /// @name Individual$evaluate
    /// @description Compute individual score (evaluate on its associated data).
    /// @return R object containing the individual score
    pub fn evaluate(&self) -> Robj {
        self.intern.evaluate(&self.data).into_robj()
    }

    /// @title Predict classes and scores
    /// @name Individual$predict
    /// @description Predict classes and scores on the provided Data object.
    /// @param data The Data object to predict on
    /// @return A list with two elements: class (predicted classes) and score (predicted scores)
    pub fn predict(&self, data: &Data) -> Robj {
        let (classes, scores) = self.intern.evaluate_class_and_score(&data.intern);
        list!(class=(classes.iter().map(|x| {*x as i32})).collect::<Vec<i32>>().into_robj(), score=scores.into_robj()).into()
    }

    /// @title Individual to string
    /// @name Individual$to_string
    /// @description Return a string representation of the Individual (debug representation).
    /// @return String representation of the Individual
    pub fn to_string(&self) -> String {
        format!("{:?}",&self.intern)
    }

    /// @title Individual memory address
    /// @name Individual$address
    /// @description Get memory address of this Individual object.
    /// @return String representing the memory address
    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern as *const _)
    }

    // @title Print as Gpredomics style
    // @name Individual$print
    /// @title Print individual (formatted)
    /// @name Individual$print
    /// @description Print the individual in gpredomics style to R console.
    pub fn print(&self) {
        let addr = self.address();
        rprintln!("@{}", addr);
        let s = &self.intern.display(&self.data, None, &self.param.general.algo, self.param.general.threshold_ci_alpha);
        for line in s.lines() {
            rprintln!("{}", line);
        }
    }

    /// @title Set threshold
    /// @name Individual$set_threshold
    /// @description Set the threshold of the individual used for binary predictions.
    /// @param threshold The new threshold value
    pub fn set_threshold(&self, threshold: f64) -> Individual {
        let mut new_intern = self.intern.clone();

        if let Some(ci) = &self.intern.threshold_ci {
            r_warning(&format!(
                "Overriding existing threshold CI [{:.3}, {:.3}] with new threshold {:.3}.",
                ci.lower, ci.upper, threshold
            ));
            new_intern.threshold_ci = None;
        }
        
        new_intern.threshold = threshold;
        (new_intern.accuracy, new_intern.sensitivity, new_intern.specificity, _, new_intern.metrics) = new_intern.compute_metrics(&self.data);

        Individual {
            intern: new_intern,
            data: Arc::clone(&self.data),
            param: Arc::clone(&self.param),
        }
    }

    /// @title Compute feature importance
    /// @name Individual$compute_importance
    /// @description Compute feature importance for this individual on the provided Data using permutation (MDA-like).
    /// @param data The Data object to compute importance on
    /// @param n_perm Number of permutations to perform (default 1000)
    /// @param seed Optional seed for random number generation
    /// @param used_only Whether to compute importance only for features used in the individual (default true)
    /// @return A named numeric vector of feature importances
    pub fn compute_importance(&self, data: &Data, n_perm: Option<i32>, seed: Option<u64>, used_only: Option<bool>) -> Result<Robj> {
        let permutations = n_perm.unwrap_or(1000).max(1) as usize;
        let used_only = used_only.unwrap_or(true);

        let features_to_process: Vec<usize> = if used_only {
            let mut v: Vec<usize> = self.intern.features.keys().copied().collect();
            v.sort_unstable();
            v
        } else {
            (0..self.data.feature_len).collect()
        };

        let mut feature_seeds: HashMap<usize, Vec<u64>> = HashMap::with_capacity(features_to_process.len());
        let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(4815162342));
        for &f in &features_to_process {
            let mut seeds = Vec::with_capacity(permutations);
            for _ in 0..permutations {
                seeds.push(rng.next_u64());
            }
            feature_seeds.insert(f, seeds);
        }

        let ic: ImportanceCollection = self.intern.compute_oob_feature_importance(&data.intern, permutations, &features_to_process, &feature_seeds); 

        let mut feature = Vec::with_capacity(ic.importances.len());
        let mut importance = Vec::with_capacity(ic.importances.len());
        for imp in &ic.importances {
            let fname = data.intern.features
                .get(imp.feature_idx)
                .cloned()
                .unwrap_or_else(|| format!("f{}", imp.feature_idx));
            feature.push(fname);
            importance.push(imp.importance);
        }
        let mut v = Doubles::from_values(importance);
        let _ = v.set_names(feature);
        Ok(v.into_robj())
    }

    /// @title Prune individual by importance threshold
    /// @name Individual$prune_by_threshold
    /// @description Prune (remove) features from this individual based on OOB permutation importance.
    /// Features with importance below the specified threshold are removed, while ensuring at least `min_k` features remain.
    /// Metrics are recomputed after pruning using the refit() method.
    /// Returns a new pruned Individual with updated metrics.
    /// @param threshold Importance threshold - drops features with importance < threshold
    /// @param n_perm Number of permutations to perform (default 100)
    /// @param seed Seed for random number generation (default 4815162342)
    /// @param min_k Minimum number of features to keep (default 1)
    /// @return A new pruned Individual with recomputed metrics
    /// @section Methods:
    /// This method creates a new Individual with low-importance features removed.
    /// It uses OOB permutation importance (MDA) to rank features, then recalculates all metrics.
    /// @export
    pub fn prune_by_threshold(
        &self,
        threshold: f64,
        n_perm: Option<i32>,
        seed: Option<u64>,
        min_k: Option<i32>,
    ) -> Individual {
        let permutations = n_perm.unwrap_or(100).max(1) as usize;
        let rng_seed = seed.unwrap_or(4815162342);
        let min_features = min_k.unwrap_or(1).max(1) as usize;
        
        // Clone the individual
        let mut pruned = self.clone();
        
        pruned.intern.prune_by_importance(
            &self.data,
            permutations,
            rng_seed,
            Some(threshold),
            None,
            min_features,
        );
        
        pruned.refit()
    }

    /// @title Prune individual by importance quantile
    /// @name Individual$prune_by_quantile
    /// @description Prune (remove) features from this individual based on OOB permutation importance.
    /// Features with importance below quantile(q) - eps are removed, while ensuring at least `min_k` features remain.
    /// Metrics are recomputed after pruning using the refit() method.
    /// Returns a new pruned Individual with updated metrics.
    /// @param exp Experiment object (currently not used but kept for compatibility)
    /// @param quantile Quantile value between 0 and 1 (e.g., 0.25 for 25th percentile)
    /// @param eps Optional epsilon value to subtract from quantile (default 0.0)
    /// @param n_perm Number of permutations to perform (default 100)
    /// @param seed Seed for random number generation (default 4815162342)
    /// @param min_k Minimum number of features to keep (default 1)
    /// @return A new pruned Individual with recomputed metrics
    /// @section Methods:
    /// This method creates a new Individual with low-importance features removed.
    /// It uses OOB permutation importance (MDA) to rank features, then recalculates all metrics.
    /// The cutoff is computed as: quantile(importance, q) - eps
    /// @export
    pub fn prune_by_quantile(
        &self,
        quantile: f64,
        eps: Option<f64>,
        n_perm: Option<i32>,
        seed: Option<u64>,
        min_k: Option<i32>,
    ) -> Individual {
        let permutations = n_perm.unwrap_or(100).max(1) as usize;
        let rng_seed = seed.unwrap_or(4815162342);
        let min_features = min_k.unwrap_or(1).max(1) as usize;
        let epsilon = eps.unwrap_or(0.0);
        
        // Clone the individual
        let mut pruned = self.clone();
        
        pruned.intern.prune_by_importance(
            &self.data,
            permutations,
            rng_seed,
            None,
            Some((quantile, epsilon)),
            min_features,
        );
        
        pruned.refit()
    }

    /// @title Get individual genealogy
    /// @name Individual$get_genealogy
    /// @description Retrieve the genealogy (ancestry tree) of this individual across generations.
    /// Returns a list with two data.frames: `nodes` (containing complete individual information with all fields)
    /// and `edges` (parent-child relationships), plus metadata. Includes optional fields (betas, threshold_ci)
    /// when present. Designed for igraph/ggraph visualization and comprehensive genealogical analysis.
    /// @param experiment The Experiment object containing all generations
    /// @param max_depth Maximum depth to traverse in the genealogy tree (default: 10)
    /// @return A list with components:
    /// - `nodes`: data.frame with columns: id, label, hash, auc, fit, sensitivity, specificity, accuracy, 
    ///   threshold, k, generation, language, data_type, epsilon, plus optional columns: beta_a, beta_b, beta_c,
    ///   ci_lower, ci_upper, ci_rejection_rate, mcc, f1_score, npv, ppv, g_mean (when available)
    /// - `edges`: data.frame with columns: from, to, depth
    /// - `metadata`: list with statistics
    /// @export
    pub fn get_genealogy(&self, experiment: &Experiment, max_depth: Option<i32>) -> Result<Robj> {
        let depth = max_depth.unwrap_or(10).max(1) as usize;

        // Collect all populations and individuals by hash (across all generations)
        let mut pops_vec: Vec<GPopulation> = Vec::new();
        let mut ind_by_hash: HashMap<u64, &GIndividual> = HashMap::new();

        for pop_collection in &experiment.intern.collections {
            for pop in pop_collection.iter() {
                for ind in &pop.individuals {
                    ind_by_hash.entry(ind.hash).or_insert(ind);
                }
                pops_vec.push(pop.clone());
            }
        }

        let genealogy_map = self.intern.get_genealogy(&pops_vec, depth);

        // Vectors for required fields
        let mut node_id: Vec<String> = Vec::new();
        let mut node_label: Vec<String> = Vec::new();
        let mut node_hash: Vec<String> = Vec::new();
        let mut node_auc: Vec<f64> = Vec::new();
        let mut node_fit: Vec<f64> = Vec::new();
        let mut node_sensitivity: Vec<f64> = Vec::new();
        let mut node_specificity: Vec<f64> = Vec::new();
        let mut node_accuracy: Vec<f64> = Vec::new();
        let mut node_threshold: Vec<f64> = Vec::new();
        let mut node_k: Vec<i32> = Vec::new();
        let mut node_gen: Vec<i32> = Vec::new();
        let mut node_language: Vec<String> = Vec::new();
        let mut node_data_type: Vec<String> = Vec::new();
        let mut node_epsilon: Vec<f64> = Vec::new();

        // Vectors for optional fields (check existence in first individual)

        // Correction : on vérifie la présence de threshold_ci sur tous les individus
        let first_has_betas = ind_by_hash.values().any(|ind| ind.betas.is_some());
        let first_has_ci = ind_by_hash.values().any(|ind| ind.threshold_ci.is_some());
        let first_has_metrics = ind_by_hash.values().any(|ind| {
            ind.metrics.mcc.is_some() || ind.metrics.f1_score.is_some() || 
            ind.metrics.npv.is_some() || ind.metrics.ppv.is_some() || ind.metrics.g_mean.is_some()
        });

    let mut node_beta_a: Vec<Robj> = Vec::new();
    let mut node_beta_b: Vec<Robj> = Vec::new();
    let mut node_beta_c: Vec<Robj> = Vec::new();
    // Correction : colonnes optionnelles plates (f64, NA = f64::NAN)
    let mut node_ci_lower: Vec<f64> = Vec::new();
    let mut node_ci_upper: Vec<f64> = Vec::new();
    let mut node_ci_rejection: Vec<f64> = Vec::new();
        let mut node_mcc: Vec<Robj> = Vec::new();
        let mut node_f1_score: Vec<Robj> = Vec::new();
        let mut node_npv: Vec<Robj> = Vec::new();
        let mut node_ppv: Vec<Robj> = Vec::new();
        let mut node_g_mean: Vec<Robj> = Vec::new();

        // Vectors for edges
        let mut edge_from: Vec<String> = Vec::new();
        let mut edge_to: Vec<String> = Vec::new();
        let mut edge_depth: Vec<i32> = Vec::new();

    for ((hash, parents_opt), depth_set) in genealogy_map {
            let id = format!("H{:016x}", hash);
            node_id.push(id.clone());
            node_label.push(format!("H#{}", &id[1..9]));
            node_hash.push(format!("{:016x}", hash));

            // Extract fields from the individual
            if let Some(ind) = ind_by_hash.get(&hash) {
                node_auc.push(ind.auc);
                node_fit.push(ind.fit);
                node_sensitivity.push(ind.sensitivity);
                node_specificity.push(ind.specificity);
                node_accuracy.push(ind.accuracy);
                node_threshold.push(ind.threshold);
                node_k.push(ind.k as i32);
                node_gen.push(ind.epoch as i32);
                node_language.push(ind.get_language().to_string());
                node_data_type.push(ind.get_data_type().to_string());
                node_epsilon.push(ind.epsilon);

                // Optional fields - betas
                if first_has_betas {
                    if let Some(ref betas) = ind.betas {
                        node_beta_a.push(Robj::from(betas.a));
                        node_beta_b.push(Robj::from(betas.b));
                        node_beta_c.push(Robj::from(betas.c));
                    } else {
                        node_beta_a.push(Robj::from(()));
                        node_beta_b.push(Robj::from(()));
                        node_beta_c.push(Robj::from(()));
                    }
                }

                // Correction : toujours une seule colonne ci_rejection_rate, NA si absent
                if first_has_ci {
                    if let Some(ref ci) = ind.threshold_ci {
                        node_ci_lower.push(ci.lower);
                        node_ci_upper.push(ci.upper);
                        node_ci_rejection.push(ci.rejection_rate);
                    } else {
                        node_ci_lower.push(f64::NAN);
                        node_ci_upper.push(f64::NAN);
                        node_ci_rejection.push(f64::NAN);
                    }
                }

                // Additional metrics
                if first_has_metrics {
                    node_mcc.push(match ind.metrics.mcc {
                        Some(v) => Robj::from(v),
                        None => Robj::from(())
                    });
                    node_f1_score.push(match ind.metrics.f1_score {
                        Some(v) => Robj::from(v),
                        None => Robj::from(())
                    });
                    node_npv.push(match ind.metrics.npv {
                        Some(v) => Robj::from(v),
                        None => Robj::from(())
                    });
                    node_ppv.push(match ind.metrics.ppv {
                        Some(v) => Robj::from(v),
                        None => Robj::from(())
                    });
                    node_g_mean.push(match ind.metrics.g_mean {
                        Some(v) => Robj::from(v),
                        None => Robj::from(())
                    });
                }
            } else {
                // Default values if individual not found
                node_auc.push(0.0);
                node_fit.push(0.0);
                node_sensitivity.push(0.0);
                node_specificity.push(0.0);
                node_accuracy.push(0.0);
                node_threshold.push(0.0);
                node_k.push(0);
                node_gen.push(-1);
                node_language.push("Unknown".to_string());
                node_data_type.push("Unknown".to_string());
                node_epsilon.push(0.0);

                if first_has_betas {
                    node_beta_a.push(Robj::from(()));
                    node_beta_b.push(Robj::from(()));
                    node_beta_c.push(Robj::from(()));
                }
                if first_has_ci {
                    node_ci_lower.push(f64::NAN);
                    node_ci_upper.push(f64::NAN);
                    node_ci_rejection.push(f64::NAN);
                }
                if first_has_metrics {
                    node_mcc.push(Robj::from(()));
                    node_f1_score.push(Robj::from(()));
                    node_npv.push(Robj::from(()));
                    node_ppv.push(Robj::from(()));
                    node_g_mean.push(Robj::from(()));
                }
            }

            // Build edges
            if let Some(parents) = parents_opt {
                for &p in &parents {
                    let pid = format!("H{:016x}", p);
                    for &d in &depth_set {
                        edge_from.push(pid.clone());
                        edge_to.push(id.clone());
                        edge_depth.push(d as i32);
                    }
                }
            }
        }


        // Construction explicite de la liste pour éviter les colonnes dynamiques
        let mut nodes_pairs = vec![
            ("id", Robj::from(node_id.clone())),
            ("label", Robj::from(node_label)),
            ("hash", Robj::from(node_hash)),
            ("auc", Robj::from(node_auc.clone())),
            ("fit", Robj::from(node_fit)),
            ("sensitivity", Robj::from(node_sensitivity)),
            ("specificity", Robj::from(node_specificity)),
            ("accuracy", Robj::from(node_accuracy)),
            ("threshold", Robj::from(node_threshold)),
            ("k", Robj::from(node_k.clone())),
            ("generation", Robj::from(node_gen)),
            ("language", Robj::from(node_language)),
            ("data_type", Robj::from(node_data_type)),
            ("epsilon", Robj::from(node_epsilon)),
        ];

        if first_has_betas {
            nodes_pairs.push(("beta_a", Robj::from(node_beta_a)));
            nodes_pairs.push(("beta_b", Robj::from(node_beta_b)));
            nodes_pairs.push(("beta_c", Robj::from(node_beta_c)));
        }
        if first_has_ci {
            nodes_pairs.push(("ci_lower", Robj::from(node_ci_lower.clone())));
            nodes_pairs.push(("ci_upper", Robj::from(node_ci_upper.clone())));
            nodes_pairs.push(("ci_rejection_rate", Robj::from(node_ci_rejection.clone())));
        }
        if first_has_metrics {
            nodes_pairs.push(("mcc", Robj::from(node_mcc)));
            nodes_pairs.push(("f1_score", Robj::from(node_f1_score)));
            nodes_pairs.push(("npv", Robj::from(node_npv)));
            nodes_pairs.push(("ppv", Robj::from(node_ppv)));
            nodes_pairs.push(("g_mean", Robj::from(node_g_mean)));
        }

        let nodes_list = List::from_pairs(nodes_pairs);
        let nodes_df = call!("as.data.frame", nodes_list)?;

        let edges_df = data_frame!(
            from = edge_from,
            to = edge_to,
            depth = edge_depth
        );

        // Enriched metadata
        let metadata = list!(
            n_nodes = nodes_df.nrows().try_into().unwrap_or(0),
            n_edges = edges_df.nrows().try_into().unwrap_or(0),
            max_depth = depth as i32,
            root_id = format!("H{:016x}", self.intern.hash),
            mean_auc = node_auc.iter().sum::<f64>() / node_auc.len().max(1) as f64,
            mean_k = node_k.iter().sum::<i32>() as f64 / node_k.len().max(1) as f64,
            has_betas = first_has_betas,
            has_threshold_ci = first_has_ci,
            has_additional_metrics = first_has_metrics
        );

        Ok(list!(nodes = nodes_df, edges = edges_df, metadata = metadata).into())
    }

    /// @title Explain Individual Prediction for One Sample
    /// @name Individual$explain_sample
    /// @param data Data object used for evaluation.
    /// @param sample_index **1-based** index of the sample (R convention: 1 to n_samples).
    /// @return An R data.frame with columns Feature, Feature_Index, Value, 
    ///         Coefficient, Contribution, CumScore.
    /// @export
    pub fn explain_sample(&self, data: &Data, sample_index: i32) -> Result<Robj> {
        
        if sample_index <= 0 || sample_index > data.intern.sample_len as i32 {
            return Err(r_error(&format!(
                "sample_index {} out of range (1..{}).",
                sample_index,
                data.intern.sample_len
            )));
        }
        let sample_idx_0based = (sample_index - 1) as usize;

        // Sort features by index for stable display
        let mut feats: Vec<(&usize, &i8)> = self.intern.features.iter().collect();
        feats.sort_by_key(|(idx, _)| *idx);

        let mut names: Vec<String> = Vec::with_capacity(feats.len());
        let mut feat_idx: Vec<i32> = Vec::with_capacity(feats.len());
        let mut values: Vec<f64> = Vec::with_capacity(feats.len());
        let mut coefs: Vec<f64> = Vec::with_capacity(feats.len());
        let mut contribs: Vec<f64> = Vec::with_capacity(feats.len());
        let mut cums: Vec<f64> = Vec::with_capacity(feats.len());

        let mut cum_score = 0.0_f64;

        for (idx, coef) in feats {
            let j = *idx;
            let coef_f = *coef as f64;

            let x_ij = data
                .intern
                .X
                .get(&(sample_idx_0based, j))
                .copied()
                .unwrap_or(0.0);

            // Apply transformation logic
            let transformed = match self.intern.data_type {
                RAW_TYPE => x_ij,
                PREVALENCE_TYPE => {
                    if x_ij >= self.intern.epsilon { 1.0 } else { 0.0 }
                }
                LOG_TYPE => {
                    if x_ij > 0.0 {
                        (x_ij + self.intern.epsilon).ln()
                    } else {
                        0.0
                    }
                }
                _ => x_ij,
            };

            let contrib = match self.intern.language {
                RATIO_LANG => coef_f * transformed,
                _ => coef_f * transformed,
            };

            cum_score += contrib;

            let fname = if j < data.intern.features.len() {
                data.intern.features[j].clone()
            } else {
                format!("Feature_{}", j)
            };

            names.push(fname);
            feat_idx.push((j + 1) as i32);  
            values.push(x_ij);
            coefs.push(coef_f);
            contribs.push(contrib);
            cums.push(cum_score);
        }

        let df = data_frame!(
            Feature = names,
            Feature_Index = feat_idx,
            Value = values,
            Coefficient = coefs,
            Contribution = contribs,
            CumScore = cums
        );

        Ok(df.into())
    }

}

///////////////////////////////////////////////////////////////
/// Experiment object
///////////////////////////////////////////////////////////////

#[extendr]
pub struct Experiment {
    intern: GExperiment,
    train_data_arc: Arc<GData>,
    param_arc: Arc<GParam>,
}

/// @title Experiment
/// @name Experiment
/// @description Complete experiment container for Gpredomics runs
/// @details 
/// Experiment encapsulates an entire Gpredomics run including:
/// - Training and optional test data
/// - Algorithm parameters
/// - Evolution history (populations across generations)
/// - Cross-validation fold information (if applicable)
/// - Final best population and optional Jury
/// - Execution metadata (timestamp, version, execution time)
/// 
/// @section Methods:
/// \describe{
///   \item{\code{individual(generation, order)}}{Retrieves a full description of an individual from a specified generation and order. Parameters: generation (i32 generation index), order (i32 order within generation). Returns Individual object.}
///   \item{\code{test_data()}}{Retrieves the test data associated with the experiment. Returns Data object representing the test data.}
///   \item{\code{train_data()}}{Retrieves the training data associated with the experiment. Returns Data object representing the training data.}
///   \item{\code{get_data_robj(train)}}{Retrieves the data associated with the experiment as an R object. Parameter: train (logical; if TRUE, returns training data, otherwise test data). Returns R object representing the data.}
///   \item{\code{get_data(train)}}{Retrieves the data associated with the experiment as a Data object. Parameter: train (logical; if TRUE, returns training data, otherwise test data). Returns Data object.}
///   \item{\code{get_generation(generation)}}{Retrieves descriptions of all individuals from a specified generation. Parameter: generation (i32 generation index). Returns R list object encapsulating features and metrics of all individuals in the generation.}
///   \item{\code{get_fold_data(fold, train)}}{Get fold data for cross-validation. Parameters: fold (fold index 0-based), train (logical: if TRUE returns training data, if FALSE returns validation data). Returns Data object for the specified fold.}
///   \item{\code{get_fold_generation(fold, generation, train)}}{Get population from a specific CV fold and generation. Parameters: fold (fold index 0-based), generation (generation index 0-based), train (logical: if TRUE uses training data, if FALSE uses validation data, default TRUE). Returns Population object.}
///   \item{\code{get_n_folds()}}{Get the number of cross-validation folds in the experiment. Returns number of CV folds as integer.}
///   \item{\code{get_best_population()}}{Get the final/best population from the experiment. Returns Population object containing the best individuals.}
///   \item{\code{compute_cv_importance(n_perm, aggregation, scaled, seed, compact)}}{Compute cross-validated feature importance (MDA-like) aggregated across CV folds. Parameters: n_perm (number of permutations, default 1000), aggregation ("mean" or "median", default "mean"), scaled (whether to scale, default TRUE), seed (random seed, default 4815162342), compact (return vector if TRUE or data.frame if FALSE, default FALSE). Returns either data.frame or named vector.}
///   \item{\code{compute_cv_importance_matrix(n_perm, used_only, seed, aggregation, scaled)}}{Compute per-fold CV population-level importance matrix. Parameters: n_perm (number of permutations, default 1000), used_only (restrict to FBM features, default TRUE), seed (base seed for RNG, default 4815162342), aggregation (population aggregation method: "mean" or "median", default "mean"), scaled (scale importances, default TRUE). Returns numeric matrix: n_features x n_folds.}
///   \item{\code{generation_number()}}{Get the number of generations in the population. Returns number of generations.}
///   \item{\code{population_size(generation)}}{Get the size (number of individuals) of a certain generation in a Population. Parameter: generation (i32 generation index). Returns integer size.}
///   \item{\code{load_data(x_path, y_path)}}{Load an external dataset to evaluate the model. Parameters: x_path (path to X data file), y_path (path to y data file). Returns Data object containing the loaded data.}
///   \item{\code{get_param()}}{Get the param object associated with the experiment. Returns Param object containing the experiment parameters.}
///   \item{\code{get_jury()}}{Get the jury object associated with the experiment. Returns Jury object containing the jury details.}
///   \item{\code{load(path)}}{Load a serialized experiment. Parameter: path (path to the experiment file). Returns loaded Experiment object.}
///   \item{\code{save(path)}}{Save an experiment. Parameter: path (path to save the experiment).}
///   \item{\code{get_population(generation)}}{Extract population from experiment, optionally specifying generation number. Parameter: generation (optional generation number, 0-based; if None, returns final population). Returns Population object for the specified generation or final population.}
///   \item{\code{address()}}{Get memory address of this Experiment object. Returns string representing the memory address.}
///   \item{\code{print()}}{Get print of this Experiment. Returns string representing the Experiment summary.}
///   \item{\code{get_history(data, scope)}}{Get training history with specified scope. Parameters: data (Data object to evaluate on), scope ("best", "fbm", "top5", or "all", default "best"). Returns data.frame with metrics per generation (and per fold if CV).}
/// }
/// @export
#[extendr]
impl Experiment {
    /// @title Get individual from generation
    /// @name Experiment$individual
    /// @description Retrieves a full description of an individual from a specified generation and order.
    /// @param generation An i32 specifying the generation index.
    /// @param order An i32 specifying the order of the individual within the generation.
    /// @return Individual object containing individual details
    pub fn individual(&self, generation: i32, order: i32) -> Result<Individual> {
        let gen_idx = (generation - 1) as usize;
        let ord_idx = (order - 1) as usize;
        
        if self.intern.collections.is_empty() || gen_idx >= self.intern.collections[0].len() {
            return Err(r_error(&format!("Generation {} out of bounds", generation)));
        }
        
        if ord_idx >= self.intern.collections[0][gen_idx].individuals.len() {
            return Err(r_error(&format!("Order {} out of bounds for generation {}", order, generation)));
        }
        
        Ok(Individual::new(
            &self.intern.collections[0][gen_idx].individuals[ord_idx],
            Arc::clone(&self.train_data_arc),
            Arc::clone(&self.param_arc)
        ))
    }

    /// @title Get test data
    /// @name Experiment$test_data
    /// @description Retrieves the test data associated with the experiment.
    /// @return Data object representing the test data.
    pub fn test_data(&self) -> Result<Data> {
        if let Some(test_data) = &self.intern.test_data {
            Ok(Data {
                intern: test_data.clone()
            })
        } else {
            Err(r_error("No test data attached to this experiment"))
        }
    }

    /// @title Get training data
    /// @name Experiment$train_data
    /// @description Retrieves the training data associated with the experiment.
    /// @return Data object representing the training data.
    pub fn train_data(&self) -> Data {
        Data {
            intern: self.intern.train_data.clone()
        }
    }

    /// @title Get data as R object
    /// @name Experiment$get_data_robj
    /// @description Retrieves the data associated with the experiment as an R object.
    /// @param train Logical; if TRUE, returns training data, otherwise test data
    /// @return R object representing the data
    pub fn get_data_robj(&self, train:bool) -> Robj {
        if train {self.train_data().get()}
        else {self.test_data().expect("No test data attached to this experiment").get()}
    }

    /// @title Get data as Data object
    /// @name Experiment$get_data
    /// @description Retrieves the data associated with the experiment as a Data object.
    /// @param train Logical; if TRUE, returns training data, otherwise test data
    /// @return Data object
    pub fn get_data(&self, train:bool) -> Result<Data> {
        if train {Ok(self.train_data())}
        else {self.test_data()}
    }
    
    /// @title Get generation
    /// @name Experiment$get_generation
    /// @description Retrieves descriptions of all individuals from a specified generation.
    /// @param generation An i32 specifying the generation index.
    /// @return An R list object (Robj) encapsulating features and metrics of all individuals in the generation.
    pub fn get_generation(&self, generation: i32) -> Result<Population> {
        if !self.intern.cv_folds_ids.is_none() {
            return Err(r_error("CV folds available. To extract a generation, use get_fold_generation(fold, generation, on_train) instead."));
        }

        if generation <=0 || generation as usize > self.intern.collections[0].len() {
            return Err(r_error(format!("Generation index out of bounds for this experiment ({} available generations).", self.intern.collections[0].len())));
        }

        let pop = if self.intern.collections.len() > 0 {
            self.intern.collections[0][generation as usize -1].clone()
        } else {
            return Err(r_error("No collection found in this experiment"));
        };

        Ok(Population {
            intern: pop, 
            data: Arc::clone(&self.train_data_arc),
            param: Arc::clone(&self.param_arc)
        })
    } 
    
    /// Get fold data for cross-validation
    /// @title Get CV fold data
    /// @name Experiment$get_fold_data
    /// @param fold Fold index (0-based)
    /// @param train Logical: if TRUE, returns training data; if FALSE, returns validation data
    /// @return Data object for the specified fold
    /// @export
    pub fn get_fold_data(&self, fold: i32, train: bool) -> Result<Data> {
        if self.intern.cv_folds_ids.is_none() {
            return Err(r_error("No CV folds available. Run with param$cv = TRUE"));
        }
        
        let cv_folds_ids = self.intern.cv_folds_ids.as_ref().unwrap();
        
        if fold <= 0 || fold as usize > cv_folds_ids.len() {
            return Err(r_error(format!("Fold index out of bounds for this experiment ({} available folds).", cv_folds_ids.len())));
        }
        
        let fold_idx = (fold - 1) as usize;
        
        let (train_ids, valid_ids) = &cv_folds_ids[fold_idx];
        let sample_ids = if train { train_ids } else { valid_ids };
        
        let fold_data = reconstruct_fold_data(&self.intern.train_data, sample_ids);
        
        Ok(Data { intern: fold_data })
    }

    /// Get population from a specific CV fold and generation
    /// @title Get CV fold generation population
    /// @name Experiment$get_fold_generation
    /// @param fold Fold index (0-based)
    /// @param generation Generation index (0-based)
    /// @param train Logical: if TRUE, uses training data; if FALSE, uses validation data (default: TRUE)
    /// @return Population object for the specified fold and generation
    /// @export
     pub fn get_fold_generation(&self, fold: i32, generation: i32, train: Option<bool>) -> Result<Population> {
        if self.intern.collections.is_empty() {
            return Err(r_error("No CV collections available. Run with param$cv = TRUE"));
        }

        if fold <= 0 || fold as usize > self.intern.collections.len() {
            return Err(r_error(format!("Fold index out of bounds for this experiment ({} available folds).", self.intern.collections.len())));
        }
        
        let fold_idx = (fold - 1) as usize;
        
        let use_train = train.unwrap_or(true);
        
        let fold_data = self.get_fold_data(fold, !use_train)?.intern;
        let fold_data_arc = Arc::new(fold_data);
        
        let fold_collections = &self.intern.collections[fold_idx];

        if generation <= 0 || generation as usize > fold_collections.len() {
            return Err(r_error(format!("Generation index out of bounds for this fold ({} available generations).", fold_collections.len())));
        }
        
        let gen_idx = (generation - 1) as usize;
        
        Ok(Population::new(
            fold_collections[gen_idx].clone(),
            fold_data_arc,  
            Arc::clone(&self.param_arc)
        ))
    }
    
    /// Get number of CV folds
    /// @title Get number of CV folds
    /// @name Experiment$get_n_folds
    /// @description Get the number of cross-validation folds in the experiment
    /// @return Number of CV folds as integer
    /// @export
    pub fn get_n_folds(&self) -> i32 {
        self.intern.collections.len() as i32
    }

    /// Get best population
    /// @title Get best population
    /// @name Experiment$get_best_population
    /// @description Get the final/best population from the experiment
    /// @return Population object containing the best individuals
    /// @export
    pub fn get_best_population(&self) -> Result<Population> {
        if let Some(final_pop) = &self.intern.final_population {
            Ok(Population::new(
                final_pop.clone(),
                Arc::clone(&self.train_data_arc),  
                Arc::clone(&self.param_arc)
            ))
        } else {
            return Err(r_error("No final population found in this experiment"));
        }
    }

    /// @title Compute cross-validated feature importance
    /// @name Experiment$compute_cv_importance
    /// @description
    /// Compute cross-validated feature importance (MDA-like) aggregated across CV folds.  
    /// This uses the native CV::compute_cv_oob_feature_importance logic on the original
    /// training data and CV folds stored in the Experiment.
    ///
    /// @param n_perm Number of permutations (default: 1000)
    /// @param aggregation Aggregation method across folds: "mean" (default) or "median"
    /// @param scaled Whether to scale importances (default: TRUE)
    /// @param seed Random seed for reproducibility (default: 4815162342)
    /// @param compact Whether to return a compact named vector (TRUE) or full data.frame (FALSE, default)
    /// @return Either:
    ///   - data.frame with columns: feature, importance, dispersion, prevalence
    ///   - or named numeric vector if compact = TRUE
    /// @export
    pub fn compute_cv_importance(
        &self,
        n_perm: Option<i32>,
        aggregation: Option<String>,
        scaled: Option<bool>,
        seed: Option<u64>,
        compact: Option<bool>,
    ) -> extendr_api::Result<Robj> {
        // 1) Sanity checks: CV must be available
        if self.intern.cv_folds_ids.is_none() || self.intern.collections.is_empty() {
            return Err(r_error("No CV folds available. Run with param$cv = TRUE"));
        }

        let permutations = n_perm.unwrap_or(1000).max(1) as usize;
        let use_median = matches!(aggregation.as_deref(), Some("median"));
        let agg = match aggregation.as_deref() {
            Some("median") => ImportanceAggregation::median,
            Some("mean") | None => ImportanceAggregation::mean,
            Some(x) => return Err(r_error(format!("Unknown aggregation '{}'", x))),
        };
        let scaled = scaled.unwrap_or(true);
        let compact = compact.unwrap_or(false);
        let seed = seed.unwrap_or(4815162342);

        // 2) Reconstruct CV object from train_data + cv_folds_ids + collections
        let cv = gpredomics::cv::CV::reconstruct(
            &self.intern.train_data,
            self.intern.cv_folds_ids.as_ref().unwrap().clone(),
            self.intern.collections.clone(),
        ).map_err(r_error)?;

        // 3) Native CV importance (agrégé across folds)
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut ic: ImportanceCollection = cv.compute_cv_oob_feature_importance(
            &self.param_arc,
            permutations,
            &mut rng,
            &agg,
            scaled,
            /* cascade */ false,
        ).map_err(r_error)?;

        // 4) Filtre de sécurité : scope = Collection, type = MDA
        ic = ic.filter(
            Some(ImportanceScope::Collection),
            Some(ImportanceType::MDA),
        );

        // 5) Construction du résultat comme pour Population$compute_importance
        let mut feature = Vec::with_capacity(ic.importances.len());
        let mut importance = Vec::with_capacity(ic.importances.len());
        let mut dispersion = Vec::with_capacity(ic.importances.len());
        let mut prevalence = Vec::with_capacity(ic.importances.len());

        for imp in &ic.importances {
            let fname = self
                .intern
                .train_data
                .features
                .get(imp.feature_idx)
                .cloned()
                .unwrap_or_else(|| format!("f{}", imp.feature_idx));
            feature.push(fname);
            importance.push(imp.importance);
            dispersion.push(imp.dispersion);
            prevalence.push(imp.scope_pct);
        }

        if compact {
            // Named vector
            let mut v = Doubles::from_values(importance);
            let _ = v.set_names(feature);
            Ok(v.into_robj())
        } else {
            // data.frame : feature, importance, dispersion, prevalence
            let mut df = data_frame!(
                feature = feature,
                importance = importance,
                dispersion = dispersion,
                prevalence = prevalence
            );

            let imp_label = if scaled { "scaled_mda" } else { "mda" };
            let disp_label = if use_median { "MAD" } else { "sd" };
            let names = Strings::from_values(vec!["feature", imp_label, disp_label, "prevalence"]);
            df.set_attrib(names_symbol(), names)?;
            Ok(df)
        }
    }

    /// @title Compute per-fold CV population-level importance matrix
    /// @name Experiment$compute_cv_importance_matrix
    /// @description
    /// Compute out-of-bag permutation feature importance at the *population* level
    /// for each CV fold separately, then return a matrix with
    /// rows = features (from FBM of each fold) and columns = folds.
    ///
    /// @param n_perm Number of permutations per feature (default: 1000)
    /// @param used_only If TRUE (default), restrict to features used in FBM of folds
    /// @param seed Base seed for RNG (default: 4815162342)
    /// @param aggregation Population-level aggregation: "mean" (default) or "median"
    /// @param should Population-level be scaled
    /// @return Numeric matrix: n_features x n_folds, with rownames = features and
    ///         colnames = Fold_1, Fold_2, ...
    /// @export
    pub fn compute_cv_importance_matrix(
        &self,
        n_perm: Option<i32>,
        used_only: Option<bool>,
        seed: Option<u64>,
        aggregation: Option<String>,
        scaled: Option<bool>,
    ) -> Result<Robj> {
        use std::collections::{HashMap, HashSet};
        use std::sync::Arc;
        use rand_chacha::ChaCha8Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rayon::prelude::*;
        use rayon::ThreadPoolBuilder;

        // ----- Sanity checks -----
        if self.intern.cv_folds_ids.is_none() || self.intern.collections.is_empty() {
            return Err(r_error("No CV folds available. Run with param$cv = TRUE"));
        }

        let permutations = n_perm.unwrap_or(1000).max(1) as usize;
        let used_only = used_only.unwrap_or(true);
        let base_seed = seed.unwrap_or(4815162342);
        let scaled = scaled.unwrap_or(true);

        let agg = match aggregation.as_deref() {
            Some("median") => ImportanceAggregation::median,
            Some("mean") | None => ImportanceAggregation::mean,
            Some(x) => return Err(r_error(format!("Unknown aggregation '{}'", x))),
        };

        let cv_folds_ids = self
            .intern
            .cv_folds_ids
            .as_ref()
            .ok_or_else(|| r_error("No CV fold IDs available in experiment"))?;
        let n_folds = self.intern.collections.len();

        let feature_set: HashSet<usize> = if used_only {
            let mut hs = HashSet::new();
            let fbm_alpha = self.param_arc.cv.cv_best_models_ci_alpha;

            for fold_pops in &self.intern.collections {
                if let Some(last_pop) = fold_pops.last() {
                    let pop = Population {
                        intern: last_pop.clone(),
                        data: Arc::clone(&self.train_data_arc),
                        param: Arc::clone(&self.param_arc),
                    };
                    let fbm = pop.get_fbm(fbm_alpha);

                    for ind in &fbm.intern.individuals {
                        for &f in ind.features.keys() {
                            hs.insert(f);
                        }
                    }
                }
            }
            hs
        } else {
            (0..self.intern.train_data.feature_len).collect()
        };

        if feature_set.is_empty() {
            return Err(r_error("No features found to compute CV importance on."));
        }

        let mut rows: Vec<usize> = feature_set.into_iter().collect();
        rows.sort_unstable();
        let n_feat = rows.len();

        let mut row_index: HashMap<usize, usize> = HashMap::with_capacity(n_feat);
        let mut rownames: Vec<String> = Vec::with_capacity(n_feat);
        for (rpos, fidx) in rows.iter().copied().enumerate() {
            row_index.insert(fidx, rpos);
            let fname = self
                .intern
                .train_data
                .features
                .get(fidx)
                .cloned()
                .unwrap_or_else(|| format!("f{}", fidx));
            rownames.push(fname);
        }

        let threads = resolve_threads()
            .or_else(|| std::thread::available_parallelism().ok().map(|n| n.get()))
            .unwrap_or(1);

        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| r_error(format!("Failed to build thread pool: {}", e)))?;

        let row_index_arc = Arc::new(row_index);
        let _rows_arc = Arc::new(rows);
        let train_data = &self.intern.train_data;
        let _param = &self.param_arc;

        let columns: Vec<Vec<f64>> = pool.install(|| {
            (0..n_folds)
                .into_par_iter()
                .map(|fold_idx| {
                    let (train_ids, _valid_ids) = &cv_folds_ids[fold_idx];
                    let mut train_indices: Vec<usize> = Vec::with_capacity(train_ids.len());

                    for id in train_ids {
                        match train_data.samples.iter().position(|s| s == id) {
                            Some(pos) => train_indices.push(pos),
                            None => {
                                panic!(
                                    "Sample '{}' from CV fold {} not found in training data",
                                    id, fold_idx
                                );
                            }
                        }
                    }

                    let fold_data: GData = train_data.subset(train_indices);

                    let fold_pops = &self.intern.collections[fold_idx];
                    let last_pop = fold_pops
                        .last()
                        .unwrap_or_else(|| {
                            panic!("Fold {} has no populations in collections", fold_idx);
                        })
                        .clone();

                    let mut rng = ChaCha8Rng::seed_from_u64(base_seed);
                    let mut ic: ImportanceCollection = last_pop
                        .compute_pop_oob_feature_importance(
                            &fold_data,
                            permutations,
                            &mut rng,
                            &agg,
                            scaled,
                            /* cascade */ false,
                            Some(fold_idx),
                        );

                    ic = ic.filter(
                        Some(ImportanceScope::Population { id: fold_idx }),
                        Some(ImportanceType::MDA),
                    );

                    let mut col = vec![0.0_f64; n_feat];
                    for imp in &ic.importances {
                        if let Some(&r) = row_index_arc.get(&imp.feature_idx) {
                            col[r] = imp.importance;
                        }
                    }

                    col
                })
                .collect()
        });

        let mut values = vec![0.0_f64; n_feat * n_folds];
        for (j, col) in columns.iter().enumerate() {
            for (r, v) in col.iter().enumerate() {
                values[r + n_feat * j] = *v;
            }
        }

        let colnames: Vec<String> = (0..n_folds)
            .map(|j| format!("Fold_{}", j + 1))
            .collect();

        let mut mat = Doubles::from_values(values).into_robj();
        mat.set_attrib("dim", r!([n_feat as i32, n_folds as i32]))?;

        let rn = Strings::from_values(rownames);
        let cn = Strings::from_values(colnames);
        mat.set_attrib("dimnames", list!(rn, cn))?;

        Ok(mat)
    }


/*
    /// list all individuals at generation #generation
    pub fn get_all_individuals(&self, generation: i32) -> Robj {
        // Create a vector of columns, initializing each column with Vec<i8>
        let mut columns: Vec<Vec<i8>> = vec![Vec::new(); self.train_data.feature_selection.len()];
        let feature_map: HashMap<usize, usize> = self.train_data.feature_selection
                    .iter()
                    .enumerate()
                    .map(|(index, &feature)| (feature, index))
                    .collect();
        println!("features {:?}", self.train_data.feature_selection);
        println!("feature map {:?}",feature_map);

        // Iterate over the rows
        for individual in self.generations[generation as usize].individuals.iter() {
            // Initialize a vector for this row with default values (e.g., 0)
            let mut current_row = vec![0; self.train_data.feature_selection.len()];
            
            
            
            println!("individual {:?}",individual.features);
            // Fill in the row based on the HashMap
            for (feature, &value) in &individual.features {
                current_row[feature_map[feature]] = value;
            }

            // Append the row's values to their respective columns
            for (i, val) in current_row.into_iter().enumerate() {
                columns[i].push(val);
            }
        }

        // Create a named list (similar to DataFrame in R)
        let pairs: Vec<(&str, Robj)> = self.train_data.feature_selection
            .iter()
            .zip(columns.into_iter())
            .map(|(feature_index, column)| (self.train_data.features[*feature_index].as_str(), column.into_robj()))
            .collect();


        // Create the List from pairs
        List::from_pairs(pairs).into_robj()
    }

     */

    /// @title Get number of generations
    /// @description Get the number of generations in the experiment/population (per fold if CV)
    /// @return Integer (no CV) or integer vector (CV)
    pub fn generation_number(&self) -> Robj {
        let n_folds = self.intern.collections.len();
        if n_folds == 0 {
            return r!(0);
        }
        if n_folds == 1 {
            r!(self.intern.collections[0].len() as i32)
        } else {
            let gens: Vec<i32> = self.intern.collections.iter().map(|c| c.len() as i32).collect();
            r!(gens)
        }
    }

    /// @title Get population size for a generation
    /// @param generation An i32 specifying the generation index.
    /// @return Size (number of individuals) of the requested generation
    pub fn population_size(&self, generation: i32) -> Result<i32> {
        if self.intern.collections.len() > 0 {
            Ok(self.intern.collections[0][generation as usize].individuals.len() as i32)
        }
        else {
            Err(r_error("Cannot size any generation in this experiment because it has none"))
        }
    }

    /// @title Load external dataset
    /// @name Experiment$load_data
    /// @description Load an external dataset to evaluate the model
    /// @param x_path Path to the X data file
    /// @param y_path Path to the y data file
    /// @param features_in_columns Logical: TRUE if features are in columns, FALSE if in rows (default: FALSE)
    /// @return Data object containing the loaded data
    pub fn load_data(&self, x_path: String, y_path: String, features_in_columns: bool) -> Result<Data> {
        // Verify files exist before attempting to load
        if !std::path::Path::new(&x_path).exists() {
            return Err(r_error(&format!("X data file not found: {}", x_path)));
        }
        if !std::path::Path::new(&y_path).exists() {
            return Err(r_error(&format!("y data file not found: {}", y_path)));
        }
        
        let mut gdata = GData::new();
        
        // Wrap the load_data call in a panic-catching block
        let load_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            gdata.load_data(&x_path, &y_path, features_in_columns)
        }));
        
        match load_result {
            Ok(Ok(())) => {
                // Successfully loaded
                if !self.intern.train_data.check_compatibility(&gdata) {
                    return Err(r_error("Data not compatible with training data"));
                }
                gdata.set_classes(self.intern.train_data.classes.clone());
                Ok(Data {intern: gdata})
            },
            Ok(Err(e)) => {
                Err(r_error(&format!("Failed to load data: {}", e)))
            },
            Err(panic_err) => {
                let panic_msg = if let Some(s) = panic_err.downcast_ref::<&str>() {
                    s.to_string()
                } else if let Some(s) = panic_err.downcast_ref::<String>() {
                    s.clone()
                } else {
                    "Unknown panic in load_data".to_string()
                };
                Err(r_error(&format!("Panic caught while loading data: {}", panic_msg)))
            }
        }
    }

    /// @title Get the param object
    /// @description Get the param object associated with the experiment
    /// @return Param object containing the experiment parameters
    pub fn get_param(&self) -> Param {
        Param {
            intern: self.intern.parameters.clone()
        }
    }

    /// @title Get the jury object
    /// @description Get the jury object associated with the experiment
    /// @return Jury object containing the jury details
    pub fn get_jury(&self) -> extendr_api::Result<Jury> {
        match &self.intern.others {
            Some(ExperimentMetadata::Jury { jury }) => {
                Ok(Jury {
                    intern: jury.clone(),
                    data: Arc::clone(&self.train_data_arc),
                    param: Arc::clone(&self.param_arc),
                })
            },
            _ =>  Err(r_error("No Jury computed during this experiment"))
        }
    }

    /// @title Load serialized experiment
    /// @name Experiment$load
    /// @description Load a serialized experiment from disk
    /// @param path Path to the experiment file
    /// @return Loaded Experiment object
    pub fn load(path: String) -> Result<Self> {
        if let Ok(experiment) = GExperiment::load_auto(&path) {
            let train_data_arc = Arc::new(experiment.train_data.clone());
            let param_arc = Arc::new(experiment.parameters.clone());    
            Ok(Experiment {
                intern: experiment,
                train_data_arc,
                param_arc,
            })
        }
        else {
            return Err(r_error(&format!("Could not read an experiment from this file: {}",&path)));
        }
    }

    /// @title Save experiment
    /// @name Experiment$save
    /// @description Save an experiment to disk
    /// @param path Path to save the experiment
    pub fn save(&self, path: String) {
        match self.intern.save_auto(&path) {
            Ok(_) => {
                r_message(&format!("Experiment saved in {}", path));
            }
            Err(e) => {
                r_print_error(format!("Could not save an experiment to file {}. Error: {}", path, e));
            }
        }
    }

    /// @title Extract population from experiment
    /// @description Extract population from experiment, optionally specifying generation number
    /// @param generation Optional generation number (0-based). If None, returns final population
    /// @return Population object for the specified generation or final population
    pub fn get_population(&self, generation: Option<i32>) -> Result<Population> {
        match generation {
            None => {
                if let Some(ref final_pop) = self.intern.final_population {
                    Ok(Population {
                        intern: final_pop.clone(),
                        data: Arc::clone(&self.train_data_arc),
                        param: Arc::clone(&self.param_arc),
                    })
                } else {
                    return Err(r_error("No final population available in this experiment. Run the algorithm first."));
                }
            },
            Some(gen_idx) => {
                let gen_idx = gen_idx as usize;

                if !self.intern.collections.is_empty() {
                    if let Some(first_collection) = self.intern.collections.get(0) {
                        if let Some(population) = first_collection.get(gen_idx) {
                            return Ok(Population {
                                intern: population.clone(),
                                data: Arc::clone(&self.train_data_arc),
                                param: Arc::clone(&self.param_arc),
                            });
                        }
                    }
                    return Err(r_error(&format!("Generation {} not found in CV collections", gen_idx)));
                }
                
                return Err(r_error("Generation history not available. Only final population is stored in non-CV mode."));
            }
        }
    }

    /// @title Experiment memory address
    /// @name Experiment$address
    /// @description Get memory address of this Experiment object
    /// @return String representing the memory address
    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern as *const _)
    }
    
    /// @title Experiment print summary
    /// @name Experiment$print
    /// @description Get a formatted summary string for this Experiment
    /// @return String representing the Experiment summary
    pub fn print(&self) -> String {
        let addr = self.address();
        
        // === ALGORITHM INFO ===
        let algorithm = &self.intern.parameters.general.algo;
        let version = &self.intern.gpredomics_version;
        
        // === EXPERIMENT METRICS ===
        let total_generations = if !self.intern.collections.is_empty() {
            self.intern.collections[0].len()
        } else { 0 };
        
        let final_population_size = self.intern.final_population
            .as_ref()
            .map(|p| p.individuals.len())
            .unwrap_or(0);
        
        let best_auc = self.intern.final_population
            .as_ref()
            .and_then(|p| p.individuals.iter()
                .map(|ind| ind.auc)
                .fold(None, |acc, x| Some(acc.map_or(x, |y| x.max(y)))))
            .unwrap_or(0.0);
        
        // === CV INFO ===
        let cv_info = if let Some(ref folds) = self.intern.cv_folds_ids {
            format!("{} folds", folds.len())
        } else {
            "No CV".to_string()
        };
        
        // === EXECUTION TIME ===
        let exec_time = self.intern.execution_time;
        
        format!(
            "@{} \n\
            Experiment: {}:v{}\n\
             Performance Metrics\n\
             \tBest AUC  : {:.3}\n\
             \tExec Time : {:.1}s\n\
             \tTimestamp : {}\n\
             Algorithm Configuration\n\
             \tMethod    : {}\n\
             \tCV Setup  : {}\n\
             \tThreads   : {}\n\
             \tGPU       : {}\n\
             Results\n\
             \tGenerations    : {}\n\
             \tFinal Pop Size : {}\n\
             \tData Samples   : {}",
            addr, algorithm, version,
            best_auc,
            exec_time,
            self.intern.timestamp,
            algorithm,
            cv_info,
            self.intern.parameters.general.thread_number,
            if self.intern.parameters.general.gpu { "Enabled" } else { "Disabled" },
            total_generations,
            final_population_size,
            self.intern.train_data.sample_len
        )
    }

    /// Get training history with scope option: best, fbm, or all population
    /// @title Get Training History with Scope
    /// @name Experiment$get_history
    /// @param data Data object to evaluate performance on.
    /// @param scope Character: "best", "fbm", "top5", or "all" (default "best").
    /// @return A data.frame with columns: Fold, Generation, AUC, Fit, Sensitivity, Specificity,
    ///   Threshold, Rejection_Rate, MCC, F1, NPV, PPV, GMeans, Size_Mean, Size_Best.
    ///   If no CV: Fold column will be 1 for all rows.
    ///   If CV: One row per (Fold, Generation) combination.
    /// @export
    pub fn get_history(&self, data: &Data, scope: &str) -> Result<Robj> {
        if self.intern.collections.is_empty() {
            return Err(r_error("No generations available. Run fit() first."));
        }
        
        let n_folds = self.intern.collections.len();
        let estimated_rows = n_folds * self.intern.collections[0].len();
        
        let mut folds = Vec::with_capacity(estimated_rows);
        let mut generations = Vec::with_capacity(estimated_rows);
        let mut aucs = Vec::with_capacity(estimated_rows);
        let mut fits = Vec::with_capacity(estimated_rows);
        let mut sensitivities = Vec::with_capacity(estimated_rows);
        let mut specificities = Vec::with_capacity(estimated_rows);
        let mut thresholds = Vec::with_capacity(estimated_rows);
        let mut rejection_rates = Vec::with_capacity(estimated_rows);
        let mut k = Vec::with_capacity(estimated_rows);
        let mut mccs = Vec::with_capacity(estimated_rows);
        let mut f1s = Vec::with_capacity(estimated_rows);
        let mut npvs = Vec::with_capacity(estimated_rows);
        let mut ppvs = Vec::with_capacity(estimated_rows);
        let mut gmeans_vec = Vec::with_capacity(estimated_rows);
        
        for (fold_idx, fold_pops) in self.intern.collections.iter().enumerate() {
            let fold_number = (fold_idx + 1) as i32;
            
            for (gen_idx, pop) in fold_pops.iter().enumerate() {
                folds.push(fold_number);
                generations.push((gen_idx + 1) as i32);
                
                if pop.individuals.is_empty() {
                    aucs.push(f64::NAN);
                    fits.push(f64::NAN);
                    sensitivities.push(f64::NAN);
                    specificities.push(f64::NAN);
                    thresholds.push(f64::NAN);
                    rejection_rates.push(f64::NAN);
                    k.push(0.0);
                    mccs.push(f64::NAN);
                    f1s.push(f64::NAN);
                    npvs.push(f64::NAN);
                    ppvs.push(f64::NAN);
                    gmeans_vec.push(f64::NAN);
                    continue;
                }
                
                // Choose individuals to evaluate based on scope
                let mut hist_pop: GPopulation = match scope {
                    "best" => GPopulation { individuals: vec![ pop.individuals.first().unwrap().clone() ] },
                    "fbm" => {
                        let alpha = 0.05;
                        pop.select_best_population(alpha)
                    },
                    "top5" => {
                        pop.select_first_pct(5.0).0
                    },
                    "all" => pop.clone(),
                    _ => return Err(r_error(format!("Unknown scope '{}', should be 'best', 'fbm' or 'all'", scope))),
                };

                hist_pop.fit(&data.intern, &mut None, &None, &None, &self.intern.parameters);
                
                // Initialize accumulators
                let mut auc_sum = 0.0;
                let mut fit_sum = 0.0;
                let mut sensitivity_sum = 0.0;
                let mut specificity_sum = 0.0;
                let mut threshold_sum = 0.0;
                let mut rejection_sum = 0.0;
                let mut k_sum = 0.0;
                
                // Additional metrics accumulators
                let mut mcc_sum = 0.0;
                let mut f1_sum = 0.0;
                let mut npv_sum = 0.0;
                let mut ppv_sum = 0.0;
                let mut gmeans_sum = 0.0;
                
                let n = hist_pop.individuals.len() as f64;
                if n == 0.0 {
                    return Err(r_error("Selected population scope is empty"));
                }
                
                for ind in &hist_pop.individuals {
                    auc_sum += ind.auc;
                    fit_sum += ind.fit;
                    sensitivity_sum += ind.sensitivity;
                    specificity_sum += ind.specificity;
                    threshold_sum += ind.threshold;
                    k_sum += ind.k as f64;

                    let rej = ind
                        .threshold_ci
                        .as_ref()
                        .map(|ci| ci.rejection_rate)
                        .unwrap_or(0.0);
                    rejection_sum += rej;
                    
                    mcc_sum += ind.metrics.mcc.unwrap_or(0.0);
                    f1_sum += ind.metrics.f1_score.unwrap_or(0.0);
                    npv_sum += ind.metrics.npv.unwrap_or(0.0);
                    ppv_sum += ind.metrics.ppv.unwrap_or(0.0);
                    gmeans_sum += ind.metrics.g_mean.unwrap_or(0.0);
                }
                
                aucs.push(auc_sum / n);
                fits.push(fit_sum / n);
                sensitivities.push(sensitivity_sum / n);
                specificities.push(specificity_sum / n);
                thresholds.push(threshold_sum / n);
                rejection_rates.push(rejection_sum / n);
                k.push(k_sum / n);
                mccs.push(mcc_sum / n);
                f1s.push(f1_sum / n);
                npvs.push(npv_sum / n);
                ppvs.push(ppv_sum / n);
                gmeans_vec.push(gmeans_sum / n);
            }
        }
    
        let df = if n_folds > 1 {
             data_frame!(
                Fold = folds,
                Generation = generations,
                AUC = aucs,
                Fit = fits,
                Sensitivity = sensitivities,
                Specificity = specificities,
                Threshold = thresholds,
                Rejection_Rate = rejection_rates,
                MCC = mccs,
                F1 = f1s,
                NPV = npvs,
                PPV = ppvs,
                GMean = gmeans_vec,
                k = k)
            } else {
                data_frame!(
                    Generation = generations,
                    AUC = aucs,
                    Fit = fits,
                    Sensitivity = sensitivities,
                    Specificity = specificities,
                    Threshold = thresholds,
                    Rejection_Rate = rejection_rates,
                    MCC = mccs,
                    F1 = f1s,
                    NPV = npvs,
                    PPV = ppvs,
                    GMean = gmeans_vec,
                    k = k)
            };
        
        Ok(df.into())
    }
}

///////////////////////////////////////////////////////////////
/// Population object
///////////////////////////////////////////////////////////////

#[extendr]
#[derive(Debug, Clone)]
pub struct Population {
    intern: GPopulation,
    data: Arc<GData>,
    param: Arc<GParam>,
}

impl Population {
    /// Internal constructor - creates Population with shared data and param references.
    /// 
    /// This method is used internally to wrap a GPopulation with its associated
    /// data and parameters for use in the R interface.
    /// 
    /// # Arguments
    /// * `intern` - The internal GPopulation object
    /// * `data` - Shared reference to the data
    /// * `param` - Shared reference to algorithm parameters
    pub(crate) fn new(intern: GPopulation, data: Arc<GData>, param: Arc<GParam> ) -> Self {
        Self { intern, data, param }
    }
}

/// @title Population
/// @name Population
/// @description gpredomics Population object
/// @details 
/// Population represents a collection of Individuals.
/// It provides methods to filter, analyze, and manipulate sets of individuals, as well as compute
/// aggregate predictions and feature importances across the population.
/// 
/// @section Methods:
/// \describe{
///   \item{\code{get()}}{Get the population associated with the experiment. Returns R object representing the Population.}
///   \item{\code{display_feature_prevalence(data, nb_features)}}{Display the prevalence of features in the population. Parameters: data (Data object), nb_features (number of top features to display).}
///   \item{\code{predict_score_matrix(data)}}{Predict all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts). Parameter: data (Data object to predict on). Returns dataframe with predicted scores.}
///   \item{\code{predict_class_matrix(data)}}{Predict classes for all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts, Values = predicted classes 0 or 1). Parameter: data (Data object to predict on). Returns dataframe with predicted classes.}
///   \item{\code{filter_by_auc(min_auc)}}{Filter population by AUC threshold. Parameter: min_auc (minimum AUC threshold). Returns filtered Population object.}
///   \item{\code{filter_by_fit(min_fit)}}{Filter population by fitness threshold. Parameter: min_fit (minimum fit threshold). Returns filtered Population object.}
///   \item{\code{filter_by_diversity(min_diversity_pct, by_niche)}}{Filter population by diversity using Jaccard dissimilarity. Parameters: min_diversity_pct (minimum diversity percentage 0-100), by_niche (whether to compute diversity within niches). Returns filtered Population object.}
///   \item{\code{filter_by_sensitivity(min_sensitivity)}}{Filter population by sensitivity threshold. Parameter: min_sensitivity (minimum sensitivity threshold). Returns filtered Population object.}
///   \item{\code{filter_by_specificity(min_specificity)}}{Filter population by specificity threshold. Parameter: min_specificity (minimum specificity threshold). Returns filtered Population object.}
///   \item{\code{filter_by_mask(mask)}}{Filter population using a logical vector (1/0). Parameter: mask (integer vector 1/0 indicating which individuals to keep). Returns filtered Population object.}
///   \item{\code{filter_by_k(min_k, max_k)}}{Filter population by number of features (k). Parameters: min_k (minimum number of features), max_k (maximum number of features). Returns filtered Population object.}
///   \item{\code{get_fbm(alpha)}}{Get Family of Best Models (FBM) using confidence interval selection. This method selects models with performance statistically equivalent to the best model. Parameters: alpha (confidence level, default 0.05 for 95% confidence). If FBM selection fails, alpha% to keep). Returns Population object containing the FBM.}
///   \item{\code{get_first_pct(pct)}}{Get first percentage of individuals sorted by fitness. Parameter: pct (percentage 0-100). Returns Population object containing the selected individuals.
///   \item{\code{fit(data, param)}}{Compute fitness metrics for all individuals on new data or parameters and sort it. Parameters: data (new Data object to fit on), param (Param object containing fit function and penalties).}
///   \item{\code{prune_by_threshold(threshold, n_perm, seed, min_k)}}{Prune all individuals by importance threshold. Parameters: threshold (importance threshold), n_perm (number of permutations, default 100), seed (base seed for RNG, default 4815162342), min_k (minimum features to keep, default 1). Returns new Population with pruned individuals.}
///   \item{\code{prune_by_quantile(quantile, eps, n_perm, seed, min_k)}}{Prune all individuals by importance quantile. Parameters: quantile (quantile value 0-1), eps (epsilon value, default 0.0), n_perm (number of permutations, default 100), seed (base seed for RNG, default 4815162342), min_k (minimum features to keep, default 1). Returns new Population with pruned individuals.}
///   \item{\code{address()}}{Get memory address of this Population object. Returns string representing the memory address.}
///   \item{\code{get_individual(index)}}{Get an Individual of a population by index. Parameter: index (index of the individual to retrieve). Returns Individual object at the specified index.}
///   \item{\code{print_as_gpredomics()}}{Print the Population in gpredomics style to R console.}
///   \item{\code{print()}}{Get comprehensive print information about the population. Returns string representing the Population summary.}
///   \item{\code{from_individuals(individuals)}}{Create a Population from a vector or list of R Individual objects. Parameter: individuals (R vector or list of Individual objects, must have at least one). Returns Population object created from the individuals.}
///   \item{\code{extend(other)}}{Extend this population with another population. Parameter: other (another Population object to add).}
///   \item{\code{add_individuals(individuals)}}{Add individuals from a vector or list to this population. Parameter: individuals (R vector or list of Individual objects).}
///   \item{\code{compute_importance(data, n_perm, aggregation, scaled, seed, compact)}}{Compute full Population-level MDA importances for this Population on given Data. Parameters: data (Data object to compute importances on), n_perm (number of permutations, default 1000), aggregation (aggregation method: "mean" (default) or "median"), scaled (whether to scale importances, default TRUE), seed (random seed for reproducibility, default 4815162342), compact (whether to return a compact vector (TRUE) or full data.frame (FALSE, default)). Returns DataFrame with columns: feature, importance, dispersion, prevalence.}
///   \item{\code{compute_importance_matrix(data, n_perm, used_only, seed)}}{Compute full Population-level MDA importance matrix for this Population on given Data. Parameters: data (Data object to compute importances on), n_perm (number of permutations, default 1000), used_only (whether to compute importances only for features used in the population, default TRUE), seed (random seed for reproducibility, default 4815162342). Returns Matrix (data.frame) with rows = features, columns = individuals.}
/// }
/// @export
#[extendr]
impl Population {

    /// @title Get population
    /// @description Get the population associated with the experiment
    /// @return R object representing the Population
    /// @export
    pub fn get(&self) -> Robj {
        let individuals = self.intern.individuals
            .iter()
            .map(|gi: &GIndividual| Individual::new(gi, 
                                    Arc::clone(&self.data),
                                    Arc::clone(&self.param)).get())
            .collect::<Vec<Robj>>();
        List::from_pairs(vec![
            ("individuals", Robj::from(individuals)),
            ("rust", self.clone().into_robj())
        ]).into_robj()
    }

    /// @title Display feature prevalence
    /// @name Population$display_feature_prevalence
    /// @description Display the prevalence of features in the population
    /// @param data The Data object to compute feature prevalence on
    /// @param nb_features Number of top features to display
    /// @export
    pub fn display_feature_prevalence(&self, data: &Data, nb_features: usize) {
        println!("{}",self.intern.display_feature_prevalence(&data.intern, nb_features));
    }

    /// @title Predict scores matrix
    /// @name Population$predict_score_matrix
    /// @description Predict all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts)
    /// @param data The Data object to predict on
    /// @return Dataframe with predicted scores
    /// @export
    pub fn predict_score_matrix(&self, data: &Data) -> Result<Robj> {
        let num_samples = data.intern.sample_len;
        let num_individuals = self.intern.individuals.len();
        
        if num_individuals == 0 {
            return Err(r_error("Population is empty - cannot predict"));
        }
        
        // Create matrix of predictions: rows=samples, cols=individuals
        let mut prediction_matrix: Vec<Vec<f64>> = vec![Vec::with_capacity(num_samples); num_individuals];
        
        // Fill predictions for each individual
        for (individual_idx, individual) in self.intern.individuals.iter().enumerate() {
            let predictions = individual.evaluate(&data.intern);
            prediction_matrix[individual_idx] = predictions;
        }
        
        // Create column names (Individual_1, Individual_2, ...)
        let column_names: Vec<String> = (1..=num_individuals)
            .map(|i| format!("Individual_{}", i))
            .collect();
        
        // Convert to R dataframe
        let mut df_list = List::from_values(
            prediction_matrix.iter().map(|col| col.into_robj())
        );
        
        let _ = df_list.set_names(&column_names);
        df_list.set_attrib(row_names_symbol(), Robj::from(&data.intern.samples)).ok();
        df_list.set_class(&["data.frame"]).ok();
        Ok(df_list.into())
    }
    
    /// @title Predict classes matrix
    /// @name Population$predict_class_matrix
    /// @description Predict classes for all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts, Values = predicted classes 0 or 1)
    /// @param data The Data object to predict on
    /// @return Dataframe with predicted classes
    /// @export  
    pub fn predict_class_matrix(&self, data: &Data) -> Result<Robj> {
        let num_samples = data.intern.sample_len;
        let num_individuals = self.intern.individuals.len();
        
        if num_individuals == 0 {
            return Err(r_error("Population is empty - cannot predict"));
        }
        
        // Create matrix of class predictions: rows=samples, cols=individuals  
        let mut prediction_matrix: Vec<Vec<i32>> = vec![Vec::with_capacity(num_samples); num_individuals];
        
        // Fill class predictions for each individual
        for (individual_idx, individual) in self.intern.individuals.iter().enumerate() {
            let scores = individual.evaluate(&data.intern);
            let classes: Vec<i32> = scores.iter()
                .map(|&score| if score >= individual.threshold { 1 } else { 0 })
                .collect();
            prediction_matrix[individual_idx] = classes;
        }
        
        // Create column names (Individual_1, Individual_2, ...)
        let column_names: Vec<String> = (1..=num_individuals)
            .map(|i| format!("Individual_{}", i))
            .collect();
        
        // Convert to R dataframe
        let mut df_list = List::from_values(
            prediction_matrix.iter().map(|col| col.into_robj())
        );
        
        let _ = df_list.set_names(&column_names);
        df_list.set_attrib(row_names_symbol(), Robj::from(&data.intern.samples)).ok();
        df_list.set_class(&["data.frame"]).ok();
        Ok(df_list.into())
    }

    /// @title Filter population by AUC threshold
    /// @param min_auc Minimum AUC threshold
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_auc(&self, min_auc: f64) -> Population {
        let mut gpopulation = GPopulation::new();
        
        for individual in &self.intern.individuals {
            if individual.auc >= min_auc {
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Population::new(gpopulation, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Filter population by fitness threshold
    /// @param min_fit Minimum fit threshold
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_fit(&self, min_fit: f64) -> Population {
        let mut gpopulation = GPopulation::new();
        
        for individual in &self.intern.individuals {
            if individual.fit >= min_fit {
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Population::new(gpopulation, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Filter population by diversity using Jaccard dissimilarity
    /// @param min_diversity_pct Minimum diversity percentage (0-100)
    /// @param by_niche Whether to compute diversity within niches
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_diversity(&self, min_diversity_pct: f64, by_niche: bool) -> Population {
        let filtered_pop = self.intern.filter_by_signed_jaccard_dissimilarity(min_diversity_pct, by_niche);
        
        Population::new(filtered_pop, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Filter population by sensitivity threshold
    /// @param min_sensitivity Minimum sensitivity threshold
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_sensitivity(&self, min_sensitivity: f64) -> Population {
        let mut gpopulation = GPopulation::new();
        
        for individual in &self.intern.individuals {
            if individual.sensitivity >= min_sensitivity {
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Population::new(gpopulation, Arc::clone(&self.data), Arc::clone(&self.param))
    }
    
    /// @title Filter population by specificity threshold
    /// @param min_specificity Minimum specificity threshold
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_specificity(&self, min_specificity: f64) -> Population {
        let mut gpopulation = GPopulation::new();
        
        for individual in &self.intern.individuals {
            if individual.specificity >= min_specificity {
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Population::new(gpopulation, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Filter population by mask
    /// @name Population$filter_by_mask
    /// @description Filter population using a logical vector (1/0)
    /// @param mask Integer vector (1/0) indicating which individuals to keep
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_mask(&self, mask: Vec<i32>) -> Result<Population> {
        if mask.len() != self.intern.individuals.len() {
            Err(r_error(&format!("Mask length ({}) must match population size ({})", 
                   mask.len(), self.intern.individuals.len())))?;
        }
        
        let mut gpopulation = GPopulation::new();
        
        for (individual, keep) in self.intern.individuals.iter().zip(mask.iter()) {
            if *keep != 0 { 
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Ok(Population::new(gpopulation, Arc::clone(&self.data), Arc::clone(&self.param)))
    }

    /// @title Filter population by number of features (k)
    /// @param min_k Minimum number of features
    /// @param max_k Maximum number of features
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_k(&self, min_k: i32, max_k: i32) -> Population {
        let min_k_usize = min_k as usize;
        let max_k_usize = max_k as usize;
        let mut gpopulation = GPopulation::new();
        
        for individual in &self.intern.individuals {
            if individual.k >= min_k_usize || individual.k <= max_k_usize {
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Population::new(gpopulation, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Get Family of Best Models (FBM)
    /// @name Population$get_fbm
    /// @description Get Family of Best Models (FBM) using confidence interval selection. This method selects models with performance statistically equivalent to the best model
    /// @param alpha Confidence level (default: 0.05 for 95% confidence)
    /// @param min_pct_fallback If FBM selection fails, minimum percentage to keep (default: 5.0)
    /// @return Population object containing the FBM
    /// @export
    pub fn get_fbm(&self, alpha: f64) -> Population {
        let mut sorted_population = self.clone();
        sorted_population.intern = sorted_population.intern.sort();
        
        let fbm_population = sorted_population.intern.select_best_population(alpha);
        
        if fbm_population.individuals.is_empty() {
            r_warning(format!("FBM selection failed or returned empty population, using {}% fallback", alpha));
            let (fallback_pop, _) = sorted_population.intern.select_first_pct(alpha);
            Population::new(fallback_pop, Arc::clone(&self.data), Arc::clone(&self.param))
        } else {
            Population::new(fbm_population, Arc::clone(&self.data), Arc::clone(&self.param))
        }
    }

    /// @title Get first percentage of population
    /// @name Population$get_first_pct
    /// @description Get the first percentage of the population based on fitness
    /// @param pct Percentage of the population to select (0-100)
    /// @return Population object containing the selected individuals
    /// @export
    pub fn get_first_pct(&self, pct: f64) -> Population {
        let (first_pct_pop, _) = self.intern.select_first_pct(pct);
        Population::new(first_pct_pop, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Recompute fitness metrics
    /// @name Population$fit
    /// @description Recompute fitness metrics for all individuals on new data using param settings
    /// @param data New Data object to fit on
    /// @param param Param object containing fit function and penalties
    /// @export
    pub fn fit(&self, data: &Data, param: &Param) -> Population {
        let mut pop = self.intern.clone();

        pop.fit(&data.intern, &mut None, &None, &None, &param.intern);
        pop = pop.sort();

        Population {
            intern: pop,
            data: Arc::new(data.intern.clone()),
            param: Arc::new(param.intern.clone()),
        }
    }

    /// @title Prune all individuals by importance threshold
    /// @name Population$prune_by_threshold
    /// @description Prune (remove) features from all individuals in the population based on OOB permutation importance.
    /// This is performed in parallel across all individuals. Features with importance below the specified threshold
    /// are removed from each individual, while ensuring at least `min_k` features remain.
    /// Metrics are recomputed for all individuals after pruning.
    /// Returns a new Population with pruned individuals.
    /// @param threshold Importance threshold - drops features with importance < threshold
    /// @param n_perm Number of permutations to perform for each individual (default 100)
    /// @param seed Base seed for random number generation (default 4815162342)
    /// @param min_k Minimum number of features to keep in each individual (default 1)
    /// @return A new Population with all individuals pruned and metrics recomputed
    /// @section Methods:
    /// This method creates a new Population where all individuals have low-importance features removed.
    /// It uses parallel processing (Rayon) for efficient computation across the population.
    /// Each individual uses a deterministic seed derived from the base seed and individual hash.
    /// After pruning, all metrics are recomputed using the fit() method.
    /// @export
    pub fn prune_by_threshold(
        &self,
        threshold: f64,
        n_perm: Option<i32>,
        seed: Option<u64>,
        min_k: Option<i32>,
    ) -> Population {
        let permutations = n_perm.unwrap_or(100).max(1) as usize;
        let base_rng_seed = seed.unwrap_or(4815162342);
        let min_features = min_k.unwrap_or(1).max(1) as usize;
        
        // Clone the population
        let mut pruned_pop = self.intern.clone();
        
        pruned_pop.prune_all_by_importance(
            &self.data,
            permutations,
            base_rng_seed,
            Some(threshold),
            None,
            min_features,
        );
        
        // Recompute metrics after pruning
        pruned_pop.fit(&self.data, &mut None, &None, &None, &self.param);
        pruned_pop = pruned_pop.sort();
        
        Population::new(pruned_pop, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Prune all individuals by importance quantile
    /// @name Population$prune_by_quantile
    /// @description Prune (remove) features from all individuals in the population based on OOB permutation importance.
    /// This is performed in parallel across all individuals. Features with importance below quantile(q) - eps
    /// are removed from each individual, while ensuring at least `min_k` features remain.
    /// Metrics are recomputed for all individuals after pruning.
    /// Returns a new Population with pruned individuals.
    /// @param quantile Quantile value between 0 and 1 (e.g., 0.25 for 25th percentile)
    /// @param eps Optional epsilon value to subtract from quantile (default 0.0)
    /// @param n_perm Number of permutations to perform for each individual (default 100)
    /// @param seed Base seed for random number generation (default 4815162342)
    /// @param min_k Minimum number of features to keep in each individual (default 1)
    /// @return A new Population with all individuals pruned and metrics recomputed
    /// @section Methods:
    /// This method creates a new Population where all individuals have low-importance features removed.
    /// It uses parallel processing (Rayon) for efficient computation across the population.
    /// Each individual uses a deterministic seed derived from the base seed and individual hash.
    /// The cutoff for each individual is computed as: quantile(importance, q) - eps
    /// After pruning, all metrics are recomputed using the fit() method.
    /// @export
    pub fn prune_by_quantile(
        &self,
        quantile: f64,
        eps: Option<f64>,
        n_perm: Option<i32>,
        seed: Option<u64>,
        min_k: Option<i32>,
    ) -> Population {
        let permutations = n_perm.unwrap_or(100).max(1) as usize;
        let base_rng_seed = seed.unwrap_or(4815162342);
        let min_features = min_k.unwrap_or(1).max(1) as usize;
        let epsilon = eps.unwrap_or(0.0);
        
        // Clone the population
        let mut pruned_pop = self.intern.clone();

        pruned_pop.prune_all_by_importance(
            &self.data,
            permutations,
            base_rng_seed,
            None,
            Some((quantile, epsilon)),
            min_features,
        );
        
        // Recompute metrics after pruning
        pruned_pop.fit(&self.data, &mut None, &None, &None, &self.param);
        pruned_pop = pruned_pop.sort();
        
        Population::new(pruned_pop, Arc::clone(&self.data), Arc::clone(&self.param))
    }

    /// @title Get Population memory address
    /// @name Population$address
    /// @description Get memory address of this Population object
    /// @return String representing the memory address
    /// @export
    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern as *const _)
    }

    /// @title Get individual from population
    /// @name Population$get_individual
    /// @description Get an Individual of a population by index
    /// @param index Index of the individual to retrieve
    /// @return Individual object at the specified index
    /// @export
    pub fn get_individual(&self, index:i32) -> Result<Individual> {
        if index <= 0 || index as usize > self.intern.individuals.len() {
            return Err(r_error(&format!("Index {} out of bounds for population of size {}", index, self.intern.individuals.len())));
        }

        if self.intern.individuals.len() > 0 {
            Ok(Individual::new(&self.intern.individuals[index as usize - 1], Arc::clone(&self.data), Arc::clone(&self.param)))
        } else {
            Err(r_error("Cannot extract an individual from an empty Population"))
        }
    }

    // @title Print as Gpredomics style
    // @name Population$print_as_gpredomics
    /// @title Print Population (formatted)
    /// @name Population$print
    /// @description Print the Population in gpredomics style to R console.
    /// @export
    pub fn print_as_gpredomics(&self) {
        let addr = self.address();
        rprintln!("@{}", addr);
        let s = &self.intern.clone().display(&self.data, None, &self.param);
        for line in s.lines() {
            rprintln!("{}", line);
        }
    }
    
    /// @title Get comprehensive print information
    /// @name Population$print
    /// @description Get comprehensive print information about the population
    /// @return String representing the Population summary
    /// @export
    pub fn print(&self) -> String {
        let count = self.intern.individuals.len();
        let addr = self.address();
        
        if count == 0 {
            return format!("Population@{}: Empty (0 individuals)", addr);
        }
        
        // === PERFORMANCE METRICS ===
        let aucs: Vec<f64> = self.intern.individuals.iter().map(|ind| ind.auc).collect();
        let fits: Vec<f64> = self.intern.individuals.iter().map(|ind| ind.fit).collect();
        let sensitivities: Vec<f64> = self.intern.individuals.iter().map(|ind| ind.sensitivity).collect();
        let specificities: Vec<f64> = self.intern.individuals.iter().map(|ind| ind.specificity).collect();
        let accuracies: Vec<f64> = self.intern.individuals.iter().map(|ind| ind.accuracy).collect();

        // === STRUCTURAL METRICS ===
        let ks: Vec<usize> = self.intern.individuals.iter().map(|ind| ind.k).collect();
        let epochs: Vec<usize> = self.intern.individuals.iter().map(|ind| ind.epoch).collect();
        
        // === LANGUAGE/data_type ANALYSIS ===
        let mut language_counts = std::collections::HashMap::new();
        let mut data_type_counts = std::collections::HashMap::new();
        
        for individual in &self.intern.individuals {
            let lang = match individual.language {
                BINARY_LANG => "Binary",
                TERNARY_LANG => "Ternary",
                POW2_LANG => "Pow2",
                RATIO_LANG => "Ratio",
                MCMC_GENERIC_LANG => "MCMC",
                _ => "Unknown",
            };
            *language_counts.entry(lang).or_insert(0) += 1;
            
            let dtype = match individual.data_type {
                RAW_TYPE => "Raw",
                PREVALENCE_TYPE => "Prevalence",
                LOG_TYPE => "Log",
                _ => "Unknown",
            };
            *data_type_counts.entry(dtype).or_insert(0) += 1;
        }
        
        // === STATISTICAL CALCULATIONS ===
        fn stats(values: &[f64]) -> (f64, f64, f64, f64) {
            if values.is_empty() { return (0.0, 0.0, 0.0, 0.0); }
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            let min = values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let max = values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let variance = if values.len() > 1 {
                values.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64
            } else { 0.0 };
            let std = variance.sqrt();
            (mean, min, max, std)
        }
        
        fn stats_usize(values: &[usize]) -> (f64, usize, usize, f64) {
            if values.is_empty() { return (0.0, 0, 0, 0.0); }
            let mean = values.iter().sum::<usize>() as f64 / values.len() as f64;
            let min = *values.iter().min().unwrap();
            let max = *values.iter().max().unwrap();
            let variance = if values.len() > 1 {
                values.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64
            } else { 0.0 };
            let std = variance.sqrt();
            (mean, min, max, std)
        }
        
        let (auc_mean, auc_min, auc_max, auc_std) = stats(&aucs);
        let (fit_mean, fit_min, fit_max, fit_std) = stats(&fits);
        let (sens_mean, sens_min, sens_max, sens_std) = stats(&sensitivities);
        let (spec_mean, spec_min, spec_max, spec_std) = stats(&specificities);
        let (acc_mean, acc_min, acc_max, acc_std) = stats(&accuracies);
        let (k_mean, k_min, k_max, k_std) = stats_usize(&ks);
        let (epoch_mean, epoch_min, epoch_max, epoch_std) = stats_usize(&epochs);
        
        // === DIVERSITY METRICS ===
        let unique_hashes: std::collections::HashSet<_> = 
            self.intern.individuals.iter().map(|ind| ind.hash).collect();
        let diversity_pct = (unique_hashes.len() as f64 / count as f64) * 100.0;
        
        // === LANGUAGE/data_type print ===
        let lang_print: Vec<String> = language_counts.iter()
            .map(|(lang, count)| format!("{}: {}", lang, count))
            .collect();
        let dtype_print: Vec<String> = data_type_counts.iter()
            .map(|(dtype, count)| format!("{}: {}", dtype, count))
            .collect();
        
        // === FORMAT OUTPUT ===
        format!(
            "@{}\n\
            Population: {} individuals ({}% unique)\n\
             Performance Metrics\n\
             \tAUC       : {:.3}±{:.3} [{:.3}, {:.3}]\n\
             \tFitness   : {:.3}±{:.3} [{:.3}, {:.3}]\n\
             \tAccuracy  : {:.3}±{:.3} [{:.3}, {:.3}]\n\
             \tSensitivity: {:.3}±{:.3} [{:.3}, {:.3}]\n\
             \tSpecificity: {:.3}±{:.3} [{:.3}, {:.3}]\n\
             Structural Properties\n\
             \tFeatures  : {:.1}±{:.1} [{}, {}]\n\
             \tGeneration: {:.1}±{:.1} [{}, {}]\n\
             Composition\n\
             \tLanguages : {}\n\
             \tData Types: {}\n",

            addr, count, diversity_pct as u32,
            auc_mean, auc_std, auc_min, auc_max,
            fit_mean, fit_std, fit_min, fit_max,
            acc_mean, acc_std, acc_min, acc_max,
            sens_mean, sens_std, sens_min, sens_max,
            spec_mean, spec_std, spec_min, spec_max,
            k_mean, k_std, k_min, k_max,
            epoch_mean, epoch_std, epoch_min, epoch_max,
            lang_print.join(", "),
            dtype_print.join(", ")
        )
    }

    /// @title Create Population from individuals
    /// @name Population$from_individuals
    /// @description Create a Population from a vector or list of R Individual objects
    /// @param individuals R vector or list of Individual objects (must have at least one)
    /// @return Population object created from the individuals
    /// @export
    pub fn from_individuals(individuals: Robj) -> Result<Population> {
        let mut gpopulation = GPopulation::new();
        let mut data_arc: Option<Arc<GData>> = None;
        let mut param_arc: Option<Arc<GParam>> = None;

        let mut push_individual = |robj: &Robj| -> Result<()> {
            match <&Individual>::try_from(robj) { 
                Ok(r_individual) => {
                    // Collect or verify Data Arc
                    if let Some(ref existing_data) = data_arc {
                        if !Arc::ptr_eq(existing_data, &r_individual.data) {
                            return Err(r_error("All Individuals must share the same Data object to build a Population"));
                        }
                    } else {
                        data_arc = Some(Arc::clone(&r_individual.data));
                    }

                    // Collect or verify Param Arc
                    if let Some(ref existing_param) = param_arc {
                        if !Arc::ptr_eq(existing_param, &r_individual.param) {
                            return Err(r_error("All Individuals must share the same Param object to build a Population"));
                        }
                    } else {
                        param_arc = Some(Arc::clone(&r_individual.param));
                    }

                    gpopulation.individuals.push(r_individual.intern.clone());
                    Ok(())
                }
                Err(_) => {
                    Err(r_error("Invalid Individual object passed to Population::from_individuals"))
                }
            }
        };
        
        // Process the input (vector or list)
        if let Some(r_list) = individuals.as_list() {
            for individual_robj in r_list.values() {
                push_individual(&individual_robj)?;
            }
        } else if individuals.is_vector() {
            let len = individuals.len();
            for i in 0..len {
                match individuals.index(i) {
                    Ok(individual_robj) => {
                        push_individual(&individual_robj)?;
                    }
                    Err(e) => {
                        return Err(r_error(format!("Cannot access element at index {}: {:?}", i, e)));
                    }
                }
            }
        } else {
            return Err(r_error("Input must be a vector or list of Individual objects"));
        }

        if gpopulation.individuals.is_empty() {
            return Err(r_error("Cannot create Population from empty Individual list"));
        }

        // Ensure both data_arc and param_arc were initialized
        let data = data_arc.ok_or_else(|| r_error("BUG: data_arc not initialized despite non-empty individual list"))?;
        let param = param_arc.ok_or_else(|| r_error("BUG: param_arc not initialized despite non-empty individual list - ensure all Individual objects have valid param references"))?;

        Ok(Population::new(gpopulation, data, param))
    }

    /// @title Extend population
    /// @name Population$extend
    /// @description Extend this population with another population (in-place modification)
    /// @param other Another Population object to add
    /// @export
    pub fn extend(&mut self, other: &Population) {
        self.intern.individuals.extend(other.intern.individuals.clone());
    }
    
    /// @title Add individuals to population
    /// @name Population$add_individuals
    /// @description Add individuals from a vector or list to this population
    /// @param individuals R vector or list of Individual objects
    /// @export
    pub fn add_individuals(&mut self, individuals: Robj) {
        let individuals_pop = Population::from_individuals(individuals);
        
        self.extend(&individuals_pop.expect("Failed to add individuals to population"));
    }

    /// @title Compute Population-level MDA importances
    /// @name Population$compute_importance
    /// @description Compute full Population-level MDA importances for this Population on given Data
    /// @param data Data object to compute importances on
    /// @param n_perm Number of permutations (default: 1000)
    /// @param aggregation Aggregation method: "mean" (default) or "median"
    /// @param scaled Whether to scale importances (default: TRUE)
    /// @param seed Random seed for reproducibility (default: 4815162342)
    /// @param compact Whether to return a compact vector (TRUE) or full data.frame (FALSE, default) 
    /// @return DataFrame with columns: feature, importance, dispersion, prevalence
    /// @export
    pub fn compute_importance(&self, data: &Data, n_perm: Option<i32>, aggregation: Option<String>, scaled: Option<bool>, 
        seed: Option<u64>, compact: Option<bool> ) -> Result<Robj> {
        let permutations = n_perm.unwrap_or(1000).max(0) as usize;
        let use_median = matches!(aggregation.as_deref(), Some("median"));
        let agg = match aggregation.as_deref() {
            Some("median") => ImportanceAggregation::median,
            Some("mean") | None => ImportanceAggregation::mean,
            Some(x) => return Err(r_error(format!("Unknown aggregation '{}'", x))),
        };
        let scaled = scaled.unwrap_or(true);
        let compact = compact.unwrap_or(false);

        let seed = seed.unwrap_or(4815162342);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);

        // Configure thread pool using resolve_threads
        let threads = resolve_threads()
            .or_else(|| std::thread::available_parallelism().ok().map(|n| n.get()))
            .unwrap_or(1);

        let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();

        let mut ic: ImportanceCollection = pool.install(|| {
            self.intern.compute_pop_oob_feature_importance(
                &data.intern,
                permutations,
                &mut rng,
                &agg,
                scaled,
                /* cascade */ false,
                Some(1),
            )
        }); // importance_type = MDA; scope = Population { id: pid } 

        // 4) Safety filters: keep only Population MDA (should already be the case)
        ic = ic.filter(Some(ImportanceScope::Population { id: 0 }), Some(ImportanceType::MDA)); // 

        let mut feature = Vec::with_capacity(ic.importances.len());
        let mut importance = Vec::with_capacity(ic.importances.len());
        let mut dispersion = Vec::with_capacity(ic.importances.len());
        let mut prevalence = Vec::with_capacity(ic.importances.len());

        for imp in &ic.importances {
            let fname = data.intern.features
                .get(imp.feature_idx)
                .cloned()
                .unwrap_or_else(|| format!("f{}", imp.feature_idx));
            feature.push(fname);
            importance.push(imp.importance);
            dispersion.push(imp.dispersion);
            prevalence.push(imp.scope_pct); 
        }

        if compact {
            let mut v = Doubles::from_values(importance);
            let _ = v.set_names(feature);
            return Ok(v.into_robj());
        }

        let mut df = data_frame!(
            feature = feature,
            importance = importance,
            dispersion = dispersion,
            prevalence = prevalence
        ); 

        let imp_label = if scaled { "scaled_mda" } else { "mda" };
        let disp_label = if use_median { "MAD" } else { "sd" };
        let names = Strings::from_values(vec!["feature", imp_label, disp_label, "prevalence"]);
        df.set_attrib(names_symbol(), names)?;
        Ok(df)
    }

    // @title Compute importance matrix
    // @name Population$compute_importance_matrix
    // @description Compute full Population-level MDA importance matrix for this Population on given Data
    /// @param data Data object to compute importances on
    /// @param n_perm Number of permutations (default: 1000)
    /// @param used_only Whether to compute importances only for features used in the population (default: TRUE)
    /// @param seed Random seed for reproducibility (default: 4815162342)
    /// @return Matrix (data.frame) with rows = features, columns = individuals
    /// @export
    pub fn compute_importance_matrix(&self, data: &Data, n_perm: Option<i32>, used_only: Option<bool>, seed: Option<u64>) -> Result<Robj> {
        let permutations = n_perm.unwrap_or(1000).max(1) as usize;
        let used_only = used_only.unwrap_or(true);

        let feature_set: HashSet<usize> = if used_only {
            let mut hs = HashSet::new();
            for ind in &self.intern.individuals {
                for &f in ind.features.keys() { hs.insert(f); }
            }
            hs
        } else {
            (0..data.intern.feature_len).collect()
        };
        let mut rows: Vec<usize> = feature_set.into_iter().collect();
        rows.sort_unstable();

        let n_feat = rows.len();
        let mut row_index: HashMap<usize, usize> = HashMap::with_capacity(n_feat);
        let mut rownames: Vec<String> = Vec::with_capacity(n_feat);
        for (rpos, fidx) in rows.iter().copied().enumerate() {
            row_index.insert(fidx, rpos);
            let fname = data.intern.features.get(fidx)
                .cloned().unwrap_or_else(|| format!("f{}", fidx));
            rownames.push(fname);
        }

        let base_seed = seed.unwrap_or(4815162342);
        let mut rng = ChaCha8Rng::seed_from_u64(base_seed);
        let mut feature_seeds_global: HashMap<usize, Vec<u64>> = HashMap::with_capacity(n_feat);
        for &f in &rows {
            let mut s = Vec::with_capacity(permutations);
            for _ in 0..permutations { s.push(rng.next_u64()); }
            feature_seeds_global.insert(f, s);
        }

        let n_ind = self.intern.individuals.len();
        let threads = resolve_threads()
            .or_else(|| std::thread::available_parallelism().ok().map(|n| n.get()))
            .unwrap_or(1);

        let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        let row_index = std::sync::Arc::new(row_index);
        let rows_arc = std::sync::Arc::new(rows);
        let seeds_arc = std::sync::Arc::new(feature_seeds_global);

        let columns: Vec<Vec<f64>> = pool.install(|| {
            self.intern.individuals.par_iter().map(|ind| {
                let feats_for_ind: Vec<usize> = if used_only {
                    let mut v: Vec<usize> = ind.features.keys().copied()
                        .filter(|f| row_index.contains_key(f))
                        .collect();
                    v.sort_unstable();
                    v
                } else {
                    rows_arc.as_ref().clone()
                };

                let mut col = vec![0.0_f64; n_feat];

                let ic = ind.compute_oob_feature_importance(
                    &data.intern,
                    permutations,
                    &feats_for_ind,
                    &seeds_arc
                );

                for imp in &ic.importances {
                    if let Some(&r) = row_index.get(&imp.feature_idx) {
                        col[r] = imp.importance;
                    }
                }
                col
            }).collect()
        });

        let mut values = vec![0.0_f64; n_feat * n_ind];
        for j in 0..n_ind {
            for r in 0..n_feat {
                values[r + n_feat * j] = columns[j][r];
            }
        }
        let colnames: Vec<String> = (0..n_ind).map(|j| format!("Individual {}", j)).collect();

        let mut mat = Doubles::from_values(values).into_robj();
        mat.set_attrib("dim", r!([n_feat as i32, n_ind as i32]))?;
        let rn = Strings::from_values(rownames);
        let cn = Strings::from_values(colnames);
        mat.set_attrib("dimnames", list!(rn, cn))?;
        Ok(mat)
    }

}

///////////////////////////////////////////////////////////////
/// Jury object
///////////////////////////////////////////////////////////////

#[extendr]
#[derive(Debug, Clone)]
pub struct Jury {
    intern: GJury,
    data: Arc<GData>,
    param: Arc<GParam>,
}

/// @title Jury
/// @name Jury
/// @description gpredomics Jury object
/// @details 
/// Jury represents an ensemble of expert models (a calibrated population) that make predictions
/// through voting and weighting schemes. It implements various voting methods (majority, consensus)
/// to aggregate predictions from multiple experts.
/// 
/// @section Methods:
/// \describe{
///   \item{\code{new_from_param(population, param)}}{Constructs a Jury object from a population and parameters. Parameters: population (Population of experts to form the Jury), param (Parameters for the Jury). Returns Jury object.}
///   \item{\code{from_population(population, threshold, window)}}{Create a Jury directly from a Population without applying any filters. Parameters: population (Population to use as experts), threshold (Voting threshold, default 0.5), window (Threshold window percentage, default 5.0). Returns Jury object.}
///   \item{\code{get()}}{Returns an R object containing all Jury fields for R interface. Returns list with Jury fields.}
///   \item{\code{get_metrics()}}{Get the base metrics already computed and stored in this Jury (AUC, accuracy, sensitivity, specificity, rejection_rate). No data required, no computation. Returns a list with base metrics only.}
///   \item{\code{compute_metrics(data)}}{Compute all metrics (base + additional: MCC, NPV, PPV, F1-score, G-mean) on the provided Data object. Returns a list with all computed metrics.}
///   \item{\code{predict(data)}}{Predict classes and scores on the provided Data object. Returns a list with two elements: class (predicted classes) and score (predicted scores).}
///   \item{\code{fit(data)}}{Fit/calibrate the Jury on new data by recomputing weights, thresholds, and all metrics. Parameter: data (Data object used for calibration).}
///   \item{\code{get_population()}}{Extract the population from the jury (experts). Returns Population object containing the experts.}
///   \item{\code{print_self_report()}}{Display the Jury report with self-evaluation metrics to console.}
///   \item{\code{print_report(test_data)}}{Display the Jury report with both training and test data evaluation. Parameter: test_data (Data object for test evaluation).}
///   \item{\code{address()}}{Get memory address of this Jury object. Returns string representing the memory address.}
///   \item{\code{print()}}{Get summary of this Jury. Returns string representing the Jury summary.}
/// }
/// @export
#[extendr]
impl Jury {

    /// @title Constructs a Jury object
    /// @name Jury$new_from_param
    /// @description Constructs a Jury object from a population and parameters
    /// @param population Population of experts to form the Jury
    /// @param param Parameters for the Jury
    /// @return Jury object
    /// @export
    pub fn new_from_param(population: &Population, param: &Param) -> Jury {
        // Create jury using the voting parameters from Param
        let intern_jury = GJury::new(
            &population.intern.clone(),
            &param.intern.voting.min_perf,
            &param.intern.voting.min_diversity,
            &param.intern.voting.method,
            &param.intern.voting.method_threshold,
            &param.intern.voting.threshold_windows_pct,
            &gpredomics::voting::WeightingMethod::Uniform  // Default weighting
        );
        Jury {
            intern: intern_jury,
            data: Arc::clone(&population.data),
            param: Arc::clone(&population.param),   
        }
    }

    /// @title Create Jury from Population
    /// @name Jury$from_population
    /// @description Create a Jury directly from a Population without applying any filters
    /// @param population Population to use as experts (no filtering applied)
    /// @param threshold Voting threshold (default: 0.5)
    /// @param window Threshold window percentage (default: 5.0)
    /// @return Jury object
    /// @export
    pub fn from_population(population: &Population, threshold: f64, window: f64) -> Jury {
        let mut internal_jury = GJury::new(
            &population.intern.clone(),
            &0.0,        // min_perf = 0.0 
            &0.0,        // min_diversity = 0.0 
            &gpredomics::voting::VotingMethod::Majority,
            &threshold,  // voting_threshold
            &window,     // threshold_window
            &gpredomics::voting::WeightingMethod::Uniform
        );

        internal_jury.evaluate(&population.data);
        
        Jury {
            intern: internal_jury,
            data: Arc::clone(&population.data),
            param: Arc::clone(&population.param),
        }
    }

    /// @title Get Jury fields
    /// @name Jury$get
    /// @description Returns an R object containing all Jury fields for R interface
    /// @return List with Jury fields
    /// @export
    pub fn get(&self) -> Robj {

        let weights: Vec<f64> = match &self.intern.weights {
            Some(weights) => weights.clone(),
            None => vec![]
        };

        let predicted_classes = match &self.intern.predicted_classes {
            Some(classes) => {
                let classes_i32: Vec<i32> = classes.iter().map(|&x| x as i32).collect();
                Robj::from(classes_i32)
            },
            None => Robj::from(Vec::<i32>::new())
        };

        let experts_individuals = {
            let data_arc = Arc::new(self.data.clone());
            let param_arc = Arc::new(self.param.clone());
            self.intern.experts.individuals
                .iter()
                .map(|gi| Individual::new(gi, 
                                        Arc::clone(&data_arc), 
                                        Arc::clone(&param_arc)).get())
                .collect::<Vec<Robj>>()
        };

        let jury_fields = vec![
            ("experts", Robj::from(experts_individuals)),
            ("voting_method", Robj::from(format!("{:?}", self.intern.voting_method))),
            ("voting_threshold", Robj::from(self.intern.voting_threshold)),
            ("threshold_window", Robj::from(self.intern.threshold_window)),
            ("weights", Robj::from(weights)),
            ("weighting_method", Robj::from(format!("{:?}", self.intern.weighting_method))),
            ("auc", Robj::from(self.intern.auc)),
            ("accuracy", Robj::from(self.intern.accuracy)),
            ("sensitivity", Robj::from(self.intern.sensitivity)),
            ("specificity", Robj::from(self.intern.specificity)),
            ("rejection_rate", Robj::from(self.intern.rejection_rate)),
            ("predicted_classes", predicted_classes),
            ("rust", self.clone().into_robj())
        ];

        let jury_robj = List::from_pairs(jury_fields);
    
        jury_robj.into_robj()
    }

    /// @title Get already computed base metrics
    /// @name Jury$get_metrics
    /// @description Get the already computed base metrics stored in this Jury (no computation).
    /// Returns only the core metrics that are always calculated during calibration: AUC, accuracy, 
    /// sensitivity, specificity, and rejection_rate. For additional metrics (MCC, NPV, PPV, F1-score, G-mean),
    /// use compute_metrics(data).
    /// @return A list with base metrics (auc, accuracy, sensitivity, specificity, rejection_rate)
    pub fn get_metrics(&self) -> Robj {
        list!(
            auc = self.intern.auc.into_robj(),
            accuracy = self.intern.accuracy.into_robj(),
            sensitivity = self.intern.sensitivity.into_robj(),
            specificity = self.intern.specificity.into_robj(),
            rejection_rate = self.intern.rejection_rate.into_robj(),
        ).into()
    }

    /// @title Compute all metrics on new dataset
    /// @name Jury$compute_metrics
    /// @description Compute all metrics including base metrics (AUC, accuracy, sensitivity, specificity, rejection rate)
    /// and additional metrics (MCC, NPV, PPV, F1-score, G-mean) on the provided Data object.
    /// @param data Data object used for metric computation
    /// @return A list with all computed metrics (base + additional)
    /// @export
    pub fn compute_metrics(&self, data: &Data) -> Robj {

        let additional_metrics = AdditionalMetrics { 
            mcc: Some(0.0),
            npv: Some(0.0),
            ppv: Some(0.0),
            f1_score: Some(0.0),
            g_mean: Some(0.0),
        };

        let mut jury = self.intern.clone();
        jury.metrics = additional_metrics;

        let (auc, accuracy, sensitivity, specificity, rejection_rate, additional) = jury.compute_new_metrics(&data.intern);
        
        list!(
            auc = auc.into_robj(),
            accuracy = accuracy.into_robj(),
            sensitivity = sensitivity.into_robj(),
            specificity = specificity.into_robj(),
            rejection_rate = rejection_rate.into_robj(),
            mcc = additional.mcc.into_robj(),
            npv = additional.npv.into_robj(),
            ppv = additional.ppv.into_robj(),
            f1_score = additional.f1_score.into_robj(),
            g_mean = additional.g_mean.into_robj(),
        ).into()

    }

    /// @title Predict classes and scores
    /// @name Jury$predict
    /// @description Predict classes and scores on the provided Data object.
    /// @param data Data object used for prediction
    /// @return A list with two elements: class (predicted classes) and score (predicted scores)
    /// @export
    pub fn predict(&self, data: &Data) -> Robj {
        let (classes, scores) = self.intern.predict(&data.intern);
        
        list!(class=(classes.iter().map(|x| {*x as i32})).collect::<Vec<i32>>().into_robj(), score=scores.into_robj()).into()
    }

    /// @title Fit/calibrate the Jury
    /// @name Jury$fit
    /// @description Fit/calibrate the Jury on new data by recomputing weights, thresholds, and all metrics.
    /// @param data Data object used for calibration
    /// @export
    pub fn fit(&mut self, data: &Data) {
        self.intern.evaluate(&data.intern);
    }

    /// @title Display Jury with training data
    /// @name Jury$display_train
    /// @description Display of the Jury with only training data
    /// @param data Data object used for display
    /// @param param Parameters for the Jury
    /// @return R object with training display
    /// @export
    pub fn print_self_report(&self) {
        let display_text = self.intern.display(&self.data, None, &self.param);
        rprintln!("{}", display_text);
    }

    /// @title Display Jury with training and test data
    /// @name Jury$display_train_and_test
    /// @description Display of the Jury with training and test data
    /// @param data Data object used for training display
    /// @param test_data Data object used for test display
    /// @param param Parameters for the Jury
    /// @return R object with display
    /// @export
    pub fn print_report(&mut self, test_data: &Data) {
        let display_text = self.intern.display(&self.data, Some(&test_data.intern), &self.param);
        rprintln!("{}", display_text);
    }

    /// @title Extract population from jury
    /// @name Jury$get_population
    /// @description Extract the population from the jury (experts)
    /// @return Population object containing the experts
    /// @export
    pub fn get_population(&self) -> Population {
        Population::new(
            self.intern.experts.clone(),
            Arc::clone(&self.data),
            Arc::clone(&self.param),
        )
    }
    
    /// @title Get Jury memory address
    /// @name Jury$address
    /// @description Get memory address of this Jury object
    /// @return String representing the memory address
    /// @export
    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern as *const _)
    }
    
    /// @title Get Jury summary
    /// @name Jury$print
    /// @description Get summary of this Jury
    /// @return String representing the Jury summary
    /// @export
    pub fn print(&self) -> String {
        let addr = self.address();
        
        let voting_method = match self.intern.voting_method {
            gpredomics::voting::VotingMethod::Majority => "Majority",
            gpredomics::voting::VotingMethod::Consensus => "Consensus",
        };
        
        let weighting_method = match self.intern.weighting_method {
            gpredomics::voting::WeightingMethod::Uniform => "Uniform",
            gpredomics::voting::WeightingMethod::Specialized { .. } => "Specialized",
        };
        
        let expert_count = self.intern.experts.individuals.len();
        let effective_experts = if let Some(ref weights) = self.intern.weights {
            weights.iter().filter(|&&w| w > 0.0).count()
        } else {
            expert_count
        };
        
        format!(
            "@{}\n\
            Jury: {}:{} (Experts: {})\n\
             Performance Metrics\n\
             \tAUC       : {:.3}\n\
             \tAccuracy  : {:.3}\n\
             \tSensitivity: {:.3}\n\
             \tSpecificity: {:.3}\n\
             \tRejection : {:.3}\n\
             Configuration\n\
             \tVoting    : {}\n\
             \tWeighting : {}\n\
             \tThreshold : {:.3}\n\
             \tWindow    : {:.1}%\n\
             Population\n\
             \tTotal Experts: {}\n\
             \tEffective    : {}",
            addr, voting_method, weighting_method, expert_count,
            self.intern.auc,
            self.intern.accuracy,
            self.intern.sensitivity,
            self.intern.specificity,
            self.intern.rejection_rate,
            voting_method,
            weighting_method,
            self.intern.voting_threshold,
            self.intern.threshold_window,
            expert_count,
            effective_experts
        )
    }
}

/// Custom log format function for flexi_logger.
/// 
/// Formats log entries with timestamp, level, and message.
/// Format: "YYYY-MM-DD HH:MM:SS [LEVEL] message"
fn custom_format(
    w: &mut dyn std::io::Write,
    now: &mut flexi_logger::DeferredNow,
    record: &log::Record,
) -> std::io::Result<()> {
    write!(
        w,
        "{} [{}] {}",
        now.now().format("%Y-%m-%d %H:%M:%S"), // Format timestamp
        record.level(),
        record.args()
    )
}

/// Resolve the number of threads to use from R options.
/// 
/// Checks the R option 'gpredomics.threads.number' and returns the value if valid.
/// Accepts integer, numeric, or string representations of positive integers.
/// 
/// # Returns
/// - `Some(usize)` if a valid thread count is found in R options
/// - `None` if no option is set or the value is invalid
fn resolve_threads() -> Option<usize> {
    if let Ok(robj) = R!("getOption('gpredomics.threads.number', NULL)") {
        if let Some(i) = robj.as_integer() { if i > 0 { return Some(i as usize); } }
        if let Some(f) = robj.as_real() { if f > 0.0 { return Some(f as usize); } }
        if let Some(s) = robj.as_str() { if let Ok(v) = s.parse::<usize>() { if v > 0 { return Some(v); } } }
    }
    None
}



///////////////////////////////////////////////////////////////
/// Glogger object
///////////////////////////////////////////////////////////////


// === Global logger handle & helpers (idempotent init) ===

/// Global logger handle using OnceCell for thread-safe initialization.
static LOGGER_HANDLE: OnceCell<LoggerHandle> = OnceCell::new();

/// Initialize the global logger or return existing handle.
/// 
/// This function ensures the logger is initialized only once, even if called
/// from multiple threads. Subsequent calls return a clone of the existing handle.
/// 
/// # Arguments
/// * `builder` - Closure that constructs the Logger
/// 
/// # Returns
/// A LoggerHandle that can be used to control the logger
fn init_or_get_with<F>(builder: F) -> LoggerHandle
where
    F: FnOnce() -> Logger,
{
    if let Some(h) = LOGGER_HANDLE.get() {
        return h.clone();
    }
    let handle = builder()
        .start()
        .expect("Failed to initialize logger");
    // If two threads race, this will be Err for one; both handles are clones of the same Arc internally.
    let _ = LOGGER_HANDLE.set(handle.clone());
    handle
}

/// Get the current logger handle if it has been initialized.
/// 
/// # Returns
/// - `Some(LoggerHandle)` if logger was previously initialized
/// - `None` if logger has not been initialized yet
fn current_handle() -> Option<LoggerHandle> {
    LOGGER_HANDLE.get().cloned()
}

#[extendr]
pub struct GLogger {
    handle: flexi_logger::LoggerHandle
}

/// @title GLogger
/// @name GLogger
/// @description An object to handle Logger
/// @details 
/// GLogger provides a configurable logging interface for the gpredomics package.
/// It supports different logging levels (info, debug, error, etc.) and can output
/// to screen or files with customizable formatting.
/// 
/// @section Methods:
/// \describe{
///   \item{\code{new()}}{Create a new screen logger with default 'info' level. Returns GLogger object.}
///   \item{\code{level(level)}}{Create a new screen logger with specified logging level. Parameter: level (logging level string, e.g., "info", "debug", "error"). Returns GLogger object.}
///   \item{\code{get(param)}}{Create a new logger from a Param object containing logging configuration. Parameter: param (Param object containing logging configuration). Returns GLogger object.}
///   \item{\code{set_level(level)}}{Change logging level. Parameter: level (new logging level string, e.g., "info", "debug", "error").}
/// }
/// @export
#[extendr]
impl GLogger {
    
    /// @title Create new screen logger
    /// @name GLogger$new
    /// @description Create a new screen logger
    /// @return GLogger object
    /// @export
    pub fn new() -> Self {
        let handle = init_or_get_with(|| {
            Logger::try_with_str("info")
                .expect("invalid log spec")
                .write_mode(WriteMode::BufferAndFlush)
        });
        Self { handle }
    }

    /// @title Create logger with level
    /// @name GLogger$level
    /// @description Create a new screen logger with specified logging level
    /// @param level Logging level string (e.g., "info", "debug", "error")
    /// @return GLogger object
    /// @export
    pub fn level(level: String) -> Self {
        if let Some(h) = current_handle() {
            if let Ok(spec) = LogSpecification::parse(&level) {
                h.set_new_spec(spec);
            }
            return Self { handle: h };
        }
        let handle = init_or_get_with(|| {
            Logger::try_with_str(&level)
                .expect("invalid log spec")
                .write_mode(WriteMode::BufferAndFlush)
        });
        Self { handle }
    }

    /// @title Create logger from Param
    /// @name GLogger$get
    /// @description Create a new logger from a Param
    /// @param param Param object containing logging configuration
    /// @return GLogger object
    /// @export  
    pub fn get(param: &Param) -> Self {
        let level = &param.intern.general.log_level;
        let log_base = &param.intern.general.log_base;
        let log_suffix = &param.intern.general.log_suffix;

        if let Some(h) = current_handle() {
            if let Ok(spec) = LogSpecification::parse(level) {
                h.set_new_spec(spec);
            }
            return Self { handle: h };
        }

        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();

        let handle = if !log_base.is_empty() {
            init_or_get_with(|| {
                Logger::try_with_str(level)
                    .expect("invalid log spec")
                    .log_to_file(
                        FileSpec::default()
                            .basename(log_base)
                            .suffix(log_suffix)
                            .discriminant(&timestamp),
                    )
                    .write_mode(WriteMode::BufferAndFlush)
                    .format_for_files(custom_format)
                    .format_for_stderr(custom_format)
            })
        } else {
            init_or_get_with(|| {
                Logger::try_with_str(level)
                    .expect("invalid log spec")
                    .write_mode(WriteMode::BufferAndFlush)
            })
        };

        Self { handle }
    }

    /// @title Change logging level
    /// @name GLogger$set_level
    /// @description Change logging level
    /// @param level New logging level string (e.g., "info", "debug", "error")
    /// @export
    pub fn set_level(&mut self, level: String) {
        if let Ok(new_spec) = LogSpecification::parse(level) {
            self.handle.set_new_spec(new_spec);
        }
    }

}

/// @title Run genetic algorithm
/// @name fit
/// @description Run the Gpredomics genetic algorithm to produce an Experiment from parameters
/// @details
/// This is the main entry point for running Gpredomics algorithms (GA, BEAM, MCMC).
/// The function loads data specified in the param object, runs the algorithm, and
/// returns a complete Experiment object with results. If param$voting$vote is TRUE,
/// a Jury ensemble is automatically computed.
/// 
/// @param param Param object containing experiment configuration (data paths, algorithm, CV settings)
/// @param running_flag RunningFlag object to control and monitor execution
/// @return Experiment object containing populations, final results, and optional Jury
/// @export
#[extendr]
pub fn fit(param: &Param, running_flag: &RunningFlag) -> Experiment {
    let mut p = param.intern.clone();
    let n_threads = resolve_threads();
    p.general.thread_number = if n_threads.is_some() {
        n_threads.unwrap()
    } else {
        p.general.thread_number
    };

    let intern = run(&p, running_flag.get_arc());
    
    let train_data_arc = Arc::new(intern.train_data.clone());
    let param_arc = Arc::new(param.intern.clone());
    let mut exp = Experiment {
        intern,
        train_data_arc,
        param_arc,
    };
    
    if param.intern.voting.vote {
        exp.intern.compute_voting();
    }
    exp
}

/// @title Run genetic algorithm with Data object
/// @name fit_on
/// @description Run the Gpredomics algorithm using a pre-loaded Data object
/// @details
/// This function runs Gpredomics algorithms (GA, BEAM, MCMC) on data already loaded in R.
/// Unlike fit(), this function uses an existing Data object instead of loading from file paths.
/// If param$voting$vote is TRUE, a Jury ensemble is automatically computed.
/// 
/// @param data Data object containing the training data (features and labels)
/// @param param Param object containing experiment configuration (algorithm, CV settings, etc.)
/// @param running_flag RunningFlag object to control and monitor execution
/// @return Experiment object containing populations, final results, and optional Jury
/// @export
#[extendr]
pub fn fit_on(data: &Data, param: &Param, running_flag: &RunningFlag) -> Experiment {
    let mut p = param.intern.clone();
    let n_threads = resolve_threads();
    p.general.thread_number = if n_threads.is_some() {
        n_threads.unwrap()
    } else {
        p.general.thread_number
    };

    p.data.holdout_ratio = 0.0; // No holdout when data is provided directly

    let intern = run_on_data(data.intern.clone(), None, &p, running_flag.get_arc());
    
    let train_data_arc = Arc::new(intern.train_data.clone());
    let param_arc = Arc::new(param.intern.clone());
    let mut exp = Experiment {
        intern,
        train_data_arc,
        param_arc,
    };
    
    if param.intern.voting.vote {
        exp.intern.compute_voting();
    }
    exp
}

// Macro to expose the struct and methods to R
extendr_module! {
    mod gpredomicsR;
    impl RunningFlag;
    impl Experiment;
    impl Param;
    impl GLogger;
    impl Individual;
    impl Data;
    impl Population;
    impl Jury;
    fn fit;
    fn fit_on;
    fn as_gpredomics_data;
}

/*#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
*/