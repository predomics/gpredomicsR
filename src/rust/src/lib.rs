#![allow(non_snake_case)]

// Note: extendr only export the documentation above impl blocks
// It's then necessary to repeat the struct documentation above impl blocks
// to have them appear in the R documentation.

///////////////////////////////////////////////////////////////
/// Main dependencies
///////////////////////////////////////////////////////////////

use extendr_api::prelude::*;
use gpredomics::param::Param as GParam;
use gpredomics::data::Data as GData;
use gpredomics::param::get as GParam_get;
use gpredomics::population::Population as GPopulation;
use gpredomics::individual::Individual as GIndividual;
use gpredomics::experiment::Experiment as GExperiment;
use gpredomics::experiment::ExperimentMetadata;
use gpredomics::experiment::{ImportanceAggregation, ImportanceCollection, ImportanceScope, ImportanceType};
use gpredomics::experiment::Jury  as GJury;
use gpredomics::{ run_ga, run_beam, run_mcmc };
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

/// Print message to R console and return as Error
#[inline]
fn r_print(msg: impl AsRef<str>) -> extendr_api::Error {
    let s = msg.as_ref();
    rprintln!("[gpredomicsR] {}", s);
    Error::Other(s.into())
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
    pub fn load(file_path: String) -> Self {
        if let Ok(this_param) = GParam_get(file_path.clone()) {
            Self {
                intern: this_param
            }
        }
        else {
            panic!("Unsuitable param file: {}",file_path)
        }
    } 


    /// @title Get Individual description
    /// @name Individual$get
    /// @description Retrieves a full description of an individual, including features and related statistics.
    /// @param self The Individual instance
    /// @return R object containing individual details
    pub fn get(&self) -> Robj {
        // Convert General fields
        let general = List::from_pairs(vec![
            ("seed", Robj::from(self.intern.general.seed)),
            ("algo", Robj::from(self.intern.general.algo.clone())),
            ("language", Robj::from(self.intern.general.language.clone())),
            ("data_type", Robj::from(self.intern.general.data_type.clone())),
            ("data_type_epsilon", Robj::from(self.intern.general.data_type_epsilon.clone())),
            ("thread_number", Robj::from(self.intern.general.thread_number)),
            ("log_base", Robj::from(self.intern.general.log_base.clone())),
            ("log_suffix", Robj::from(self.intern.general.log_suffix.clone())),
            ("log_level", Robj::from(self.intern.general.log_level.clone())),
            ("fit", Robj::from(format!("{:?}",self.intern.general.fit))),
            ("k_penalty", Robj::from(self.intern.general.k_penalty.clone())),
            ("fr_penalty", Robj::from(self.intern.general.fr_penalty.clone())),
            ("gpu", Robj::from(self.intern.general.gpu.clone())),
        ]);

        // Convert Data fields
        let data = List::from_pairs(vec![
            ("X", Robj::from(self.intern.data.X.clone())),
            ("y", Robj::from(self.intern.data.y.clone())),
            ("Xtest", Robj::from(self.intern.data.Xtest.clone())),
            ("ytest", Robj::from(self.intern.data.ytest.clone())),
            ("feature_maximal_number_per_class", Robj::from(self.intern.data.feature_maximal_number_per_class.clone())),
            ("feature_selection_method", Robj::from(self.intern.data.feature_selection_method.clone())),
            ("feature_minimal_prevalence_pct", Robj::from(self.intern.data.feature_minimal_prevalence_pct)),
            ("feature_maximal_pvalue", Robj::from(self.intern.data.feature_maximal_pvalue)),
            ("feature_minimal_feature_value", Robj::from(self.intern.data.feature_minimal_feature_value)),
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
            ("method", Robj::from(self.intern.beam.method.clone())),
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
            ("overfit_penalty", Robj::from(self.intern.cv.overfit_penalty.clone())),
            ("inner_folds", Robj::from(self.intern.cv.inner_folds)),
            ("outer_folds", Robj::from(self.intern.cv.outer_folds)),
            ("fit_on_valid", Robj::from(self.intern.cv.fit_on_valid)),
            ("cv_best_models_ci_alpha", Robj::from(self.intern.cv.cv_best_models_ci_alpha))
        ]);
        
        // Convert Voting fields
        let voting = List::from_pairs(vec![
            ("vote", Robj::from(self.intern.voting.vote)),
            ("use_fbm", Robj::from(self.intern.voting.use_fbm)),
            ("min_perf", Robj::from(self.intern.voting.min_perf)),
            ("min_diversity", Robj::from(self.intern.voting.min_diversity)),
            ("method", Robj::from(format!("{:?}",self.intern.voting.method))),
            ("method_threshold", Robj::from(self.intern.voting.method_threshold)),
            ("threshold_windows_pct", Robj::from(self.intern.voting.threshold_windows_pct)),
            ("complete_display", Robj::from(self.intern.voting.complete_display)),
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
            "epsilon" => self.intern.general.data_type_epsilon = value,
            "thread_number" => self.intern.general.thread_number = value as usize,
            "k_penalty" => self.intern.general.k_penalty = value,
            "fr_penalty" => self.intern.general.fr_penalty = value,
            "n_model_to_display" => self.intern.general.n_model_to_display = value as u32,
            "display_level" => self.intern.general.display_level = value as usize,
            
            // Data parameters
            "feature_minimal_prevalence_pct" => self.intern.data.feature_minimal_prevalence_pct = value,
            "feature_maximal_pvalue" => self.intern.data.feature_maximal_pvalue = value,
            "feature_minimal_feature_value" => self.intern.data.feature_minimal_feature_value = value,
            "feature_minimal_log_abs_bayes_factor" => self.intern.data.feature_minimal_log_abs_bayes_factor = value,
            "feature_maximal_number_per_class" => self.intern.data.feature_maximal_number_per_class = value as usize,
            
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
            "outer_folds" => self.intern.cv.outer_folds = value as usize,
            "cv_best_models_ci_alpha" => self.intern.cv.cv_best_models_ci_alpha = value,
            
            // Voting parameters
            "min_perf" => self.intern.voting.min_perf = value,
            "min_diversity" => self.intern.voting.min_diversity = value,
            "method_threshold" => self.intern.voting.method_threshold = value,
            "threshold_windows_pct" => self.intern.voting.threshold_windows_pct = value,
            
            // Importance parameters
            "n_permutations_oob" => self.intern.importance.n_permutations_oob = value as usize,
            
            // GPU parameters
            "max_total_memory_mb" => self.intern.gpu.max_total_memory_mb = value as u64,
            "max_buffer_size_mb" => self.intern.gpu.max_buffer_size_mb = value as u32,

            "X"|"y"|"Xtest"|"ytest"|"pvalue_method"|"algo"|"language"|"data_type"|"fit" => panic!("Use param$set_string() for {}",variable),
            "keep_all_generations"|"gpu" => panic!("Use param$set_bool() for {}",variable),
            "log_level"|"log_base"|"log_suffix" => panic!("Cannot set logs this way, create or get back your GLogger object"),
            _ => panic!("Unknown variable: {} ", variable)
        }
    }

    pub fn address(&self) -> String {
        format!("0x{:p}", &self.intern)
    }

    /// @description Set a string parameter by name.
    /// @param variable Name of the parameter to set
    /// @param string New string value for the parameter
    /// @export
    pub fn set_string(&mut self, variable: &str, string: String) {
        match variable {
            // Data parameters
            "X" => self.intern.data.X = string,
            "y" => self.intern.data.y = string,
            "Xtest" => self.intern.data.Xtest = string,
            "ytest" => self.intern.data.ytest = string,
            "feature_selection_method" => self.intern.data.feature_selection_method = string,
            
            // General parameters
            "algo" => self.intern.general.algo = string,
            "language" => self.intern.general.language = string,
            "data_type" => self.intern.general.data_type = string,
            
            // BEAM parameters
            "method" => self.intern.beam.method = string,
            
            // MCMC parameters
            "fit" => self.intern.general.fit = match string.to_lowercase().as_str() { 
                "auc" => gpredomics::param::FitFunction::auc,
                "specificity" => gpredomics::param::FitFunction::specificity,
                "sensitivity" => gpredomics::param::FitFunction::sensitivity,
                "mcc" =>  gpredomics::param::FitFunction::mcc,
                _ => panic!("Unknown fit function: {}", string)
             }, //self.intern.general.fit = string,
            "log_level"|"log_base"|"log_suffix" => panic!("Cannot set logs this way, create or get back your GLogger object"),
            _ => panic!("Variable unknown or not settable byt set_string: {} ", variable)
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
            "inverse_classes" => self.intern.data.inverse_classes = value,
            
            // CV parameters
            "fit_on_valid" => self.intern.cv.fit_on_valid = value,
            
            // Voting parameters
            "vote" => self.intern.voting.vote = value,
            "use_fbm" => self.intern.voting.use_fbm = value,
            "complete_display" => self.intern.voting.complete_display = value,
            "specialized" => self.intern.voting.specialized = value,
            
            // Importance parameters
            "compute_importance" => self.intern.importance.compute_importance = value,
            "scaled_importance" => self.intern.importance.scaled_importance = value,
            
            // GPU parameters
            "fallback_to_cpu" => self.intern.gpu.fallback_to_cpu = value,
            
            "log_level"|"log_base"|"log_suffix" => panic!("Cannot set logs this way, create or get back your GLogger object"),
            _ => panic!("Unknown boolean variable: {}", variable),
            }
    }

}

///////////////////////////////////////////////////////////////
/// Data object
///////////////////////////////////////////////////////////////

/// Convert a sparse matrix into an R data frame.
/// 
/// - `matrix`: HashMap<(row_index, col_index), value> 
/// - `column_names`: names for each column (length = number of columns)
/// - `row_names`: names for each row (length = number of rows)
/// 
/// Any missing (row, col) in `matrix` is treated as 0.0 in the data frame.
///
/// Example usage from R side (after wrapping with `.Call`):
/// ```r
/// df <- .Call("my_module_sparse_matrix_to_dataframe", sparse_mat, col_names, row_names)
/// print(df)
/// ```

/// Convert a sparse matrix to a DataFrame
/// @param matrix Sparse matrix represented as a HashMap with keys (col_index, row_index) and values as f64
/// @param column_names Names for each column (length = number of columns)
/// @param row_names Names for each row (length = number of rows)
/// @return DataFrame as an R object
/// @export
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

/// Convert a sparse vector to a dense vector
/// @param vector Sparse vector represented as a HashMap with keys as indices and values as u8
/// @param length Length of the resulting dense vector
/// @return Dense vector as an R object
/// @export
fn sparse_vector_to_vector(vector: &HashMap<usize, u8>, length: usize) -> Robj {
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
/// @return Gpredomics Data object with aligned and processed X and y.
/// @export
#[extendr]
pub fn as_gpredomics_data(df: Robj, y_vec: Robj, features_in_columns: bool) -> extendr_api::Result<Data> {
    // DataFrame as list
    let df_list = df.as_list().ok_or_else(|| r_print("df must be a data.frame (list-like)."))?;
    let ncol = df_list.len();

    // First column to infer nrow
    let first_col = df_list.values().next().ok_or_else(|| r_print("data.frame has zero columns."))?;
    let nrow = if let Some(v) = first_col.as_real_vector() {
        v.len()
    } else if let Some(v) = first_col.as_integer_vector() {
        v.len()
    } else if let Some(v) = first_col.as_logical_vector() {
        v.len()
    } else if let Some(v) = first_col.as_string_vector() {
        v.len()
    } else {
        return Err(r_print("Unsupported first column type; expected numeric/integer/logical/character."));
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
        r_print("y has no names(); set names(y) to chosen sample names before calling.");
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
                return Err(r_print(&format!("y factor must have at least two valid levels; found {}.", valid_codes.len())));
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
                r_print(&format!("Warning: {} samples with values outside the two main factor levels will be classified as 'unknown' (class 2).", unknown_count));
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
                        r_print(&format!("Warning: {} samples with values outside {{0,1}} will be classified as 'unknown' (class 2).", unknown_count));
                    }
                    (mapped, vec!["0".into(), "1".into(), "unknown".into()])
                },
                [a, b] => {
                    let (lo, hi) = (*a, *b);
                    r_print(&format!("y has two distinct integer values {{{},{}}}; mapping {}→0, {}→1, others→2.", lo, hi, lo, hi));
                    let mapped: Vec<u8> = vals.iter().map(|&z| {
                        if z == lo { 0 }
                        else if z == hi { 1 }
                        else { 2 } // Unknown
                    }).collect();
                    let unknown_count = mapped.iter().filter(|&&x| x == 2).count();
                    if unknown_count > 0 {
                        r_print(&format!("Warning: {} samples with values outside {{{},{}}} will be classified as 'unknown' (class 2).", unknown_count, lo, hi));
                    }
                    (mapped, vec![lo.to_string(), hi.to_string(), "unknown".into()])
                },
                _ => {
                    // More than 2 distinct values: use first two, rest→unknown
                    r_print(&format!("y has {} distinct integer values; using first two ({}, {}), others→unknown.", uniq.len(), uniq[0], uniq[1]));
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
            return Err(r_print(&format!("y must have at least two string classes; found {}.", uniq.len())));
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
            r_print(&format!("Warning: {} samples with values outside {{{}, {}}} will be classified as 'unknown' (class 2).", unknown_count, class0, class1));
        }
        
        (mapped, vec![class0, class1, "unknown".into()])
    } else {
        return Err(r_print("Unsupported y type; provide integer, character, or 2-level factor."));
    };

    if y_names_raw.len() != y_labels_u8.len() {
        return Err(r_print("Length mismatch between names(y) and y values."));
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
                .ok_or_else(|| r_print(&format!(
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
                return Err(r_print(&format!(
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
                return Err(r_print(&format!(
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
                return Err(r_print(&format!(
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
            return Err(r_print(&format!(
                "Character column at {} not supported; encode to numeric first.",
                col_idx
            )));
        } else {
            return Err(r_print(&format!(
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
        let data = List::from_pairs(vec![
            ("X", sparse_matrix_to_dataframe(&self.intern.X, &self.intern.samples, &self.intern.features)),
            ("y", Robj::from(self.intern.y.iter().map(|x|{*x as i32}).collect::<Vec<i32>>())),
            ("features", Robj::from(&self.intern.features)),
            ("samples", Robj::from(&self.intern.samples)),
            ("feature_class", sparse_vector_to_vector(&self.intern.feature_class, self.intern.feature_len)),
            // indexes start by 1 in R and need to be adapted
            ("feature_selection", Robj::from(&self.intern.feature_selection
                .iter()
                .map(|x| {x+1})
                .collect::<Vec<usize>>()
            )),
            ("feature_len", Robj::from(&self.intern.feature_len)),
            ("sample_len", Robj::from(&self.intern.sample_len)),
            ("classes", Robj::from(&self.intern.classes))
            
        ]);

        // Return the data as an R list object
        Robj::from(data)
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
             \tSparse Elements: {}",
            addr, total_samples, total_features,
            total_samples,
            total_features,
            selected_features, selection_status,
            sparsity_pct,
            class0_count,
            class1_count,
            other_count,
            sparse_elements
        )
    }

}

///////////////////////////////////////////////////////////////
/// Individual object
///////////////////////////////////////////////////////////////

#[derive(Debug, Clone)]
#[extendr]
pub struct Individual {
    intern: GIndividual,
    data: Arc<GData>
}


impl Individual {
    
    /// Create a new Individual object filled with a GIndividual
    pub fn new(ind: &GIndividual, data: Arc<GData>) -> Self {
        Self {
            intern: ind.clone(),
            data,
        }
    }

}

/// @title Individual
/// @name Individual
/// @description gpredomicsR proxy object for Individual
/// @details 
/// Individual represents a single model from Gpredomics.
/// It contains features, coefficients, and various performance metrics (AUC, accuracy, sensitivity, specificity).
/// 
/// @section Methods:
/// \describe{
///   \item{\code{get()}}{Retrieves a full description of an individual, including features and related statistics. Returns an R object containing individual details.}
///   \item{\code{compute_auc(data)}}{Compute AUC for this individual on the provided Data object. Returns a new Individual with updated AUC.}
///   \item{\code{compute_metrics(data)}}{Compute threshold/accuracy/sensitivity/specificity for this individual on the provided Data object. Returns a new Individual with updated metrics.}
///   \item{\code{compute_all(data)}}{Compute auc/threshold/accuracy/sensitivity/specificity for this individual on the provided Data object. Returns a new Individual with updated AUC and metrics.}
///   \item{\code{evaluate()}}{Compute individual score. Returns R object containing the individual score.}
///   \item{\code{predict_class_and_score(data)}}{Use individual on a data object to provide predicted class and score. Returns a list with two elements: class (predicted class) and score (predicted score).}
///   \item{\code{predict(data)}}{Return a list of predicted class for the samples in the data. Returns a vector of predicted classes (0 or 1).}
///   \item{\code{to_string()}}{Print the individual. Returns string representation of the Individual.}
///   \item{\code{address()}}{Get memory address of this Individual object. Returns string representing the memory address.}
///   \item{\code{print()}}{Print as Gpredomics style with detailed formatting.}
///   \item{\code{set_threshold(threshold)}}{Set the threshold of the individual.}
///   \item{\code{compute_importance(data, n_perm, seed, used_only)}}{Compute feature importance for this individual on the provided data. Parameters: n_perm (number of permutations, default 1000), seed (optional seed for random number generation), used_only (whether to compute importance only for features used in the individual, default true). Returns a named numeric vector of feature importances.}
///   \item{\code{get_genealogy(experiment, max_depth, include_metrics)}}{Retrieve the genealogy (ancestry tree) of this individual across generations. Parameters: experiment (Experiment object containing all generations), max_depth (maximum depth to traverse, default 10), include_metrics (whether to include AUC/k/generation, default TRUE). Returns a list with nodes and edges data.frames ready for igraph/ggraph visualization. Use with plot_genealogy() helper.}
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

    /// @title Compute AUC
    /// @name Individual$compute_auc
    /// @description Compute AUC for this individual on the provided Data object.
    /// @param data The Data object to compute AUC on
    /// @return A new Individual with updated AUC
    pub fn compute_auc(&self, data: &Data) -> Individual {
        let mut gi = self.clone();
        gi.intern.compute_auc(&data.intern);
        gi
    }    

    /// @title Compute metrics
    /// @name Individual$compute_metrics
    /// @description Compute threshold/accuracy/sensitivity/specificity for this individual on the provided Data object.
    /// @param data The Data object to compute metrics on
    /// @return A new Individual with updated metrics
    pub fn compute_metrics(&self, data: &Data) -> Individual {
        let mut gi  = self.clone();
        let i = &mut gi.intern;
        (i.accuracy, i.sensitivity, i.specificity) = i.compute_metrics(&data.intern);
        gi
    }

    /// @title Compute all metrics
    /// @name Individual$compute_all
    /// @description Compute auc, threshold, accuracy, sensitivity and specificity for this individual on the provided Data object.
    /// @param data The Data object to compute all metrics on
    /// @return A new Individual with updated AUC and metrics
    pub fn compute_all(&self, data: &Data) -> Individual {
        let mut gi  = self.clone();
        gi.intern.compute_auc(&data.intern);
        let i = &mut gi.intern;
        (i.accuracy, i.sensitivity, i.specificity) = i.compute_metrics(&data.intern);
        gi
    }

    /// @title Evaluate individual
    /// @name Individual$evaluate
    /// @description Compute individual score (evaluate on its associated data).
    /// @return R object containing the individual score
    pub fn evaluate(&self) -> Robj {
        self.intern.evaluate(&self.data).into_robj()
    }

    /// @title Predict class and score
    /// @name Individual$predict_class_and_score
    /// @description Use individual on a Data object to provide predicted class and score.
    /// @param data The Data object to predict on
    /// @return A list with two elements: class (predicted class) and score (predicted score)
    pub fn predict_class_and_score(&self, data: &Data) -> Robj {
        let (classes, scores) = self.intern.evaluate_class_and_score(&data.intern);
        list!(class=(classes.iter().map(|x| {*x as i32})).collect::<Vec<i32>>().into_robj(), score=scores.into_robj()).into()
    }

    /// @title Predict classes
    /// @name Individual$predict
    /// @description Return predicted classes for samples in the provided Data object.
    /// @param data The Data object to predict on
    /// @return A vector of predicted classes (0 or 1)
    pub fn predict(&self, data: &Data) -> Robj {
        self.intern.evaluate(&data.intern).into_iter().map(|x| {if x>=self.intern.threshold {1} else {0}}).collect::<Vec<i32>>()
        .into_robj()
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
        let s = &self.intern.display(&self.data, None, &"ga".to_string(), 2, false);
        for line in s.lines() {
            rprintln!("{}", line);
        }
    }

    /// @title Set threshold
    /// @name Individual$set_threshold
    /// @description Set the threshold of the individual used for binary predictions.
    /// @param threshold The new threshold value
    pub fn set_threshold(&mut self, threshold: f64) {
        self.intern.threshold = threshold;
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

    /// @title Get individual genealogy
    /// @name Individual$get_genealogy
    /// @description Retrieve the genealogy (ancestry tree) of this individual across generations.
    /// Returns a list with two data.frames: `nodes` (containing node information) and `edges` (parent-child relationships),
    /// plus metadata. The output is designed to work with igraph/ggraph for visualization.
    /// @param experiment The Experiment object containing all generations
    /// @param max_depth Maximum depth to traverse in the genealogy tree (default: 10)
    /// @param include_metrics Whether to include performance metrics (AUC, k, generation) for each node (default: TRUE)
    /// @return A list with components: `nodes` (data.frame with columns: id, label, auc, k, generation),
    /// `edges` (data.frame with columns: from, to, depth), and `metadata` (list with depth info).
    /// Use with `plot_genealogy()` for visualization.
    /// @export
    pub fn get_genealogy(&self, experiment: &Experiment, max_depth: Option<i32>, include_metrics: Option<bool>) -> Result<Robj>  {
        let depth = max_depth.unwrap_or(10).max(1) as usize;
        let with_metrics = include_metrics.unwrap_or(true);

        let mut pops_vec: Vec<GPopulation> = Vec::new();
        let mut meta_by_hash: HashMap<u64, (f64, i32, i32)> = HashMap::new();
        for (gen, pop) in experiment.intern.collections[0].iter().enumerate() {
            if with_metrics {
                for ind in &pop.individuals {
                    meta_by_hash.entry(ind.hash).or_insert((ind.auc, ind.k as i32, gen as i32));
                }
            }
            pops_vec.push(pop.clone());
        }
        

        let genealogy_map = self.intern.get_genealogy(&pops_vec, depth);

        let mut node_id: Vec<String> = Vec::new();
        let mut node_label: Vec<String> = Vec::new();
        let mut node_auc: Vec<f64> = Vec::new();
        let mut node_k: Vec<i32> = Vec::new();
        let mut node_gen: Vec<i32> = Vec::new();

        let mut edge_from: Vec<String> = Vec::new();
        let mut edge_to: Vec<String> = Vec::new();
        let mut edge_depth: Vec<i32> = Vec::new();

        for ((hash, parents_opt), depth_set) in genealogy_map {
            let id = format!("H{:016x}", hash);
            node_id.push(id.clone());
            node_label.push(format!("H#{}", &id[1..9]));

            let (auc, k, gen) = meta_by_hash.get(&hash).cloned().unwrap_or((0.0, 0, -1));
            node_auc.push(auc);
            node_k.push(k);
            node_gen.push(gen);

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

        let nodes_df = data_frame!(
            id = node_id,
            label = node_label,
            auc = node_auc,
            k = node_k,
            generation = node_gen
        );
        let edges_df = data_frame!(
            from = edge_from,
            to = edge_to,
            depth = edge_depth
        );
        let metadata = list!(
            n_nodes = nodes_df.nrows() as i32,
            n_edges = edges_df.nrows() as i32,
            max_depth = depth as i32,
            root_id = format!("H{:016x}", self.intern.hash)
        );

        Ok(list!(nodes = nodes_df, edges = edges_df, metadata = metadata).into())
    }
}

///////////////////////////////////////////////////////////////
/// Experiment object
///////////////////////////////////////////////////////////////

#[extendr]
pub struct Experiment {
    intern: GExperiment,
    train_data_arc: Arc<GData>,
}

/// @title Experiment
/// @name Experiment
/// @description A global Experiment object that proxies all the different Rust gpredomics objects under the hood
/// @details 
/// Experiment encapsulates the entire genetic programming workflow including training data, test data,
/// parameters, and algorithm results (populations across generations, final population, etc.).
/// 
/// @section Methods:
/// \describe{
///   \item{\code{individual(generation, order)}}{Retrieves a full description of an individual from a specified generation and order. Parameters: generation (i32 generation index), order (i32 order within generation). Returns Individual object.}
///   \item{\code{test_data()}}{Retrieves the test data associated with the experiment. Returns Data object representing the test data.}
///   \item{\code{train_data()}}{Retrieves the training data associated with the experiment. Returns Data object representing the training data.}
///   \item{\code{get_data_robj(train)}}{Retrieves the data associated with the experiment as an R object. Parameter: train (logical; if TRUE, returns training data, otherwise test data). Returns R object representing the data.}
///   \item{\code{get_data(train)}}{Retrieves the data associated with the experiment as a Data object. Parameter: train (logical; if TRUE, returns training data, otherwise test data). Returns Data object.}
///   \item{\code{get_generation(generation)}}{Retrieves descriptions of all individuals from a specified generation. Parameter: generation (i32 generation index). Returns R list object encapsulating features and metrics of all individuals in the generation.}
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
    pub fn individual(&self, generation: i32, order:i32) -> Individual {
        if self.intern.collections.len() > 0 {
            Individual::new(
                &self.intern.collections[0][generation as usize].individuals[order as usize], 
                Arc::clone(&self.train_data_arc)
            )
        } else {
            panic!("Cannot extract an individual from an experiment without collection")
        }
    }

    /// @title Get test data
    /// @name Experiment$test_data
    /// @description Retrieves the test data associated with the experiment.
    /// @return Data object representing the test data.
    pub fn test_data(&self) -> Data {
        if let Some(test_data) = &self.intern.test_data {
            Data {
                intern: test_data.clone()
            }
        } else {
            panic!("No test data attached to this experiment")
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
        else {self.test_data().get()}
    }

    /// @title Get data as Data object
    /// @name Experiment$get_data
    /// @description Retrieves the data associated with the experiment as a Data object.
    /// @param train Logical; if TRUE, returns training data, otherwise test data
    /// @return Data object
    pub fn get_data(&self, train:bool) -> Data {
        if train {self.train_data()}
        else {self.test_data()}
    }
    
    /// @title Get generation
    /// @name Experiment$get_generation
    /// @description Retrieves descriptions of all individuals from a specified generation.
    /// @param generation An i32 specifying the generation index.
    /// @return An R list object (Robj) encapsulating features and metrics of all individuals in the generation.
    pub fn get_generation(&self, generation: i32) -> Robj {
        if self.intern.collections.len() > 0 {
            self.intern.collections[0][generation as usize].individuals.iter()
            .map(|i| Individual::new(i, Arc::clone(&self.train_data_arc)).get())
            .collect::<Vec<Robj>>()
            .into_robj()
        }
        else {
            panic!("No collection found in this experiment")
        }
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
    /// @description Get the number of generations in the experiment/population
    /// @return Number of generations
    pub fn generation_number(&self) -> i32 {
        if self.intern.collections.len() > 0 {
            self.intern.collections[0].len() as i32
        } else {
            0
        }
    }

    /// @title Get population size for a generation
    /// @param generation An i32 specifying the generation index.
    /// @return Size (number of individuals) of the requested generation
    pub fn population_size(&self, generation: i32) -> i32 {
        if self.intern.collections.len() > 0 {
            self.intern.collections[0][generation as usize].individuals.len() as i32
        }
        else {
            panic!("Cannot size any generation in this experiment because it has none")
        }
    }

    /// @title Load external dataset
    /// @name Experiment$load_data
    /// @description Load an external dataset to evaluate the model
    /// @param x_path Path to the X data file
    /// @param y_path Path to the y data file
    /// @return Data object containing the loaded data
    pub fn load_data(&self, x_path: String, y_path: String) -> Data {
        let mut gdata = GData::new();
        let _ = gdata.load_data(&x_path, &y_path);
        if !self.intern.train_data.check_compatibility(&gdata) {
            panic!("Data not compatible with training data");
        }
        gdata.set_classes(self.intern.train_data.classes.clone());
        Data {intern:gdata}
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
    pub fn get_jury(&self) -> Jury {
        match &self.intern.others {
            Some(ExperimentMetadata::Jury { jury }) => {
                Jury {
                    intern: jury.clone(),
                    data: Arc::clone(&self.train_data_arc),
                }
            },
            _ => panic!("No Jury computed during this experiment"),
        }
    }

    /// @title Load serialized experiment
    /// @name Experiment$load
    /// @description Load a serialized experiment from disk
    /// @param path Path to the experiment file
    /// @return Loaded Experiment object
    pub fn load(path: String) -> Self {
        if let Ok(experiment) = GExperiment::load_auto(&path) {
            let train_data_arc = Arc::new(experiment.train_data.clone());
            Experiment {
                intern: experiment,
                train_data_arc,
            }
        }
        else {
            panic!("Could not read an experiment from this file: {}",&path)
        }
    }

    /// @title Save experiment
    /// @name Experiment$save
    /// @description Save an experiment to disk
    /// @param path Path to save the experiment
    pub fn save(&self, path: String) {
        match self.intern.save_auto(&path) {
            Ok(_) => {
                println!("Experiment saved in {}", path);
            }
            Err(e) => {
                // On error, propagate as a panic so R receives an error instead of a silent stderr message.
                panic!("Could not save an experiment to file {}. Error: {}", path, e);
            }
        }
    }

    /// @title Extract population from experiment
    /// @description Extract population from experiment, optionally specifying generation number
    /// @param generation Optional generation number (0-based). If None, returns final population
    /// @return Population object for the specified generation or final population
    pub fn get_population(&self, generation: Option<i32>) -> Population {
        match generation {
            None => {
                if let Some(ref final_pop) = self.intern.final_population {
                    Population {
                        intern: final_pop.clone(),
                        data: Arc::clone(&self.train_data_arc),
                    }
                } else {
                    panic!("No final population available in this experiment. Run the algorithm first.")
                }
            },
            Some(gen_idx) => {
                let gen_idx = gen_idx as usize;

                if !self.intern.collections.is_empty() {
                    if let Some(first_collection) = self.intern.collections.get(0) {
                        if let Some(population) = first_collection.get(gen_idx) {
                            return Population {
                                intern: population.clone(),
                                data: Arc::clone(&self.train_data_arc),
                            };
                        }
                    }
                    panic!("Generation {} not found in CV collections", gen_idx);
                }
                
                panic!("Generation history not available. Only final population is stored in non-CV mode.")
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

}

///////////////////////////////////////////////////////////////
/// Population object
///////////////////////////////////////////////////////////////

#[extendr]
#[derive(Debug, Clone)]
pub struct Population {
    intern: GPopulation,
    data: Arc<GData>,
}

impl Population {
    /// Internal constructor - creates Population with shared data reference
    pub(crate) fn new(intern: GPopulation, data: Arc<GData>) -> Self {
        Self { intern, data }
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
///   \item{\code{predict_scores_matrix(data)}}{Predict all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts). Parameter: data (Data object to predict on). Returns dataframe with predicted scores.}
///   \item{\code{predict_classes_matrix(data)}}{Predict classes for all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts, Values = predicted classes 0 or 1). Parameter: data (Data object to predict on). Returns dataframe with predicted classes.}
///   \item{\code{filter_by_auc(min_auc)}}{Filter population by AUC threshold. Parameter: min_auc (minimum AUC threshold). Returns filtered Population object.}
///   \item{\code{filter_by_fit(min_fit)}}{Filter population by fitness threshold. Parameter: min_fit (minimum fit threshold). Returns filtered Population object.}
///   \item{\code{filter_by_diversity(min_diversity_pct, by_niche)}}{Filter population by diversity using Jaccard dissimilarity. Parameters: min_diversity_pct (minimum diversity percentage 0-100), by_niche (whether to compute diversity within niches). Returns filtered Population object.}
///   \item{\code{filter_by_sensitivity(min_sensitivity)}}{Filter population by sensitivity threshold. Parameter: min_sensitivity (minimum sensitivity threshold). Returns filtered Population object.}
///   \item{\code{filter_by_specificity(min_specificity)}}{Filter population by specificity threshold. Parameter: min_specificity (minimum specificity threshold). Returns filtered Population object.}
///   \item{\code{filter_by_mask(mask)}}{Filter population using a logical vector (1/0). Parameter: mask (integer vector 1/0 indicating which individuals to keep). Returns filtered Population object.}
///   \item{\code{filter_by_k(min_k, max_k)}}{Filter population by number of features (k). Parameters: min_k (minimum number of features), max_k (maximum number of features). Returns filtered Population object.}
///   \item{\code{get_fbm(alpha, min_pct_fallback)}}{Get Family of Best Models (FBM) using confidence interval selection. This method selects models with performance statistically equivalent to the best model. Parameters: alpha (confidence level, default 0.05 for 95% confidence), min_pct_fallback (if FBM selection fails, minimum percentage to keep, default 5.0). Returns Population object containing the FBM.}
///   \item{\code{fit_on(data, fit_function, k_penalty, thread_number)}}{Recompute fitness metrics for all individuals on new data. Parameters: data (new Data object to fit on), fit_function (fitness function to use: "auc", "mcc", "sensitivity", "specificity"), k_penalty (penalty coefficient for model complexity, default 0.0), thread_number (number of threads to use, default 4).}
///   \item{\code{address()}}{Get memory address of this Population object. Returns string representing the memory address.}
///   \item{\code{get_individual(index)}}{Get an Individual of a population by index. Parameter: index (index of the individual to retrieve). Returns Individual object at the specified index.}
///   \item{\code{print()}}{Get comprehensive print information about the population. Returns string representing the Population summary.}
///   \item{\code{from_individuals(individuals)}}{Create a Population from a vector or list of R Individual objects. Parameter: individuals (R vector or list of Individual objects, must have at least one). Returns Population object created from the individuals.}
///   \item{\code{extend(other)}}{Extend this population with another population (in-place modification). Parameter: other (another Population object to add).}
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
            .map(|gi: &GIndividual| Individual::new(gi, Arc::clone(&self.data)).get())
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
    /// @name Population$predict_scores_matrix
    /// @description Predict all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts)
    /// @param data The Data object to predict on
    /// @return Dataframe with predicted scores
    /// @export
    pub fn predict_scores_matrix(&self, data: &Data) -> Robj {
        let num_samples = data.intern.sample_len;
        let num_individuals = self.intern.individuals.len();
        
        if num_individuals == 0 {
            panic!("Population is empty - cannot predict");
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
        df_list.set_attrib(row_names_symbol(), Robj::from(&self.data.samples)).ok();
        df_list.set_class(&["data.frame"]).ok();
        df_list.into()
    }
    
    /// @title Predict classes matrix
    /// @name Population$predict_classes_matrix
    /// @description Predict classes for all individuals of the population on data and return a dataframe (Rows = samples, Columns = individuals/experts, Values = predicted classes 0 or 1)
    /// @param data The Data object to predict on
    /// @return Dataframe with predicted classes
    /// @export  
    pub fn predict_classes_matrix(&self, data: &Data) -> Robj {
        let num_samples = data.intern.sample_len;
        let num_individuals = self.intern.individuals.len();
        
        if num_individuals == 0 {
            panic!("Population is empty - cannot predict");
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
        df_list.into()
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
        
        Population::new(gpopulation, Arc::clone(&self.data))
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
        
        Population::new(gpopulation, Arc::clone(&self.data))
    }

    /// @title Filter population by diversity using Jaccard dissimilarity
    /// @param min_diversity_pct Minimum diversity percentage (0-100)
    /// @param by_niche Whether to compute diversity within niches
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_diversity(&self, min_diversity_pct: f64, by_niche: bool) -> Population {
        let filtered_pop = self.intern.filter_by_signed_jaccard_dissimilarity(min_diversity_pct, by_niche);
        
        Population::new(filtered_pop, Arc::clone(&self.data))
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
        
        Population::new(gpopulation, Arc::clone(&self.data))
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
        
        Population::new(gpopulation, Arc::clone(&self.data))
    }

    /// @title Filter population by mask
    /// @name Population$filter_by_mask
    /// @description Filter population using a logical vector (1/0)
    /// @param mask Integer vector (1/0) indicating which individuals to keep
    /// @return Filtered Population object
    /// @export
    pub fn filter_by_mask(&self, mask: Vec<i32>) -> Population {
        if mask.len() != self.intern.individuals.len() {
            panic!("Mask length ({}) must match population size ({})", 
                   mask.len(), self.intern.individuals.len());
        }
        
        let mut gpopulation = GPopulation::new();
        
        for (individual, keep) in self.intern.individuals.iter().zip(mask.iter()) {
            if *keep != 0 { 
                gpopulation.individuals.push(individual.clone());
            }
        }
        
        Population::new(gpopulation, Arc::clone(&self.data))
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
        
        Population::new(gpopulation, Arc::clone(&self.data))
    }

    /// @title Get Family of Best Models (FBM)
    /// @name Population$get_fbm
    /// @description Get Family of Best Models (FBM) using confidence interval selection. This method selects models with performance statistically equivalent to the best model
    /// @param alpha Confidence level (default: 0.05 for 95% confidence)
    /// @param min_pct_fallback If FBM selection fails, minimum percentage to keep (default: 5.0)
    /// @return Population object containing the FBM
    /// @export
    pub fn get_fbm(&self, alpha: f64, min_pct_fallback: f64) -> Population {
        let mut sorted_population = self.clone();
        sorted_population.intern = sorted_population.intern.sort();
        
        let fbm_population = sorted_population.intern.select_best_population(alpha);
        
        if fbm_population.individuals.is_empty() {
            warn!("FBM selection failed or returned empty population, using {}% fallback", min_pct_fallback);
            let (fallback_pop, _) = sorted_population.intern.select_first_pct(min_pct_fallback);
            Population::new(fallback_pop, Arc::clone(&self.data))
        } else {
            Population::new(fbm_population, Arc::clone(&self.data))
        }
    }

    /// @title Recompute fitness metrics
    /// @name Population$fit_on
    /// @description Recompute fitness metrics for all individuals on new data
    /// @param data New Data object to fit on
    /// @param fit_function Fitness function to use ("auc", "mcc", "sensitivity", "specificity")
    /// @param k_penalty Penalty coefficient for model complexity (default: 0.0)
    /// @param thread_number Number of threads to use (default: 4)
    /// @export
    pub fn fit_on(&mut self, data: &Data, fit_function: String, k_penalty: f64) {
        let threads = resolve_threads()
        .or_else(|| std::thread::available_parallelism().ok().map(|n| n.get()))
        .unwrap_or(1);
        
        match fit_function.to_lowercase().as_str() {
            "auc" => {
                self.intern.auc_fit(&data.intern, k_penalty, threads, true);
            },
            "mcc" => {
                // MCC fit requires additional penalty for sensitivity/specificity
                self.intern.mcc_fit(&data.intern, k_penalty, 1.0, threads);
            },
            "sensitivity" | "specificity" => {
                // Use objective fit with appropriate penalties
                let (fpr_penalty, fnr_penalty) = match fit_function.as_str() {
                    "sensitivity" => (1.0, 1.0), // Focus on sensitivity
                    "specificity" => (1.0, 1.0), // Focus on specificity  
                    _ => (1.0, 1.0),
                };
                self.intern.objective_fit(&data.intern, fpr_penalty, fnr_penalty, k_penalty, threads, true);
            },
            _ => {
                panic!("Unknown fit function: {}. Use 'auc', 'mcc', 'sensitivity', or 'specificity'", fit_function);
            }
        }
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
    pub fn get_individual(&self, index:i32) -> Individual {
        if self.intern.individuals.len() > 0 {
            Individual::new(&self.intern.individuals[index as usize], Arc::clone(&self.data))
        } else {
            panic!("Cannot extract an individual from an empty Population")
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
                0 => "Binary",  // BINARYLANG
                1 => "Ternary", // TERNARYLANG  
                2 => "Pow2",    // POW2LANG
                3 => "Ratio",   // RATIOLANG
                101 => "MCMC",  // MCMCGENERICLANG
                _ => "Unknown",
            };
            *language_counts.entry(lang).or_insert(0) += 1;
            
            let dtype = match individual.data_type {
                0 => "Raw",        // RAWTYPE
                1 => "Prevalence", // PREVALENCETYPE
                2 => "Log",        // LOGTYPE
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
    pub fn from_individuals(individuals: Robj) -> Population {
        let mut gpopulation = GPopulation::new();
        let mut data_arc: Option<Arc<GData>> = None;
        
        if let Some(r_list) = individuals.as_list() {
            for individual_robj in r_list.values() {
                match <&Individual>::try_from(&individual_robj) {
                    Ok(r_individual) => {
                        if data_arc.is_none() {
                            data_arc = Some(Arc::clone(&r_individual.data));
                        }
                        gpopulation.individuals.push(r_individual.intern.clone());
                    },
                    Err(_) => {
                        panic!("Invalid Individual object in list");
                    }
                }
            }
        }
        else if individuals.is_vector() {
            let len = individuals.len();
            for i in 0..len {
                match individuals.index(i) {
                    Ok(individual_robj) => {
                        match <&Individual>::try_from(&individual_robj) {
                            Ok(r_individual) => {
                                if data_arc.is_none() {
                                    data_arc = Some(Arc::clone(&r_individual.data));
                                }
                                gpopulation.individuals.push(r_individual.intern.clone());
                            },
                            Err(_) => {
                                panic!("Invalid Individual object at index {}", i);
                            }
                        }
                    },
                    Err(e) => {
                        panic!("Cannot access element at index {}: {:?}", i, e);
                    }
                }
            }
        }
        else {
            panic!("Input must be a vector or list of Individual objects");
        }
        
        if let Some(data) = data_arc {
            Population::new(gpopulation, data)
        } else {
            panic!("Cannot create Population from empty Individual list");
        }
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
        
        self.extend(&individuals_pop);
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
            Some(x) => return Err(r_print(format!("Unknown aggregation '{}'", x))),
        };
        let scaled = scaled.unwrap_or(true);
        let compact = compact.unwrap_or(false);

        let seed = seed.unwrap_or(4815162342);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);

        let mut ic: ImportanceCollection = self.intern.compute_pop_oob_feature_importance(
            &data.intern,
            permutations,
            &mut rng,
            &agg,
            scaled,
            /* cascade */ false,
            Some(1),
        ); // importance_type = MDA; scope = Population { id: pid } 

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

        let imp_label = if scaled { "scaled_importance" } else { "importance" };
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
///   \item{\code{evaluate(data)}}{Calibrates the expert population on the training data. Parameter: data (Data object used for calibration).}
///   \item{\code{predict_class_and_score(data)}}{Compute class and scores on a new dataset. Parameter: data (Data object used for prediction). Returns a list with two elements: class (integer vector) and score (numeric vector).}
///   \item{\code{compute_new_metrics(data)}}{Compute AUC/accuracy/sensitivity/rejection rate on a new dataset. Parameter: data (Data object used for metric computation). Returns a list with computed metrics.}
///   \item{\code{display_train(data, param)}}{Display of the Jury with only training data. Parameters: data (Data object used for display), param (Parameters for the Jury). Returns dataframe with training display.}
///   \item{\code{display_train_and_test(data, test_data, param)}}{Display of the Jury with training and test data. Parameters: data (Data object used for training display), test_data (Data object used for test display), param (Parameters for the Jury). Returns dataframe with display.}
///   \item{\code{get()}}{Returns an R object containing all Jury fields for R interface. Returns list with Jury fields.}
///   \item{\code{get_population()}}{Extract the population from the jury (experts). Returns Population object containing the experts.}
///   \item{\code{address()}}{Get memory address of this Jury object. Returns string representing the memory address.}
///   \item{\code{print()}}{Get summary of this Jury. Returns string representing the Jury summary.}
///   \item{\code{from_population(population)}}{Constructs a Jury object from a population using default parameters. Parameter: population (Population of experts). Returns Jury object.}
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
        let intern_jury = GJury::new_from_param(&population.intern, &param.intern);
        Jury {
            intern: intern_jury,
            data: Arc::clone(&population.data),
        }
    }
    
    /// @title Calibrates the expert population
    /// @name Jury$evaluate
    /// @description Calibrates the expert population on the training data
    /// @param data Data object used for calibration
    /// @export
    pub fn evaluate(&mut self, data: &Data) {
        self.intern.evaluate(&data.intern);
    }

    /// @title Compute class and scores
    /// @name Jury$predict_class_and_score
    /// @description Compute class and scores on a new dataset
    /// @param data Data object used for prediction
    /// @return A list with two elements: class (integer vector) and score (numeric vector)
    /// @export
    pub fn predict_class_and_score(&self, data: &Data) -> Robj {
        let (classes, scores) = self.intern.predict(&data.intern);
        
        list!(class=(classes.iter().map(|x| {*x as i32})).collect::<Vec<i32>>().into_robj(), score=scores.into_robj()).into()
    }

    /// @title Compute metrics on new dataset
    /// @name Jury$compute_new_metrics
    /// @description Compute AUC/accuracy/sensitivity/rejection rate on a new dataset 
    /// @param data Data object used for metric computation
    /// @return A list with computed metrics
    /// @export
    pub fn compute_new_metrics(&self, data: &Data) -> Robj {
        let (auc, accuracy, sensitivity, specificity, rejection_rate) = self.intern.compute_new_metrics(&data.intern);
        
        list!(
            auc = auc.into_robj(),
            accuracy = accuracy.into_robj(),
            sensitivity = sensitivity.into_robj(),
            specificity = specificity.into_robj(),
            rejection_rate = rejection_rate.into_robj()
        ).into()
    }

    /// @title Display Jury with training data
    /// @name Jury$display_train
    /// @description Display of the Jury with only training data
    /// @param data Data object used for display
    /// @param param Parameters for the Jury
    /// @return Dataframe with training display
    /// @export
    pub fn display_train(&mut self, data: &Data, param: &Param) -> Robj {
        let display_text = self.intern.display(&data.intern, None, &param.intern);
        display_text.into_robj()
    }

    /// @title Display Jury with training and test data
    /// @name Jury$display_train_and_test
    /// @description Display of the Jury with training and test data
    /// @param data Data object used for training display
    /// @param test_data Data object used for test display
    /// @param param Parameters for the Jury
    /// @export
    pub fn display_train_and_test(&mut self, data: &Data, test_data: &Data, param: &Param) -> Robj {
        let display_text = self.intern.display(&data.intern, Some(&test_data.intern), &param.intern);
        display_text.into_robj()
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
            self.intern.experts.individuals
                .iter()
                .map(|gi| Individual::new(gi, Arc::clone(&data_arc)).get())
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

    /// @title Extract population from jury
    /// @name Jury$get_population
    /// @description Extract the population from the jury (experts)
    /// @return Population object containing the experts
    /// @export
    pub fn get_population(&self) -> Population {
        Population::new(
            self.intern.experts.clone(),
            Arc::clone(&self.data)
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
            gpredomics::experiment::VotingMethod::Majority => "Majority",
            gpredomics::experiment::VotingMethod::Consensus => "Consensus",
        };
        
        let weighting_method = match self.intern.weighting_method {
            gpredomics::experiment::WeightingMethod::Uniform => "Uniform",
            gpredomics::experiment::WeightingMethod::Specialized { .. } => "Specialized",
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

    /// @title Create Jury from Population
    /// @name Jury$from_population
    /// @description Create a Jury directly from a Population without applying any filters
    /// @param population Population to use as experts (no filtering applied)
    /// @param threshold Voting threshold (default: 0.5)
    /// @param window Threshold window percentage (default: 5.0)
    /// @return Jury object
    /// @export
    pub fn from_population(population: &Population, threshold: f64, window: f64) -> Jury {
        let internal_jury = GJury::new(
            &population.intern.clone(),
            &0.0,        // min_perf = 0.0 
            &0.0,        // min_diversity = 0.0 
            &gpredomics::experiment::VotingMethod::Majority,
            &threshold,  // voting_threshold
            &window,     // threshold_window
            &gpredomics::experiment::WeightingMethod::Uniform
        );
        
        Jury {
            intern: internal_jury,
            data: Arc::clone(&population.data),
        }
    }
}

/// custom format for logs
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
static LOGGER_HANDLE: OnceCell<LoggerHandle> = OnceCell::new();

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
/// @description The simple genetic algorithm (ga) produce a Population from a Param object. The RunningFlag object is convenient when launching ga in a subthread, it must be provided (but you can let it live its own way)
/// @param param Param object containing experiment configuration
/// @param running_flag RunningFlag object to monitor execution
/// @return Experiment object containing notably the resulting Population
/// @export
#[extendr]
pub fn fit(param: &Param, running_flag: &RunningFlag) -> Experiment {
    let algo = &param.intern.general.algo;

    let mut p = param.intern.clone();
    let n_threads = resolve_threads();
    p.general.thread_number = if n_threads.is_some() {
        n_threads.unwrap()
    } else {
        p.general.thread_number
    };

    let intern = match algo.as_str() {
        "ga" => { run_ga(&p, running_flag.get_arc()) },
        "beam" => { run_beam(&p, running_flag.get_arc()) },
        "mcmc" => { run_mcmc(&p, running_flag.get_arc()) },
        _ => panic!("No such algo {}",algo)
    };
    
    let train_data_arc = Arc::new(intern.train_data.clone());
    let mut exp = Experiment {
        intern,
        train_data_arc,
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