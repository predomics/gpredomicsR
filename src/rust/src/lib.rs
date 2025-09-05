#![allow(non_snake_case)]

///////////////////////////////////////////////////////////////
/// Main dependencies
///////////////////////////////////////////////////////////////

use extendr_api::prelude::*;
//use gpredomics::param::Param;
use gpredomics::param::Param as GParam;
use gpredomics::data::Data as GData;
use gpredomics::param::get as GParam_get;
use gpredomics::population::Population  as GPopulation;
use gpredomics::individual::Individual  as GIndividual;
use gpredomics::experiment::Experiment  as GExperiment;
use gpredomics::experiment::ExperimentMetadata;
use gpredomics::experiment::Jury  as GJury;
use gpredomics::{ run_ga, run_beam, run_mcmc };

use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use flexi_logger::{Logger, LoggerHandle, WriteMode, FileSpec, LogSpecification};
use chrono::Local;
use std::collections::HashMap;
use once_cell::sync::OnceCell;

///////////////////////////////////////////////////////////////
/// Running Flag
///////////////////////////////////////////////////////////////
/// A struct to manage the `running` flag.
/// @export 
#[derive(Debug, Clone)]
#[extendr]
pub struct RunningFlag {
    flag: Arc<AtomicBool>,
}

/// @export 
#[extendr]
impl RunningFlag {
    /// Create a new `RunningFlag` (initially set to `true`).
    /// @export 
    pub fn new() -> Self {
        Self {
            flag: Arc::new(AtomicBool::new(true)),
        }
    }

    /// Set the `running` flag to `false`.
    /// @export 
    pub fn stop(&self) {
        self.flag.store(false, Ordering::Relaxed);
    }

    /// Check the current value of the `running` flag.
    /// @export 
    pub fn is_running(&self) -> bool {
        self.flag.load(Ordering::Relaxed)
    }

    /// Reset the `running` flag to `true`.
    /// @export 
    pub fn reset(&self) {
        self.flag.store(true, Ordering::Relaxed);
    }

}

/*#[extendr]
pub fn get_param() -> param::Param {
    param::Param::new()
}*/

impl RunningFlag {
    /// Get a clone of the inner Arc (useful for passing to long-running functions).
    fn get_arc(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.flag)
    }
}


///////////////////////////////////////////////////////////////
/// Param object
///////////////////////////////////////////////////////////////

/// gpredomics param object that create all settings
/// @export 
#[extendr]
#[derive(Clone)]
pub struct Param {
    intern: GParam
}

/// @export 
#[extendr]
impl Param {
    /// Create a new empty Param object
    /// @export 
    pub fn new() -> Self {
        Self {
            intern: GParam::new()
        }
    }

    /// Load a param.yaml to file to create a new Param
    /// @export 
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


    /// Returns an R object representing the current state of Param.
    /// @export
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
            ("mutation_non_null_chance_pct", Robj::from(self.intern.ga.mutation_non_null_chance_pct))
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
            ("cv", Robj::from(cv)),
            ("voting", Robj::from(voting)),
        ]);

        Robj::from(param_list)
    }   

    /// @export 
    pub fn set(&mut self, variable: &str, value: f64) {
        match variable {
            "feature_minimal_prevalence_pct" => self.intern.data.feature_minimal_prevalence_pct = value,
            "seed" => self.intern.general.seed = value as u64,
            "data_type_epsilon" => self.intern.general.data_type_epsilon = value,
            "thread_number" => self.intern.general.thread_number = value as usize,
            "k_penalty" => self.intern.general.k_penalty = value,
            "select_elite_pct" => self.intern.ga.select_elite_pct = value,
            "select_niche_pct" => self.intern.ga.select_niche_pct = value,
            "select_random_pct" => self.intern.ga.select_random_pct = value,
            // TODO @raynald continue here
            "X"|"y"|"Xtest"|"ytest"|"pvalue_method"|"algo"|"language"|"data_type"|"fit" => panic!("Use param$set_string() for {}",variable),
            "keep_all_generations"|"gpu" => panic!("Use param$set_bool() for {}",variable),
            "log_level"|"log_base"|"log_suffix" => panic!("Cannot set logs this way, create or get back your GLogger object"),
            _ => panic!("Unknown variable: {} ", variable)
        }
    }


    pub fn set_string(&mut self, variable: &str, string: String) {
        match variable {
            "X" => self.intern.data.X = string,
            "y" => self.intern.data.y = string,
            "Xtest" => self.intern.data.Xtest = string,
            "ytest" => self.intern.data.ytest = string,
            "feature_selection_method" => self.intern.data.feature_selection_method = string,
            "algo" => self.intern.general.algo = string,
            "language" => self.intern.general.language = string,
            "data_type" => self.intern.general.data_type = string,
            "fit" => self.intern.general.fit = match string.to_lowercase().as_str() { 
                "auc" => gpredomics::param::FitFunction::auc,
                "specificity" => gpredomics::param::FitFunction::specificity,
                "sensitivity" => gpredomics::param::FitFunction::sensitivity,
                _ => panic!("Unknown fit function: {}", string)
             }, //self.intern.general.fit = string,
            "log_level"|"log_base"|"log_suffix" => panic!("Cannot set logs this way, create or get back your GLogger object"),
            _ => panic!("Variable unknown or not settable byt set_string: {} ", variable)
        }
    }

    pub fn set_bool(&mut self, variable: &str, value: bool) {
        match variable {
            "gpu" => self.intern.general.gpu = value,
            "keep_trace" => self.intern.general.keep_trace = value,
            _ => panic!("Variable unknown or not settable byt set_bool: {} ", variable)
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

/// gpredomics Data object 
/// @export 
#[extendr]
#[derive(Clone)]
pub struct Data {
    intern: GData
}

/// @export 
#[extendr]
impl Data {
    /// Create a new empty Data object
    /// @export 
    pub fn new() -> Self {
        Self {
            intern: GData::new()
        }
    }

    
    /// Converts the Data struct fields to an R object.
    /// @export
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
}


///////////////////////////////////////////////////////////////
/// Individual object
///////////////////////////////////////////////////////////////

/// gpredomicsR proxy object for Individual
/// @export
#[derive(Debug, Clone)]
#[extendr]
pub struct Individual {
    intern: GIndividual,
    features: Vec<String>
}


impl Individual {
    
/// Create a new Individual object filled with a GIndividual
    pub fn new(ind:&GIndividual, data:&GData) -> Self {
        Self {
            intern: ind.clone(),
            features: ind.features
            .iter()
            .filter_map(|i| data.features.get(*i.0).cloned()) // sécurise contre out-of-bound
            .collect()
        }
    }

}


/// @export
#[extendr]
impl Individual {

    /// Create a new Individual object
    /// @export
    //pub fn new() -> Self {
    //    Self {
    //        intern: GIndividual::new(),
    //        features: Vec::new()
    //    }
    //}

    /// Retrieves a full description of an individual.
    /// This function returns an R object that includes individual features and related statistics.
    /// # Returns
    /// An R object (Robj) encapsulating the individual's features and metrics.
    ///
    /// # Example
    /// ```
    /// let robj_individual = model.get_individual_full(1, 2, true);
    /// ```
    ///
    /// @export
    pub fn get(&self) -> Robj {
        let mut coeff = Vec::new();
        let mut indexes = Vec::new();
        for (index, coefficient) in self.intern.features.iter() {
            coeff.push(*coefficient as i32);     // R integers are 32-bit
            indexes.push(*index as i32 + 1);     // indexes start by 1 in R and need to be adapted
        }
    
        // Build all fields including optionals before constructing the List
        let mut individual_fields = vec![
            ("features", Robj::from(self.features.clone())),
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



    /// Compute auc for this individual
    /// @export
    pub fn compute_auc(&mut self, data: &Data) {
        self.intern.compute_auc(&data.intern);
    }    

    /// Compute threshold/accuracy/sensitivity/specificity for this individual
    /// @export
    pub fn compute_metrics(&mut self, data: &Data) {
        let i = &mut self.intern;

        (i.accuracy, i.sensitivity, i.specificity) = 
                    i.compute_metrics(&data.intern);
    }

    /// Compute algorithm score
    /// @export
    pub fn evaluate(&self, data: &Data) -> Robj {
        self.intern.evaluate(&data.intern).into_robj()
    }

    /// Use individual on a data object to provide predicted class and score
    /// @export
    pub fn evaluate_class_and_score(&self, data: &Data) -> Robj {
        let (classes, scores) = self.intern.evaluate_class_and_score(&data.intern);
        
        list!(class=(classes.iter().map(|x| {*x as i32})).collect::<Vec<i32>>().into_robj(), score=scores.into_robj()).into()
    }

    /// Return a list of predicted class for the samples in the data
    /// @export
    pub fn predict(&self, data: &Data) -> Robj {
        self.intern.evaluate(&data.intern).into_iter().map(|x| {if x>=self.intern.threshold {1} else {0}}).collect::<Vec<i32>>()
        .into_robj()
    }

    /// Print the individual
    /// @export
    pub fn to_string(&self) -> String {
        format!("{:?}",&self.intern)
    }
    
}


///////////////////////////////////////////////////////////////
/// Experiment object
///////////////////////////////////////////////////////////////

/// A global Experiment object that proxies all the different Rust gpredomics objects under
/// the hood
/// @export
#[extendr]
pub struct Experiment {
    intern: GExperiment,
    // param: Param,
    // train_data: Data,
    // test_data: Data,
    // generations: Vec<GPopulation>
}

/// TODO add a load function to load a new Data (and check_compatibility)

/// @export
#[extendr]
impl Experiment {
    /// Retrieves a full description of an individual from a specified generation and order.
    /// This function returns an R object that includes individual features and related statistics.
    /// 
    /// # Arguments
    /// * `generation` - An i32 specifying the generation index.
    /// * `order` - An i32 specifying the order index within the generation.
    /// * `verbose` - A boolean flag that when true, prints the individual's details to the console.
    /// 
    /// # Returns
    /// An R object (Robj) encapsulating the individual's features and metrics.
    ///
    /// # Example
    /// ```
    /// let robj_individual = model.get_individual_full(1, 2, true);
    /// ```
    ///
    /// @export
    pub fn individual(&self, generation: i32, order:i32) -> Individual {
        if self.intern.collections.len() > 0 {
                Individual {
                    intern: self.intern.collections[0][generation as usize].individuals[order as usize].clone(),
                    features: self.intern.train_data.features.clone()
            }
        } else {
            panic!("Cannot extract an individual from an experiment without collection")
        }
    }

    /// @export
    pub fn test_data(&self) -> Data {
        if let Some(test_data) = &self.intern.test_data {
            Data {
                intern: test_data.clone()
            }
        } else {
            panic!("No test data attached to this experiment")
        }
    }

    /// @export
    pub fn train_data(&self) -> Data {
        Data {
            intern: self.intern.train_data.clone()
        }
    }

    /// @export
    pub fn get_data_robj(&self, train:bool) -> Robj {
        if train {self.train_data().get()}
        else {self.test_data().get()}
    }

    /// @export
    pub fn get_data(&self, train:bool) -> Data {
        if train {self.train_data()}
        else {self.test_data()}
    }
    
    /// Retrieves descriptions of all individuals from a specified generation.
    /// This function returns an R list of objects, each including individual features and related statistics.
    /// 
    /// # Arguments
    /// * `generation` - An i32 specifying the generation index.
    /// * `verbose` - A boolean flag that when true, prints the details of all individuals to the console.
    /// 
    /// # Returns
    /// An R list object (Robj) encapsulating the features and metrics of all individuals in the generation.
    ///
    /// # Example
    /// ```
    /// let robj_generation = model.get_generation(1, false);
    /// ```
    ///
    /// @export
    pub fn get_generation(&self, generation: i32) -> Robj {
        if self.intern.collections.len() > 0 {
            self.intern.collections[0][generation as usize].individuals.iter().cloned()
            .map(|i| {(Individual::new(&i, &self.intern.train_data)).get()})
            .collect::<Vec<Robj>>()
            .into_robj()
        }
        else {
            panic!("No collection found in this experiment")
        }
    } 

/*
    /// list all individuals at generation #generation
    /// @export
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

    /// get the number of generation included in the Population object
    /// @export
    pub fn generation_number(&self) -> i32 {
        if self.intern.collections.len() > 0 {
            self.intern.collections[0].len() as i32
        } else {
            0
        }
    }

    /// get the size (number of individuals) of a certain generation in a Population
    /// @export
    pub fn population_size(&self, generation: i32) -> i32 {
        if self.intern.collections.len() > 0 {
            self.intern.collections[0][generation as usize].individuals.len() as i32
        }
        else {
            panic!("Cannot size any generation in this experiment because it has none")
        }
    }

    /// load an external dataset to evaluate the model
    /// @export
    pub fn load_data(&self, x_path: String, y_path: String) -> Data {
        let mut gdata = GData::new();
        let _ = gdata.load_data(&x_path, &y_path);
        if !self.intern.train_data.check_compatibility(&gdata) {
            panic!("Data not compatible with training data");
        }
        gdata.set_classes(self.intern.train_data.classes.clone());
        Data {intern:gdata}
    }

    /// lget a param object    
    /// @export
    pub fn get_param(&self) -> Param {
        Param {
            intern: self.intern.parameters.clone()
        }
    }
    
    /// get a jury object
    /// @export
    pub fn get_jury(&self) -> Jury {
        match &self.intern.others {
            Some(ExperimentMetadata::Jury { jury }) => {
                Jury {
                    intern: jury.clone()
                }
            },
            _ => panic!("No Jury")
        }
    }

    /// load a serialized experiment
    /// @export
    pub fn load(path: String) -> Self {
        if let Ok(experiment) = GExperiment::load_auto(&path) {
            Experiment {
                intern: experiment
            }
        }
        else {
            panic!("Could not read an experiment from this file: {}",&path)
        }
    }

    /// save an experiment
    /// @export
    pub fn save(&self, path: String) {
        match self.intern.save_auto(&path) {
            Ok(_) => {
                println!("Experiment saved in {}", path);
            }
            Err(e) => {
                eprintln!("Could not save an experiment from this file: {}. Error: {}", path, e);
                // Optionally: panic! if you still want to crash
                // panic!("Save failed: {}", e);
            }
        }
    }

}

///////////////////////////////////////////////////////////////
/// Population object
///////////////////////////////////////////////////////////////

/// @export
#[extendr]
pub struct Population {
    intern: GPopulation,
}

/// @export
#[extendr]
impl Population {
    /// Create a new empty Population object
    /// @export 
    pub fn new() -> Self {
        Self {
            intern: GPopulation::new()
        }
    }

    /// Get population
    /// @export
    pub fn get(&self, data: &Data) -> Robj {
        let individuals = self.intern.individuals
            .iter()
            .cloned()
            .map(|gi: GIndividual| Individual::new(&gi, &data.intern).get())
            .collect::<Vec<Robj>>();
        List::from_pairs(vec![
            ("individuals", Robj::from(individuals)),
        ]).into_robj()
    }

}

///////////////////////////////////////////////////////////////
/// Jury object
///////////////////////////////////////////////////////////////

/// @export
#[extendr]
pub struct Jury {
    intern: GJury,
}

/// @export
#[extendr]
impl Jury {

    /// Constructs a Jury object
    /// @export
    pub fn new_from_param(population: &Population, param: &Param) -> Jury {
        let intern_jury = GJury::new_from_param(&population.intern, &param.intern);
        Jury {
            intern: intern_jury
        }
    }
    
    /// Calibrates the expert population on the training data
    /// @export
    pub fn evaluate(&mut self, data: &Data) {
        self.intern.evaluate(&data.intern);
    }

    /// Compute class and scores on a new dataset
    /// @export
    pub fn evaluate_class_and_score(&self, data: &Data) -> Robj {
        let (classes, scores) = self.intern.predict(&data.intern);
        
        list!(class=(classes.iter().map(|x| {*x as i32})).collect::<Vec<i32>>().into_robj(), score=scores.into_robj()).into()
    }

    /// Compute AUC/accuracy/sensitivity/rejection rate on a new dataset 
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

    /// Display of the Jury with only training data
    /// @export
    pub fn display_train(&mut self, data: &Data, param: &Param) -> Robj {
        let display_text = self.intern.display(&data.intern, None, &param.intern);
        display_text.into_robj()
    }

    /// Display of the Jury with training and test data
    /// @export
    pub fn display_train_and_test(&mut self, data: &Data, test_data: &Data, param: &Param) -> Robj {
        let display_text = self.intern.display(&data.intern, Some(&test_data.intern), &param.intern);
        display_text.into_robj()
    }

    /// Returns an R object containing all Jury fields for R interface
    /// @export
    pub fn get(&self, data: &Data) -> Robj {

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

        let experts_individuals = self.intern.experts.individuals
            .iter()
            .cloned()
            .map(|gi| Individual::new(&gi, &data.intern).get())
            .collect::<Vec<Robj>>();

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
        ];

        let jury_robj = List::from_pairs(jury_fields);
    
        jury_robj.into_robj()
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

/// An object to handle Logger
/// @export
#[extendr]
pub struct GLogger {
    handle: flexi_logger::LoggerHandle
}

/// @export
#[extendr]
impl GLogger {
    
    /// Create a new screen logger
    /// @export
    pub fn new() -> Self {
        let handle = init_or_get_with(|| {
            Logger::try_with_str("info")
                .expect("invalid log spec")
                .write_mode(WriteMode::BufferAndFlush)
        });
        Self { handle }
    }

    /// Create a new screen logger
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

    /// Create a new logger from a Param
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

    /// Change logging level
    /// @export
    pub fn set_level(&mut self, level: String) {
        if let Ok(new_spec) = LogSpecification::parse(level) {
            self.handle.set_new_spec(new_spec);
        }
    }

}

/// The simple genetic algorithm (ga) produce a Population from a Param object
/// the RunningFlag object is convenient when launching ga in a subthread, it must be
/// provided (but you can let it live its own way)
/// @export
#[extendr]
pub fn fit(param: &Param, running_flag: &RunningFlag) -> Experiment {
    let algo = &param.intern.general.algo;
    let mut exp = Experiment {
        intern: match algo.as_str() {
        "ga" => { run_ga(&param.intern, running_flag.get_arc()) },
        "beam" => { run_beam(&param.intern, running_flag.get_arc()) },
        "mcmc" => { run_mcmc(&param.intern, running_flag.get_arc()) },
        _ => panic!("No such algo {}",algo) }
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