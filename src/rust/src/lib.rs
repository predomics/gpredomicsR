///////////////////////////////////////////////////////////////
/// Main dependencies
///////////////////////////////////////////////////////////////

use extendr_api::prelude::*;
//use gpredomics::param::Param;
use gpredomics::param::Param as GParam;
use gpredomics::data::Data as GData;
use gpredomics::param::get as GParam_get;

use gpredomics::population::Population  as GPopulation;
use gpredomics::{ga_run, ga_no_overfit};

use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use flexi_logger::{Logger, WriteMode, FileSpec, LogSpecification, LevelFilter};
use chrono::Local;
use std::collections::HashMap;

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
            ("thread_number", Robj::from(self.intern.general.thread_number)),
            ("log_base", Robj::from(self.intern.general.log_base.clone())),
            ("log_suffix", Robj::from(self.intern.general.log_suffix.clone())),
            ("log_level", Robj::from(self.intern.general.log_level.clone())),
        ]);

        // Convert Data fields
        let data = List::from_pairs(vec![
            ("X", Robj::from(self.intern.data.X.clone())),
            ("y", Robj::from(self.intern.data.y.clone())),
            ("Xtest", Robj::from(self.intern.data.Xtest.clone())),
            ("ytest", Robj::from(self.intern.data.ytest.clone())),
            ("pvalue_method", Robj::from(self.intern.data.pvalue_method.clone())),
            ("feature_minimal_prevalence_pct", Robj::from(self.intern.data.feature_minimal_prevalence_pct)),
            ("feature_maximal_pvalue", Robj::from(self.intern.data.feature_maximal_pvalue)),
            ("feature_minimal_feature_value", Robj::from(self.intern.data.feature_minimal_feature_value)),
        ]);

        // Convert GA fields
        let ga = List::from_pairs(vec![
            ("population_size", Robj::from(self.intern.ga.population_size)),
            ("max_epochs", Robj::from(self.intern.ga.max_epochs)),
            ("min_epochs", Robj::from(self.intern.ga.min_epochs)),
            ("max_age_best_model", Robj::from(self.intern.ga.max_age_best_model)),
            ("kmin", Robj::from(self.intern.ga.kmin)),
            ("kmax", Robj::from(self.intern.ga.kmax)),
            ("kpenalty", Robj::from(self.intern.ga.kpenalty)),
            ("select_elite_pct", Robj::from(self.intern.ga.select_elite_pct)),
            ("select_random_pct", Robj::from(self.intern.ga.select_random_pct)),
            ("mutated_children_pct", Robj::from(self.intern.ga.mutated_children_pct)),
            ("mutated_features_pct", Robj::from(self.intern.ga.mutated_features_pct)),
            ("mutation_non_null_chance_pct", Robj::from(self.intern.ga.mutation_non_null_chance_pct)),
            ("feature_importance_permutations", Robj::from(self.intern.ga.feature_importance_permutations)),
            ("keep_all_generations", Robj::from(self.intern.ga.keep_all_generations)),
        ]);

        // Convert CV fields
        let cv = List::from_pairs(vec![
            ("fold_number", Robj::from(self.intern.cv.fold_number)),
            ("overfit_penalty", Robj::from(self.intern.cv.overfit_penalty)),
        ]);

        // Combine all sections into a single R list object
        let param_list = List::from_pairs(vec![
            ("general", Robj::from(general)),
            ("data", Robj::from(data)),
            ("ga", Robj::from(ga)),
            ("cv", Robj::from(cv)),
        ]);

        Robj::from(param_list)
    }   

    /// Set the minimal prevalence of feature for feature selection
    /// @export 
    pub fn set_feature_minimal_prevalence_pct(&mut self, pct: f64) {
        self.intern.data.feature_minimal_prevalence_pct = pct;
    }

    /// Set the log level (possible value trace, debug, info, warning, error)
    /// @export 
    pub fn set_log_level(&mut self, level: String) {
        self.intern.general.log_level = level;
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
    param: Param,
    train_data: GData,
    test_data: GData,
    generations: Vec<GPopulation>
}

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
    pub fn get_individual_full(&self, generation: i32, order:i32, verbose:bool) -> Robj {

        let mut features: Vec<String> = Vec::new();
        let mut value = Vec::new();
        let individual = self.generations[generation as usize].individuals[order as usize].clone();
        for (k, v) in individual.features.iter() {
            features.push(self.train_data.features[*k].clone());      // R strings
            value.push(*v as i32);     // R integers are 32-bit
        }
        
        // Create an integer vector from 'vals'
        let value_robj: Robj = Robj::from(value);
        let features_robj: Robj = Robj::from(features);
        let k_robj: Robj = Robj::from(individual.k as i32);
        let auc_robj: Robj = Robj::from(individual.auc);
        let n_robj: Robj = Robj::from(individual.n as i32);
        let fit_robj: Robj = Robj::from(self.generations[generation as usize].fit[order as usize]);

        let individual = List::from_pairs(vec![
            ("features", features_robj),
            ("value", value_robj),
            ("k", k_robj),
            ("auc", auc_robj),
            ("n", n_robj),
            ("fit", fit_robj)
        ]).into_robj();

        if verbose {
            println!("individual: {:?}",&individual);
        }
        individual

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
        let generation_index = generation as usize;
        let individuals_count = self.generations[generation_index].individuals.len();
        let mut individuals_list = Vec::new();

        for order in 0..individuals_count {
            let individual_robj = self.get_individual_full(generation, order as i32, false);
            individuals_list.push(individual_robj);
        }

        let generation_robj = Robj::from(individuals_list);

        generation_robj
    }



    /// list an individual at generation #generation, order #order with the number
    /// of generation available with generation_number and the number of individuals
    /// of a certain generation with population_size
    /// @export
    pub fn get_individual(&self, generation: i32, order:i32) -> Robj {
        
        let mut features: Vec<String> = Vec::new();
        let mut value = Vec::new();
        for (k, v) in self.generations[generation as usize].individuals[order as usize].features.iter() {
            features.push(self.train_data.features[*k].clone());      // R strings
            value.push(*v as i32);     // R integers are 32-bit
        }

        //println!("F {:?} V {:?}", features, value);
        // Create an integer vector from 'vals'
        let mut robj = Robj::from(value);

        // Assign names to the vector
        // (set_names requires an iterator of &str, so pass &keys[..])
        robj.set_names(features).ok(); // ignore error for simplicity

        robj
    }

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


    /// get the number of generation included in the Population object
    /// @export
    pub fn generation_number(&self) -> i32 {
        self.generations.len() as i32
    }

    /// get the size (number of individuals) of a certain generation in a Population
    /// @export
    pub fn population_size(&self, generation: i32) -> i32 {
        self.generations[generation as usize].individuals.len() as i32
    }

    /// get an individual auc
    /// @export
    pub fn get_individual_train_auc(&self, generation: i32, order:i32) -> f64 {
        self.generations[generation as usize].individuals[order as usize].auc
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
/// Data object
///////////////////////////////////////////////////////////////

/*
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
    pub fn get_data(&self) -> Robj {
        // Convert Data fields to R objects
        let data = List::from_pairs(vec![
            ("X", Robj::from(&self.intern.X)),
            ("y", Robj::from(&self.intern.y)),
            ("Xtest", Robj::from(&self.intern.Xtest)),
            ("ytest", Robj::from(&self.intern.ytest)),
            ("pvalue_method", Robj::from(&self.intern.pvalue_method)),
            ("feature_minimal_prevalence_pct", Robj::from(self.intern.feature_minimal_prevalence_pct)),
            ("feature_maximal_pvalue", Robj::from(self.intern.feature_maximal_pvalue)),
            ("feature_minimal_feature_value", Robj::from(self.intern.feature_minimal_feature_value)),
        ]);

        // Return the data as an R list object
        Robj::from(data)
    }

}
 */

///////////////////////////////////////////////////////////////
/// Glogger object
///////////////////////////////////////////////////////////////

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
        println!("You can only set a logger once");
        Self {
            handle: Logger::try_with_str("info") // Set the log level (e.g., "info")
                        .unwrap()
                        .write_mode(WriteMode::BufferAndFlush) // Use buffering for smoother output
                        .start() // Start the logger
                        .unwrap_or_else(|e| panic!("Logger initialization failed with {}", e))
        }
    }

    /// Create a new screen logger
    /// @export
    pub fn level(level: String) -> Self {
        println!("You can only set a logger once");
        Self {
            handle: Logger::try_with_str(level) // Set the log level (e.g., "info")
                        .unwrap()
                        .write_mode(WriteMode::BufferAndFlush) // Use buffering for smoother output
                        .start() // Start the logger
                        .unwrap_or_else(|e| panic!("Logger initialization failed with {}", e))
        }
    }

    /// Create a new logger from a Param
    /// @export  
    pub fn get(param: &Param) -> Self {
        println!("You can only set a logger once");
        let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S").to_string();

        // Initialize the logger
        let handle = if param.intern.general.log_base.len()>0 {
            Logger::try_with_str(&param.intern.general.log_level) // Set log level (e.g., "info")
                .unwrap()
                .log_to_file(
                    FileSpec::default()
                        .basename(&param.intern.general.log_base) // Logs go into the "logs" directory
                        .suffix(&param.intern.general.log_suffix)     // Use the ".log" file extension
                        .discriminant(&timestamp), // Add timestamp to each log file
                )
                .write_mode(WriteMode::BufferAndFlush) // Control file write buffering
                .format_for_files(custom_format) // Custom format for the log file
                .format_for_stderr(custom_format) // Same format for the console
                .start()
                .unwrap_or_else(|e| panic!("Logger initialization failed with {}", e))
        }
        else {
            Logger::try_with_str(&param.intern.general.log_level) // Set the log level (e.g., "info")
                .unwrap()
                .write_mode(WriteMode::BufferAndFlush) // Use buffering for smoother output
                .start() // Start the logger
                .unwrap_or_else(|e| panic!("Logger initialization failed with {}", e))
        };

        Self {
            handle: handle
        }

    }

    /// Change logging level
    /// @export
    pub fn set_level(&mut self, level: String) {
        let new_spec = LogSpecification::parse(level).unwrap();
        self.handle.set_new_spec(new_spec);
    }
}

/// The simple genetic algorithm (ga) produce a Population from a Param object
/// the RunningFlag object is convenient when launching ga in a subthread, it must be
/// provided (but you can let it live its own way)
/// @export
#[extendr]
pub fn ga(param: &Param, running_flag: &RunningFlag) -> Experiment {
    let algo = &param.intern.general.algo;
    let (generations, train_data, test_data) = 
        if (algo=="ga")||(algo=="ga2") { ga_run(&param.intern, running_flag.get_arc()) }
        else  { if (algo=="ga_no_overfit")||(algo=="ga2_no_overfit") {
                ga_no_overfit(&param.intern, running_flag.get_arc())
              } else { panic!("No such algo {}",algo) } };

    Experiment {
        param: param.clone(),
        train_data: train_data,
        test_data: test_data,
        generations: generations
    }
}


// Macro to expose the struct and methods to R
extendr_module! {
    mod gpredomicsR;
    impl RunningFlag;
    impl Experiment;
    impl Param;
    impl GLogger;
    fn ga;
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