
use extendr_api::prelude::*;
//use gpredomics::param::Param;
use gpredomics::param::Param as GParam;
use gpredomics::data::Data as GData;
use gpredomics::param::get as GParam_get;

use gpredomics::population::Population  as GPopulation;
use gpredomics::ga_run;

use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use flexi_logger::{Logger, WriteMode, FileSpec};
use chrono::Local;

use extendr_api::prelude::*;

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


/// gpredomics param object that create all settings
/// @export 
#[extendr]
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
    pub fn get(file_path: String) -> Self {
        if let Ok(this_param) = GParam_get(file_path.clone()) {
            Self {
                intern: this_param
            }
        }
        else {
            panic!("Unsuitable param file: {}",file_path)
        }
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

/// A global Population object that proxies all the different Rust gpredomics objects under
/// the hood
/// @export
#[extendr]
pub struct Population {
    generations: Vec<GPopulation>,
    train_data: GData,
    test_data: GData
}

/// @export
#[extendr]
impl Population {

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
        log::set_max_level(
            match &level as &str {
                "trace" => log::LevelFilter::Trace,
                "debug" => log::LevelFilter::Debug,
                "info" => log::LevelFilter::Info,
                "warning" => log::LevelFilter::Warn,
                "error" => log::LevelFilter::Error,

                other => panic!("No such level {}", other)
            }
        );
    }
}

/// The simple genetic algorithm (ga) that produce a Population from a Param object
/// the RunningFlag object is convenient when launching ga in a subthread, it must be
/// provided (but you can let it live its own way)
/// @export
#[extendr]
pub fn ga(param: &Param, running_flag: &RunningFlag) -> Population {
    let (generations, train_data, test_data) = ga_run(&param.intern, running_flag.get_arc());

    Population {
        generations: generations,
        train_data: train_data,
        test_data: test_data
    }
}

//fn get_pop() -> Population {
//    Population::new()
//}

// Macro to expose the struct and methods to R
extendr_module! {
    mod gpredomicsR;
    impl RunningFlag;
    impl Population;
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