
use extendr_api::prelude::*;
//use gpredomics::param::Param;
use gpredomics::param::Param as GParam;
use gpredomics::data::Data as GData;
use gpredomics::param::get as GParam_get;

use gpredomics::population::Population  as GPopulation;
use gpredomics::ga_run;

use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use extendr_api::prelude::*;

/// A struct to manage the `running` flag.
#[derive(Debug, Clone)]
#[extendr]
pub struct RunningFlag {
    flag: Arc<AtomicBool>,
}

#[extendr]
impl RunningFlag {
    /// Create a new `RunningFlag` (initially set to `true`).
    pub fn new() -> Self {
        Self {
            flag: Arc::new(AtomicBool::new(true)),
        }
    }

    /// Set the `running` flag to `false`.
    pub fn stop(&self) {
        self.flag.store(false, Ordering::Relaxed);
    }

    /// Check the current value of the `running` flag.
    pub fn is_running(&self) -> bool {
        self.flag.load(Ordering::Relaxed)
    }

    /// Reset the `running` flag to `true`.
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


#[extendr]
pub struct Param {
    intern: GParam
}

#[extendr]
impl Param {
    pub fn new() -> Self {
        Self {
            intern: GParam::new()
        }
    }

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

    pub fn set_feature_minimal_prevalence_pct(&mut self, pct: f64) {
        self.intern.data.feature_minimal_prevalence_pct = pct;
    }

}

#[extendr]
pub struct Population {
    generations: Vec<GPopulation>,
    train_data: GData,
    test_data: GData
}

#[extendr]
impl Population {
    pub fn get_individuals(&self, generation: i32, order:i32) -> Robj {
        
        let mut features: Vec<String> = Vec::new();
        let mut value = Vec::new();
        for (k, v) in self.generations[generation as usize].individuals[order as usize].features.iter() {
            features.push(self.train_data.features[*k].clone());      // R strings
            value.push(*v as i32);     // R integers are 32-bit
        }

        // Create an integer vector from 'vals'
        let mut robj = Robj::from(value);

        // Assign names to the vector
        // (set_names requires an iterator of &str, so pass &keys[..])
        robj.set_names(features).ok(); // ignore error for simplicity

        robj
    }

    pub fn generation_number(&self) -> i32 {
        self.generations.len() as i32
    }

    pub fn population_size(&self, generation: i32) -> i32 {
        self.generations[generation as usize].individuals.len() as i32
    }
}


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