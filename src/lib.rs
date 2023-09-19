use std::marker::PhantomData;

pub use nalgebra::{DVector, dvector};

//todo
//remove cloning?
//make it size-safe!
//system lib?
//state vector builder
//receive only statevectors, inputvectors and outputvectors?
//index for statevectors

mod state_vector;
pub use state_vector::{StateVector, IntoSV};

pub trait DynamicalSystem {

    const STATE_VECTOR_SIZE: usize;
    const INPUT_SIZE      : usize;
    const OUTPUT_SIZE     : usize;
    
    //accept references or StateVectors!!!
    fn xdot(&self, t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64>;
    fn y(&self, t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64>;
}

mod series;
pub use series::Series;


mod feedback;
pub use feedback::{NegativeFeedback, UnitySystem, UnityFeedback};

mod parallel;
pub use parallel::Parallel;