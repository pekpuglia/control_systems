pub use nalgebra::{DVector, dvector};

//todo
//remove cloning?
//make it size-safe!
//system lib?
//receive only inputvectors and outputvectors?

mod state_vector;
pub use state_vector::ComposableVector;

pub trait DynamicalSystem<IN, ST, OUT> 
where
    IN: ComposableVector,
    ST: ComposableVector,
    OUT: ComposableVector
{

    //accept references or StateVectors!!!
    fn xdot(&self, t: f64, 
        x: &ST, 
        u: &IN) -> ST;
    fn y(&self, t: f64, 
        x: &ST, 
        u: &IN) -> OUT;
}

mod series;
pub use series::Series;


// mod feedback;
// pub use feedback::{NegativeFeedback, UnitySystem, UnityFeedback};

// mod parallel;
// pub use parallel::Parallel;