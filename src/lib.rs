pub use nalgebra::{DVector, dvector};

//todo
//remove cloning?
//make it size-safe! - PROC MACRO CRATE
//system lib?

mod state_vector;
pub use state_vector::{ComposableVector, VecConcat};

pub trait DynamicalSystem
{
    type IN: ComposableVector;
    type ST: ComposableVector;
    type OUT: ComposableVector;

    //accept references or StateVectors!!!
    fn xdot(&self, t: f64, 
        x: &Self::ST, 
        u: &Self::IN) -> Self::ST;
    fn y(&self, t: f64, 
        x: &Self::ST, 
        u: &Self::IN) -> Self::OUT;
}

pub mod systems;


// mod series;
// pub use series::Series;


// mod feedback;
// pub use feedback::{NegativeFeedback, UnitySystem, UnityFeedback};

// mod parallel;
// pub use parallel::Parallel;