use std::marker::PhantomData;

pub use nalgebra::{DVector, dvector};

//todo
//remove cloning?
//make it size-safe!
//system lib?
//state vector builder
//receive only statevectors, inputvectors and outputvectors?
pub trait DynamicalSystem {

    const STATE_VECTOR_SIZE: usize;
    const INPUT_SIZE      : usize;
    const OUTPUT_SIZE     : usize;
    //accept references or StateVectors!!!
    fn xdot(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64>;
    fn y(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64>;
}

#[derive(Clone, Debug)]
pub struct StateVector<System: DynamicalSystem> {
    pub data: DVector<f64>,
    _phantom: PhantomData<System>
}

impl<System: DynamicalSystem> StateVector<System>  {
    pub fn new(x: DVector<f64>) -> Self {
        assert!(x.len() == System::STATE_VECTOR_SIZE);
        StateVector { data: x, _phantom: PhantomData }
    }
}

impl<System: DynamicalSystem> From<DVector<f64>> for StateVector<System> {
    fn from(value: DVector<f64>) -> Self {
        StateVector { data: value, _phantom: PhantomData }
    }
}

impl<System: DynamicalSystem> From<&[f64]> for StateVector<System>  {
    fn from(value: &[f64]) -> Self {
        assert!(value.len() == System::STATE_VECTOR_SIZE);
        StateVector { data: DVector::from_row_slice(value), 
            _phantom: PhantomData }
    }
}

mod series;
pub use series::Series;


mod feedback;
pub use feedback::{NegativeFeedback, UnitySystem, UnityFeedback};

mod parallel;
pub use parallel::Parallel;