pub use nalgebra::{DVector, dvector};

//fix warnings
//system lib?

//problema: pode ter mais de 1 implementação por tipo p parâmetros genéricos diferentes
pub trait DynamicalSystem {

    const STATE_VECTOR_SIZE: usize;
    const INPUT_SIZE      : usize;
    const OUTPUT_SIZE     : usize;

    fn xdot(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64>;
    fn y(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64>;
}

mod series;
pub use series::Series;


mod feedback;
pub use feedback::{NegativeFeedback, UnitySystem, UnityFeedback};