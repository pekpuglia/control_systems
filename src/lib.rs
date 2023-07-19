#![feature(generic_const_exprs)]
pub use nalgebra::SVector;

//problema: pode ter mais de 1 implementação por tipo p parâmetros genéricos diferentes
pub trait DynamicalSystem<
    const STATE_VECTOR_SIZE: usize,
    const INPUT_SIZE      : usize,
    const OUTPUT_SIZE     : usize> {

    fn xdot(&self, t: f64, 
        x: SVector<f64, STATE_VECTOR_SIZE>, 
        u: SVector<f64, INPUT_SIZE>) -> SVector<f64, STATE_VECTOR_SIZE>;
    fn y(&self, t: f64, 
        x: SVector<f64, STATE_VECTOR_SIZE>, 
        u: SVector<f64, INPUT_SIZE>) -> SVector<f64, OUTPUT_SIZE>;
}

mod series;
pub use series::Series;


mod feedback;
pub use feedback::NegativeFeedback;