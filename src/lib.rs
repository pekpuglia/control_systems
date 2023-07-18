#![feature(generic_const_exprs)]
use nalgebra::SVector;
trait DynamicalSystem<
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
mod feedback;