#![feature(generic_const_exprs)]
use ode_solvers::{self, DVector, SVector};
use nalgebra::{Vector, Vector2, ArrayStorage, ToConst};

trait DynamicalSystem<
    const StateVectorSize: usize,
    const InputSize      : usize,
    const OutputSize     : usize> {

    fn xdot(&self, t: f64, 
        x: SVector<f64, StateVectorSize>, 
        u: SVector<f64, InputSize>) -> SVector<f64, StateVectorSize>;
    fn y(&self, t: f64, 
        x: SVector<f64, StateVectorSize>, 
        u: SVector<f64, InputSize>) -> SVector<f64, OutputSize>;
}

mod series;
