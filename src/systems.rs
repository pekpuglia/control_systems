use nalgebra::{dvector, vector, DVector, SVector};

use crate::{state_vector, DynamicalSystem};

#[derive(Clone, Copy)]
pub struct LinearFunc {
    a: f64
}

impl LinearFunc {
    pub const fn new(a: f64) -> LinearFunc {
        LinearFunc { a }
    }
}

impl DynamicalSystem for LinearFunc {
    fn xdot(&self, t: f64, 
        x: &state_vector::StateVector<LinearFunc>, 
        u: DVector<f64>) -> DVector<f64> {
        dvector![]
    }

    fn y(&self, t: f64, 
        x: &state_vector::StateVector<LinearFunc>, 
        u: DVector<f64>) -> DVector<f64> {
        self.a * u
    }
    
    const STATE_VECTOR_SIZE: usize = 0;
    
    const INPUT_SIZE      : usize = 1;
    
    const OUTPUT_SIZE     : usize = 1;
}

pub struct Exp {
    alpha: f64
}

impl Exp {
    pub const fn new(alpha: f64) -> Exp {
        Exp { alpha }
    }
}

impl DynamicalSystem for Exp {
    fn xdot(&self, t: f64, 
        x: &state_vector::StateVector<Exp>, 
        u: DVector<f64>) -> DVector<f64> {
        self.alpha * (u - x.data.clone())
    }
    
    fn y(&self, t: f64, 
        x: &state_vector::StateVector<Exp>, 
        u: DVector<f64>) -> DVector<f64> {
        x.data.clone_owned()
    }
    
    const STATE_VECTOR_SIZE: usize = 1;
    
    const INPUT_SIZE      : usize = 1;
    
    const OUTPUT_SIZE     : usize = 1;
}

pub struct SecondOrder {
    k: f64,
    c: f64
}

impl SecondOrder {
    pub const fn new(k: f64, c: f64) -> SecondOrder {
        SecondOrder { k, c }
    }
}

impl DynamicalSystem for SecondOrder {
    fn xdot(&self, t: f64, 
        x: &state_vector::StateVector<SecondOrder>, 
        u: DVector<f64>) -> DVector<f64> {
        dvector![
            x[1],
            -self.k * x[0] - self.c * x[1]
        ]
    }
    
    fn y(&self, t: f64, 
        x: &state_vector::StateVector<SecondOrder>, 
        u: DVector<f64>) -> DVector<f64> {
        dvector![x[0]]
    }
    
    const STATE_VECTOR_SIZE: usize = 2;
    
    const INPUT_SIZE      : usize = 1;
    
    const OUTPUT_SIZE     : usize = 1;
}