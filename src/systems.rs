use nalgebra::{vector, SVector};

use crate::DynamicalSystem;

#[derive(Clone, Copy)]
pub struct LinearFunc {
    a: f64
}

impl LinearFunc {
    pub const fn new(a: f64) -> LinearFunc {
        LinearFunc { a }
    }
}

impl DynamicalSystem<SVector<f64, 1>, SVector<f64, 0>, SVector<f64, 1>> for LinearFunc {
    fn xdot(&self, t: f64, 
        x: &SVector<f64, 0>, 
        u: &SVector<f64, 1>) -> SVector<f64, 0> {
        vector![]
    }

    fn y(&self, t: f64, 
        x: &SVector<f64, 0>, 
        u: &SVector<f64, 1>) -> SVector<f64, 1> {
        self.a * u
    }
}

pub struct Exp {
    alpha: f64
}

impl Exp {
    pub const fn new(alpha: f64) -> Exp {
        Exp { alpha }
    }
}

impl DynamicalSystem<SVector<f64, 1>, SVector<f64, 1>, SVector<f64, 1>> for Exp {
    fn xdot(&self, t: f64, 
        x: &SVector<f64, 1>, 
        u: &SVector<f64, 1>) -> SVector<f64, 1> {
        self.alpha * (u - x)
    }
    
    fn y(&self, t: f64, 
        x: &SVector<f64, 1>, 
        u: &SVector<f64, 1>) -> SVector<f64, 1> {
        *x
    }
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

impl DynamicalSystem<SVector<f64, 1>, SVector<f64, 2>, SVector<f64, 1>> for SecondOrder {
    fn xdot(&self, t: f64, 
        x: &SVector<f64, 2>, 
        u: &SVector<f64, 1>) -> SVector<f64, 2> {
        vector![
            x.y,
            -self.k * x.x - self.c * x.y
        ]
    }
    
    fn y(&self, t: f64, 
        x: &SVector<f64, 2>, 
        u: &SVector<f64, 1>) -> SVector<f64, 1> {
        vector![x.x]
    }
}