use nalgebra::{vector, SVector, Vector1, Vector2};

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

impl DynamicalSystem for LinearFunc {
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
    
    type IN = Vector1<f64>;
    
    type ST = SVector<f64, 0>;
    
    type OUT = Vector1<f64>;
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
        x: &SVector<f64, 1>, 
        u: &SVector<f64, 1>) -> SVector<f64, 1> {
        self.alpha * (u - x)
    }
    
    fn y(&self, t: f64, 
        x: &SVector<f64, 1>, 
        u: &SVector<f64, 1>) -> SVector<f64, 1> {
        *x
    }
    
    type IN = Vector1<f64>;
    
    type ST = Vector1<f64>;
    
    type OUT = Vector1<f64>;
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
    
    type IN = Vector1<f64>;
    
    type ST = Vector2<f64>;
    
    type OUT = Vector1<f64>;
}