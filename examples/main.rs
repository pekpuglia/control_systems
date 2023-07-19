use control_systems::{DynamicalSystem, SVector};

struct SecondOrder {
    k: f64,
    c: f64
}

impl DynamicalSystem<2, 1, 1> for SecondOrder {
    fn xdot(&self, _t: f64, 
        x: SVector<f64, 2>, 
        u: SVector<f64, 1>) -> SVector<f64, 2> {
        [
            x.y,
            -self.k*x.x-self.c*x.y+u.x
        ].into()
    }

    fn y(&self, _t: f64, 
        x: SVector<f64, 2>, 
        _u: SVector<f64, 1>) -> SVector<f64, 1> {
        SVector::<f64, 1>::new(x.x)
    }
}

fn main() {
    
}