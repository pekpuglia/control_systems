use control_systems::{DynamicalSystem, Series, NegativeFeedback, SVector};

struct SecondOrder {
    k: f64,
    c: f64
}

impl DynamicalSystem<2, 1, 1> for SecondOrder {
    fn xdot(&self, t: f64, 
        x: SVector<f64, 2>, 
        u: SVector<f64, 1>) -> SVector<f64, 2> {
        [
            x.y,
            -self.k*x.x-self.c*x.y+u.x
        ].into()
    }

    fn y(&self, t: f64, 
        x: SVector<f64, 2>, 
        u: SVector<f64, 1>) -> SVector<f64, 1> {
        SVector::<f64, 1>::new(x.x)
    }
}

fn main() {
    
}