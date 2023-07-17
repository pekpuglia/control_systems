use ode_solvers::{self, System, DVector};

trait DynamicalSystem {
    const STATE_VECTOR_SIZE: usize;
    type Input;
    type Output;

    fn xdot(&self, t: f64, x: DVector<f64>, u: Self::Input) -> DVector<f64>;
    fn y(&self, t: f64, x: DVector<f64>, u: Self::Input) -> Self::Output;
}

mod series;
