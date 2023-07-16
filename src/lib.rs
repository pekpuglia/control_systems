use ode_solvers::{self, System, DVector};



trait DynamicalSystem {
    type StateVector;
    type Input;
    type Output;

    fn xdot(&self, t: f64, x: Self::StateVector, u: Self::Input) -> Self::StateVector;
    fn y(&self, t: f64, x: Self::StateVector, u: Self::Input) -> Self::Output;
}

struct Series<DS1, DS2> {
    dynsys1: DS1,
    dynsys2: DS2
}

// impl<I, M, O, DS1: DynamicalSystem<Input = I, Output = M>, DS2: DynamicalSystem<Input = M, Output = O>> DynamicalSystem for Series<DS1, DS2> {
//     type StateVector = DVector<f64>;

//     type Input = I;

//     type Output = O;

//     fn xdot(&self, t: f64, x: Self::StateVector, u: Self::Input) -> Self::StateVector {
//         self.dynsys1.xdot(t, x.slice((0,0), (DS1::StateVector::)), u)
//     }

//     fn y(&self, t: f64, x: Self::StateVector, u: Self::Input) -> Self::Output {
//         todo!()
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
}
