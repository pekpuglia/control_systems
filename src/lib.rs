use ode_solvers::{self, System, DVector};



trait DynamicalSystem {
    const STATE_VECTOR_SIZE: usize;
    type Input;
    type Output;

    fn xdot(&self, t: f64, x: DVector<f64>, u: Self::Input) -> DVector<f64>;
    fn y(&self, t: f64, x: DVector<f64>, u: Self::Input) -> Self::Output;
}

struct Series<DS1, DS2> {
    dynsys1: DS1,
    dynsys2: DS2
}

impl<I:Copy, M, O, DS1: DynamicalSystem<Input = I, Output = M>, DS2: DynamicalSystem<Input = M, Output = O>> DynamicalSystem for Series<DS1, DS2> {
    type Input = I;

    type Output = O;

    fn xdot(&self, t: f64, x: DVector<f64>, u: Self::Input) -> DVector<f64> {
        let mut x1dot = self.dynsys1.xdot(
            t, 
            x.rows(0, DS1::STATE_VECTOR_SIZE).into(), 
            u
        );

        x1dot.extend(
            self.dynsys2.xdot(
                t, 
                x.rows(
                    DS1::STATE_VECTOR_SIZE, 
                    DS2::STATE_VECTOR_SIZE).into(), 
                self.dynsys1.y(t, x.rows(0, DS1::STATE_VECTOR_SIZE).into(), 
                u)).into_iter().copied()
        );
        x1dot
    }

    fn y(&self, t: f64, x: DVector<f64>, u: Self::Input) -> Self::Output {
        todo!()
    }

    const STATE_VECTOR_SIZE: usize = DS1::STATE_VECTOR_SIZE + DS2::STATE_VECTOR_SIZE;
}

#[cfg(test)]
mod tests {
    use super::*;

    
}
