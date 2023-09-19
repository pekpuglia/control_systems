
use crate::*;
#[derive(Clone, Copy, Debug)]
pub struct Series<DS1, DS2>
{
    dynsys1: DS1,
    dynsys2: DS2
}

impl<DS1: DynamicalSystem, DS2: DynamicalSystem> Series<DS1, DS2>  {
    pub fn new(ds1: DS1, ds2: DS2) -> Self {
        Series { dynsys1: ds1, dynsys2: ds2 }
    }
    pub fn ds1_ref(&self) -> &DS1 {
        &self.dynsys1
    }

    pub fn ds2_ref(&self) -> &DS2 {
        &self.dynsys2
    }

    pub fn y1(&self, t: f64, x: DVector<f64>, u:DVector<f64>) -> DVector<f64> {
        self.dynsys1.y(
            t, 
            StateVector::<Self>::new(x.clone()).x1().data, 
            u
        )
    }
}

impl<DS1, DS2> DynamicalSystem for Series<DS1, DS2> 
where
    DS1: DynamicalSystem,
    DS2: DynamicalSystem
{
    const STATE_VECTOR_SIZE: usize = DS1::STATE_VECTOR_SIZE + DS2::STATE_VECTOR_SIZE;

    const INPUT_SIZE      : usize = DS1::INPUT_SIZE;

    const OUTPUT_SIZE     : usize = DS2::OUTPUT_SIZE;

    fn xdot(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64> {
        let mut x1dot = self.dynsys1.xdot(
            t, 
            StateVector::<Self>::new(x.clone()).x1().data, 
            u.clone()
        );
        let x2dot = self.dynsys2.xdot(
            t, 
            StateVector::<Self>::new(x.clone()).x2().data, 
            self.dynsys1.y(
                t, x.rows(0, DS1::STATE_VECTOR_SIZE).into(), u)
        );

        x1dot
            .extend(x2dot.iter().copied());
        x1dot
    }

    fn y(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64> {
        let u1 = self.dynsys1.y(
            t, 
            StateVector::<Self>::new(x.clone()).x1().data, 
            u
        );

        self.dynsys2.y(
            t, 
            StateVector::<Self>::new(x.clone()).x2().data, 
            u1)
    }
}


impl<DS1: DynamicalSystem, DS2: DynamicalSystem> StateVector<Series<DS1, DS2>> {
    pub fn x1(&self) -> StateVector<DS1> {
        self.data
        .rows(
            0, 
            DS1::STATE_VECTOR_SIZE)
        .as_slice()
        .into()
    }
    pub fn x2(&self) -> StateVector<DS2> {
        self.data
        .rows(
            DS1::STATE_VECTOR_SIZE, 
            DS2::STATE_VECTOR_SIZE)
        .as_slice()
        .into()
    }
}

impl<System: DynamicalSystem> StateVector<System> {
    pub fn series<S2: DynamicalSystem>(&self, sv2: StateVector<S2>) -> StateVector<Series<System, S2>>{
        let dataiter = self.data
        .iter()
        .chain(
            sv2.data
                .iter()
        ).copied();
        StateVector::new(DVector::from_iterator(
                System::STATE_VECTOR_SIZE + S2::STATE_VECTOR_SIZE, 
                dataiter)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DVector, dvector};
    struct Exp {
        alpha: f64
    }

    impl DynamicalSystem for Exp {
        fn xdot(&self, _t: f64, x: DVector<f64>, u: DVector<f64>) -> DVector<f64> {
            self.alpha * (u - x)
        }

        fn y(&self, _t: f64, x: DVector<f64>, _u: DVector<f64>) -> DVector<f64> {
            x
        }

        const STATE_VECTOR_SIZE: usize = 1;

        const INPUT_SIZE      : usize = 1;

        const OUTPUT_SIZE     : usize = 1;
    }

    struct SecondOrder {
        k: f64,
        c: f64
    }

    impl DynamicalSystem for SecondOrder {
        fn xdot(&self, _t: f64, 
            x: DVector<f64>, 
            u: DVector<f64>) -> DVector<f64> {
            
            dvector![
                x[1],
                -self.k*x[0]-self.c*x[1]+u[0]
            ]
        }

        fn y(&self, _t: f64, 
            x: DVector<f64>, 
            _u: DVector<f64>) -> DVector<f64> {
            dvector![x[0]]
        }

        const STATE_VECTOR_SIZE: usize = 2;

        const INPUT_SIZE      : usize = 1;

        const OUTPUT_SIZE     : usize = 1;
    }

    #[test]
    fn test_series_xdot() {
        let series = Series {
            dynsys1: Exp{alpha: 0.5},
            dynsys2: SecondOrder{ k: 1.0, c: 1.0 }
        };

        let xdot = series.xdot(0.0, dvector![0.0,0.0, 1.0], dvector![1.0]);
        dbg!(&xdot);
        assert!(xdot == DVector::from_row_slice(&[0.5, 1.0, -1.0]))
    }

    #[test]
    fn test_series_output() {
        let series = Series {
            dynsys1: Exp{alpha: 0.5},
            dynsys2: SecondOrder{ k: 1.0, c: 1.0 }
        };

        let out = series.y(0.0, dvector![0.5, 0.7, 1.0], dvector![1.0]);
        assert!(out == [0.7].into())
    }

    #[test]
    fn test_state_vector_slicing() {
        let x = dvector![1.0, 2.0, 3.0];
        
        let sv = StateVector::<Series<Exp, SecondOrder>>::new(x);
        
        assert!(sv.x1().data == [1.0].into());
        assert!(sv.x2().data == dvector![2.0, 3.0]);
    }

    #[test]
    fn test_y1() {
        let series = Series {
            dynsys1: Exp{alpha: 0.5},
            dynsys2: SecondOrder{ k: 1.0, c: 1.0 }
        };

        let sys1 = Exp{ alpha: 0.5 };

        assert!(series.y1(0.0, dvector![0.5, 0.7, 1.0], dvector![1.0]) == sys1.y(0.0, dvector![0.5], dvector![1.0]))
    }

    #[test]
    fn test_series_builder() {
        let svseries = [0.5].into_sv::<Exp>()
            .series(
                [0.7, 1.0].into_sv::<SecondOrder>()
        );
        assert!(svseries.data == dvector![0.5, 0.7, 1.0])
    }
}
