
use crate::*;
pub struct Series<DS1, DS2>
{
    dynsys1: DS1,
    dynsys2: DS2
}

// impl<DS1, DS2, const SVS1: usize, const SVS2: usize, const OS1IS2: usize> Series<DS1, DS2, SVS1, SVS2, OS1IS2> {
//     pub fn new(dynsys1: DS1, dynsys2: DS2) -> Series<DS1, DS2, SVS1, SVS2, OS1IS2> {
//         Series { dynsys1, dynsys2 }
//     }
// }

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
            x.rows(0, DS1::STATE_VECTOR_SIZE).into(), 
            u.clone()
        );
        let x2dot = self.dynsys2.xdot(
            t, 
            x.rows(DS1::STATE_VECTOR_SIZE, DS2::STATE_VECTOR_SIZE).into(), 
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
            x.rows(0, DS1::STATE_VECTOR_SIZE).into(), 
            u
        );

        self.dynsys2.y(
            t, 
            x.rows(DS1::STATE_VECTOR_SIZE, DS2::STATE_VECTOR_SIZE).into(), 
            u1)
    }
}

// impl<const SVS1: usize, 
//     const SVS2: usize, 
//     const IS1: usize, 
//     const OS1IS2: usize, 
//     const OS2: usize, 
//     DS1: DynamicalSystem, 
//     DS2: DynamicalSystem> DynamicalSystem for Series<DS1, DS2, SVS1, SVS2, OS1IS2> 
// {
    
//     fn xdot(&self, t: f64, x: DVector<f64>, u: DVector<f64>) -> DVector<f64> {
//         let x1dot = self.dynsys1.xdot(
//             t, 
//             x.rows(0, 1),
//             u
//         );

//         let x2dot = self.dynsys2.xdot(
//             t, 
//             x.rows(0, 1), 
//             self.dynsys1.y(t, x.rows, u));
        
//         //não reclama por algum motivo
//         SVector::from_iterator(
//             x1dot
//             .into_iter()
//             .chain(x2dot.into_iter())
//             .copied()
//         )
//     }

//     fn y(&self, t: f64, x: SVector<f64, {SVS1 + SVS2}>, u: SVector<f64, IS1>) -> SVector<f64, OS2> {
//         self.dynsys2.y(
//             t, 
//             SVector::from_row_slice(x.fixed_rows::<SVS2>(SVS1).data.into_slice()), 
//             self.dynsys1.y(t, SVector::from_row_slice(x.fixed_rows::<SVS1>(0).data.into_slice()), u))
//     }

// }

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
}
