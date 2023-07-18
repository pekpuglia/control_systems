use crate::*;
struct Series<DS1, DS2, const SVS1: usize, const SVS2: usize, const OS1IS2: usize> {
    dynsys1: DS1,
    dynsys2: DS2
}

impl<const SVS1: usize, 
    const SVS2: usize, 
    const IS1: usize, 
    const OS1IS2: usize, 
    const OS2: usize, 
    DS1: DynamicalSystem<SVS1, IS1, OS1IS2>, 
    DS2: DynamicalSystem<SVS2, OS1IS2, OS2>> DynamicalSystem<{SVS1+SVS2}, IS1, OS2> for Series<DS1, DS2, SVS1, SVS2, OS1IS2> 
{
    
    fn xdot(&self, t: f64, x: SVector<f64, {SVS1+SVS2}>, u: SVector<f64, IS1>) -> SVector<f64, {SVS1+SVS2}> {
        let x1dot = self.dynsys1.xdot(
            t, 
            SVector::from_row_slice(x.fixed_rows::<SVS1>(0).data.into_slice()),
            u
        );

        let x2dot = self.dynsys2.xdot(
            t, 
            SVector::from_row_slice(x.fixed_rows::<SVS2>(SVS1).data.into_slice()), 
            self.dynsys1.y(t, SVector::from_row_slice(x.fixed_rows::<SVS1>(0).data.into_slice()), u));
        
        //não reclama por algum motivo
        SVector::from_iterator(
            x1dot
            .into_iter()
            .chain(x2dot.into_iter())
            .copied()
        )
    }

    fn y(&self, t: f64, x: SVector<f64, {SVS1 + SVS2}>, u: SVector<f64, IS1>) -> SVector<f64, OS2> {
        self.dynsys2.y(
            t, 
            SVector::from_row_slice(x.fixed_rows::<SVS2>(SVS1).data.into_slice()), 
            self.dynsys1.y(t, SVector::from_row_slice(x.fixed_rows::<SVS1>(0).data.into_slice()), u))
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;
    struct Exp {
        alpha: f64
    }

    impl DynamicalSystem<1, 1, 1> for Exp {
        fn xdot(&self, t: f64, x: SVector<f64, 1>, u: SVector<f64, 1>) -> SVector<f64, 1> {
            self.alpha * (u - x)
        }

        fn y(&self, t: f64, x: SVector<f64, 1>, u: SVector<f64, 1>) -> SVector<f64, 1> {
            x
        }
    }

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

    #[test]
    fn test_series_xdot() {
        let series = Series {
            dynsys1: Exp{alpha: 0.5},
            dynsys2: SecondOrder{ k: 1.0, c: 1.0 }
        };

        let xdot = series.xdot(0.0, [0.0,0.0, 1.0].into(), [1.0].into());
        assert!(xdot == DVector::from_row_slice(&[0.5, 1.0, -1.0]))
    }

    #[test]
    fn test_series_output() {
        let series = Series {
            dynsys1: Exp{alpha: 0.5},
            dynsys2: SecondOrder{ k: 1.0, c: 1.0 }
        };
        let out = series.y(0.0, [0.5, 0.7, 1.0].into(), [1.0].into());
        assert!(out == [0.7].into())
    }
}
