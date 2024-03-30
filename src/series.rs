
use std::marker::PhantomData;

use crate::*;

use self::state_vector::VecConcat;
#[derive(Clone, Copy, Debug)]
pub struct Series<DS1, IN1, ST1, OUT1IN2, DS2, ST2, OUT2>
where
    IN1: ComposableVector,
    ST1: ComposableVector,
    OUT1IN2: ComposableVector,
    ST2: ComposableVector,
    OUT2: ComposableVector,
    DS1: DynamicalSystem<IN1, ST1, OUT1IN2>,
    DS2: DynamicalSystem<OUT1IN2, ST2, OUT2>
{
    dynsys1: DS1,
    dynsys2: DS2,
    _phantom_in1: PhantomData<IN1>,
    _phantom_st1: PhantomData<ST1>,
    _phantom_out1in2: PhantomData<OUT1IN2>,
    _phantom_st2: PhantomData<ST2>,
    _phantom_out2: PhantomData<OUT2>,
}

impl<DS1, IN1, ST1, OUT1IN2, DS2, ST2, OUT2> Series<DS1, IN1, ST1, OUT1IN2, DS2, ST2, OUT2>
where
    IN1: ComposableVector,
    ST1: ComposableVector,
    OUT1IN2: ComposableVector,
    ST2: ComposableVector,
    OUT2: ComposableVector,
    DS1: DynamicalSystem<IN1, ST1, OUT1IN2>,
    DS2: DynamicalSystem<OUT1IN2, ST2, OUT2>
{
    pub fn new(ds1: DS1, ds2: DS2) -> Self {
        Series { dynsys1: ds1, dynsys2: ds2, _phantom_in1: PhantomData, _phantom_st1: PhantomData, _phantom_out1in2: PhantomData, _phantom_st2: PhantomData, _phantom_out2: PhantomData }
    }
    pub fn ds1_ref(&self) -> &DS1 {
        &self.dynsys1
    }

    pub fn ds2_ref(&self) -> &DS2 {
        &self.dynsys2
    }

    pub fn y1(&self, t: f64, x: &VecConcat<ST1, ST2>, u: &IN1) -> OUT1IN2 {
        self.dynsys1.y(
            t, 
            x.first(), 
            u
        )
    }
}

impl<DS1, IN1, ST1, OUT1IN2, DS2, ST2, OUT2> DynamicalSystem<IN1, VecConcat<ST1, ST2>, OUT2> for Series<DS1, IN1, ST1, OUT1IN2, DS2, ST2, OUT2> 
where
    IN1: ComposableVector,
    ST1: ComposableVector,
    OUT1IN2: ComposableVector,
    ST2: ComposableVector,
    OUT2: ComposableVector,
    DS1: DynamicalSystem<IN1, ST1, OUT1IN2>,
    DS2: DynamicalSystem<OUT1IN2, ST2, OUT2>
{
    fn xdot(&self, t: f64, 
        x: &VecConcat<ST1, ST2>, 
        u: &IN1) -> VecConcat<ST1, ST2> {
        self.dynsys1.xdot(t, x.first(), u).concat(
            self.dynsys2.xdot(t, x.second(), 
                &self.dynsys1.y(
                    t, x.first(), u
                )
            )
        )
    }

    fn y(&self, t: f64, 
        x: &VecConcat<ST1, ST2>, 
        u: &IN1) -> OUT2 {
        self.dynsys2.y(t, x.second(), &self.dynsys1.y(t, x.first(), u))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{dvector, vector, DVector, SVector};
    struct Exp {
        alpha: f64
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

    struct SecondOrder {
        k: f64,
        c: f64
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

    #[test]
    fn test_series_xdot() {
        let series = Series::new(
            Exp{alpha: 0.5}, SecondOrder{ k: 1.0, c: 1.0 }
        );

        let xdot = series.xdot(0.0, &vector![0.0].concat(vector![0.0, 1.0]), &vector![1.0]);
        assert!(xdot.first().x == 0.5);
        assert!(xdot.second().x == 1.0);
        assert!(xdot.second().y == -1.0);
    }

    #[test]
    fn test_series_output() {
        let series = Series::new(
            Exp{alpha: 0.5},
            SecondOrder{ k: 1.0, c: 1.0 }
        );

        let out = series.y(0.0, &vector![0.5].concat(vector![0.7, 1.0]), &vector![1.0]);
        assert!(out == [0.7].into())
    }

    #[test]
    fn test_y1() {
        let series = Series::new(
            Exp{alpha: 0.5},
            SecondOrder{ k: 1.0, c: 1.0 }
        );

        let sys1 = Exp{ alpha: 0.5 };

        assert!(series.y1(0.0, &vector![0.5].concat(vector![0.7, 1.0]), &vector![1.0]) == sys1.y(0.0, &vector![0.5], &vector![1.0]))
    }
}
