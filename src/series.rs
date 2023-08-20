
use std::marker::PhantomData;

use crate::*;
#[derive(Clone, Copy, Debug)]
pub struct Series<DS1, DS2, StateVectorEnum>
{
    dynsys1: DS1,
    dynsys2: DS2,
    _phantom: PhantomData<StateVectorEnum>
}

impl<DS1, DS2, StateVectorEnum> Series<DS1, DS2, StateVectorEnum>  {
    pub fn new(ds1: DS1, ds2: DS2) -> Self {
        Series { dynsys1: ds1, dynsys2: ds2, _phantom: PhantomData }
    }
    pub fn ds1_ref(&self) -> &DS1 {
        &self.dynsys1
    }

    pub fn ds2_ref(&self) -> &DS2 {
        &self.dynsys2
    }
}

use paste::paste;

macro_rules! sum_of_enums {
    ($base_name: ident = $($factors: ty),+) => {
        paste!{
            #[derive(Enum)]
            enum [<$base_name $($factors) +>] {
                $(
                    $factors($factors)
                ),+
            }
        }
    };
}

impl<DS1, DS2, StateVectorEnum: EnumArray<f64>> DynamicalSystem for Series<DS1, DS2, StateVectorEnum> 
where
    DS1: DynamicalSystem,
    DS2: DynamicalSystem
{
    type Input = DS1::Input;

    type StateVector = StateVectorEnum;

    type Output = DS2::Output;

    fn xdot(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
            todo!()
    }

    fn y(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
        todo!()
    }

    // fn xdot(&self, t: f64, 
    //     x: DVector<f64>, 
    //     u: DVector<f64>) -> DVector<f64> {
    //     let mut x1dot = self.dynsys1.xdot(
    //         t, 
    //         x.rows(0, DS1::STATE_VECTOR_SIZE).into(), 
    //         u.clone()
    //     );
    //     let x2dot = self.dynsys2.xdot(
    //         t, 
    //         x.rows(DS1::STATE_VECTOR_SIZE, DS2::STATE_VECTOR_SIZE).into(), 
    //         self.dynsys1.y(
    //             t, x.rows(0, DS1::STATE_VECTOR_SIZE).into(), u)
    //     );

    //     x1dot
    //         .extend(x2dot.iter().copied());
    //     x1dot
    // }

    // fn y(&self, t: f64, 
    //     x: DVector<f64>, 
    //     u: DVector<f64>) -> DVector<f64> {
    //     let u1 = self.dynsys1.y(
    //         t, 
    //         x.rows(0, DS1::STATE_VECTOR_SIZE).into(), 
    //         u
    //     );

    //     self.dynsys2.y(
    //         t, 
    //         x.rows(DS1::STATE_VECTOR_SIZE, DS2::STATE_VECTOR_SIZE).into(), 
    //         u1)
    // }
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
    use enum_map::enum_map;
    struct ExpSys {
        alpha: f64
    }

    #[derive(Enum)]
    enum Input {
        Reference
    }

    #[derive(Enum)]
    enum State {
        Position
    }

    #[derive(Enum)]
    enum Output {
        Position
    }

    impl DynamicalSystem for ExpSys {
        type Input = Input;

        type StateVector = State;

        type Output = Output;

        fn xdot(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
            enum_map!(
                State::Position => self.alpha * (u[Input::Reference] - x[State::Position])
            )
        }

        fn y(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
            enum_map! {
                Output::Position => x[State::Position]
            }
        }
    }

    struct SecondOrder {
        k: f64,
        c: f64
    }

    #[derive(Enum)]
    enum SOInp {
        Force
    }

    #[derive(Enum)]
    enum SOSV {
        Position,
        Velocity
    }

    #[derive(Enum)]
    enum SOOut {
        Position
    }

    impl DynamicalSystem for SecondOrder {
        type Input = SOInp;

        type StateVector = SOSV;

        type Output = SOOut;

        fn xdot(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
            enum_map!(
                SOSV::Position => x[SOSV::Velocity],
                SOSV::Velocity => -self.k*x[SOSV::Position]-self.c*x[SOSV::Velocity]+u[SOInp::Force]
            )
        }

        fn y(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
            enum_map! {SOOut::Position => x[SOSV::Position]}
        }
        // fn xdot(&self, _t: f64, 
        //     x: DVector<f64>, 
        //     u: DVector<f64>) -> DVector<f64> {
            
        //     dvector![
        //         x[1],
        //         -self.k*x[0]-self.c*x[1]+u[0]
        //     ]
        // }

        // fn y(&self, _t: f64, 
        //     x: DVector<f64>, 
        //     _u: DVector<f64>) -> DVector<f64> {
        //     dvector![x[0]]
        // }

        // const STATE_VECTOR_SIZE: usize = 2;

        // const INPUT_SIZE      : usize = 1;

        // const OUTPUT_SIZE     : usize = 1;
    }

    #[test]
    fn test_series_xdot() {
        let series = Series::new(ExpSys{alpha: 0.5}, SecondOrder{ k: 1.0, c: 1.0 });
        
        

        let xdot = series.xdot(0.0, dvector![0.0,0.0, 1.0], dvector![1.0]);
        dbg!(&xdot);
        assert!(xdot == DVector::from_row_slice(&[0.5, 1.0, -1.0]))
    }

    #[test]
    fn test_series_output() {
        let series = Series::new(ExpSys{alpha: 0.5}, SecondOrder{ k: 1.0, c: 1.0 });
        let out = series.y(0.0, dvector![0.5, 0.7, 1.0], dvector![1.0]);
        assert!(out == [0.7].into())
    }
}