use core::panic;
use std::marker::PhantomData;

use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;
use nalgebra::dvector;

// #[derive(Clone, Copy, Debug)]
// pub struct UnitySystem<const SIZE: usize>;

// impl<const SIZE: usize> DynamicalSystem for UnitySystem<SIZE> {
//     fn xdot(&self, _t: f64, 
//         _x: DVector<f64>, 
//         _u: DVector<f64>) -> DVector<f64> {
//         dvector![]
//     }

//     fn y(&self, _t: f64, 
//         _x: DVector<f64>, 
//         u: DVector<f64>) -> DVector<f64> {
//         u
//     }

//     const STATE_VECTOR_SIZE: usize = 0;

//     const INPUT_SIZE      : usize = SIZE;

//     const OUTPUT_SIZE     : usize = SIZE;
// }

use enum_iterator::{Sequence, all};

pub trait AsMap<T: EnumArray<f64> + Sequence + Default> {
    fn as_map(&self) -> EnumMap<T, f64>;
}

impl<T: EnumArray<f64> + Sequence + Default> AsMap<T> for DVector<f64> {
    fn as_map(&self) -> EnumMap<T, f64> {
        EnumMap::from_iter(self.iter().zip(all::<T>()).map(|(num, key)| (key, *num) ))
    }
}

pub trait AsVector {
    fn as_vector(&self) -> DVector<f64>;
}

impl<T: EnumArray<f64>> AsVector for EnumMap<T, f64> {
    fn as_vector(&self) -> DVector<f64> {
        DVector::from_row_slice(self.as_slice())
    }
}

#[derive(Debug, Clone, Copy)]
//direct dynamical system, reverse dynamical system
pub struct NegativeFeedback<DDS, RDS, StateVectorEnum> 
{
    dirsys: DDS,
    revsys: RDS,
    _phantom: PhantomData<StateVectorEnum>
}

impl<DDS, RDS, StateVectorEnum> 
    NegativeFeedback<DDS, RDS, StateVectorEnum>
where
    DDS: DynamicalSystem,
    RDS: DynamicalSystem,
    StateVectorEnum: EnumArray<f64>,
    RDS::Input: Default + Sequence,
    DDS::Input: Default + Sequence,
    StateVectorEnum: Default + Sequence,
    DVector<f64>: AsMap<StateVectorEnum> + AsMap<RDS::Input> + AsMap<DDS::Input>,
    EnumMap<DDS::Output, f64>: AsVector,
    EnumMap<RDS::Output, f64>: AsVector,
    EnumMap<StateVectorEnum, f64>: SubMapOps<First = DDS::StateVector, Second = RDS::StateVector>
{
    fn residue(&self, t: f64, y: DVector<f64>, x: DVector<f64>, u: DVector<f64>) -> DVector<f64> {
        y.clone() - self.dirsys.y(
            t, 
            x.as_map().first(), 
            (u - self.revsys.y(
                t, 
                x.as_map().second(), 
                y.as_map()).as_vector()).as_map()
            ).as_vector()
    }

    pub fn new(dirsys: DDS, revsys: RDS) -> NegativeFeedback<DDS, RDS, StateVectorEnum> {
        // NegativeFeedback::<DDS, RDS>::assert_sizes();
        NegativeFeedback { dirsys, revsys, _phantom: PhantomData }
    }

    // const fn assert_sizes() {
    //     if DDS::OUTPUT_SIZE != RDS::INPUT_SIZE || DDS::INPUT_SIZE != RDS::OUTPUT_SIZE {
    //         panic!("Wrong sizes!")
    //     }
    // }

    pub fn dir_ref(&self) -> &DDS {
        &self.dirsys
    }

    pub fn rev_ref(&self) -> &RDS {
        &self.revsys
    }
}

// pub type UnityFeedback<DDS, const SIZE: usize> = NegativeFeedback<DDS, UnitySystem<SIZE>>;

// impl<DDS: DynamicalSystem, const SIZE: usize> NegativeFeedback<DDS, UnitySystem<SIZE>> {
//     pub fn new_unity_feedback(dirsys: DDS) -> NegativeFeedback<DDS, UnitySystem<SIZE>>
//     {
//         NegativeFeedback::new(dirsys, UnitySystem{})
//     }
// }

impl<DDS, RDS, StateVectorEnum: EnumArray<f64>> DynamicalSystem for 
    NegativeFeedback<DDS, RDS, StateVectorEnum> 
where
    DDS: DynamicalSystem,
    RDS: DynamicalSystem,
    StateVectorEnum: EnumArray<f64>,
    RDS::Input: Default + Sequence,
    DDS::Input: Default + Sequence,
    DDS::Output: Default + Sequence,
    StateVectorEnum: Default + Sequence,
    DVector<f64>: AsMap<StateVectorEnum> + AsMap<RDS::Input> + AsMap<DDS::Input> + AsMap<DDS::Output>,
    EnumMap<DDS::Input, f64>: Copy,
    EnumMap<DDS::Output, f64>: AsVector + ConvertMap<RDS::Input>,
    EnumMap<RDS::Output, f64>: AsVector,
    EnumMap<StateVectorEnum, f64>: SubMapOps<First = DDS::StateVector, Second = RDS::StateVector> + Copy
{
    type Input = DDS::Input;

    type StateVector = StateVectorEnum;

    type Output = DDS::Output;

    fn xdot(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
        let output = self.y(t, x, u);
        let rev_output = self.revsys.y(t, x.second(), output.convert());

        let dirxdot = self.dirsys.xdot(t, x.first(), (u.as_vector() - rev_output.as_vector()).as_map());
        let revxdot = self.revsys.xdot(t, x.second(), output.convert());

        SubMapOps::merge(dirxdot, revxdot)
    }

    fn y(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
        MultiVarNewtonFD::new(
            |output| self.residue(t, output, x.as_vector(), u.as_vector())
        ).solve(self.dirsys.y(t, x.first(), u).as_vector()).unwrap().as_map()
    }
}

#[cfg(test)]
mod tests {


    use crate::DynamicalSystem;
    use super::*;
    use enum_map::enum_map;

    #[derive(Clone, Copy)]
    struct LinearFunc {
        a: f64
    }

    #[derive(Enum, Clone, Copy, Debug, Default, Sequence)]
    enum LFInp {
        #[default]
        X
    }

    //terrível
    #[derive(Enum, Clone, Copy, Debug, Default, Sequence)]
    enum LFSV {
        #[default]
        Phantom
    }

    #[derive(Enum, Clone, Copy, Debug, Default, Sequence)]
    enum LFOut {
        #[default]
        Y
    }

    impl DynamicalSystem for LinearFunc {
        type Input = LFInp;

        type StateVector = LFSV;

        type Output = LFOut;

        fn xdot(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
            enum_map! {LFSV::Phantom => 0.0}
        }

        fn y(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
            enum_map! {
                LFOut::Y => self.a * u[LFInp::X]
            }
        }
    }

    struct ExpSys {
        alpha: f64
    }

    #[derive(Enum)]
    enum ESInp {
        Reference
    }

    #[derive(Enum)]
    enum ESSV {
        Position
    }

    #[derive(Enum)]
    enum ESOut {
        Error
    }

    impl DynamicalSystem for ExpSys {
        type Input = ESInp;

        type StateVector = ESSV;

        type Output = ESOut;

        fn xdot(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
            enum_map!(
                ESSV::Position => self.alpha * (u[ESInp::Reference] - x[ESSV::Position])
            )
        }

        fn y(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
            enum_map! {
                ESOut::Error => u[ESInp::Reference] - x[ESSV::Position]
            }
        }
    }

    #[test]
    fn test_feedback_output() {
        let sys = LinearFunc{a: 2.0};
        let feedback_sys = NegativeFeedback::new(sys, sys);

        // let output = feedback_sys.y(0.0, dvector![], dvector![1.0]);

        // assert!((output[0] - 0.4).abs() < 1e-15);
    }

    #[test]
    fn test_feedback_xdot() {
        // let exp1 = Exp{ alpha: 1.0 };
        // let exp2 = Exp{ alpha: 1.0 };

        // let feedback_sys = NegativeFeedback { 
        //     dirsys: exp1, revsys: exp2 };

        // let xdot = feedback_sys.xdot(0.0, dvector![1.0, 2.0], dvector![3.0]);

        // assert!(xdot == dvector![0.0, -1.0]);
    }
}