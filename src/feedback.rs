use std::marker::PhantomData;

use crate::*;

#[derive(Enum)]
pub struct Empty;

impl AsMap<Empty> for DVector<f64> {
    fn as_map(&self) -> EnumMap<Empty, f64> {
        enum_map::enum_map! {
            Empty => 0.0
        }
    }
}

use eqsolver::multivariable::MultiVarNewtonFD;

#[derive(Clone, Copy, Debug)]
pub struct UnitySystem<IO> {
    _phantom: PhantomData<IO>
}

impl<IO> DynamicalSystem for UnitySystem<IO> 
where
    IO: EnumArray<f64>
{
    type Input = IO;

    type StateVector = Empty;

    type Output = IO;

    fn xdot(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
        dvector![].as_map()
    }

    fn y(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
        u
    }
}

pub trait AsMap<T: EnumArray<f64>> {
    fn as_map(&self) -> EnumMap<T, f64>;
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
    DVector<f64>: AsMap<RDS::Input> + AsMap<DDS::Input>,
    EnumMap<DDS::Output, f64>: AsVector,
    EnumMap<RDS::Output, f64>: AsVector,
    EnumMap<StateVectorEnum, f64>: SubMapOps<First = DDS::StateVector, Second = RDS::StateVector>
{
    fn residue(&self, t: f64, y: DVector<f64>, x: EnumMap<StateVectorEnum, f64>, u: DVector<f64>) -> DVector<f64> {
        y.clone() - self.dirsys.y(
            t, 
            x.first(), 
            (u - self.revsys.y(
                t, 
                x.second(), 
                y.as_map()).as_vector()).as_map()
            ).as_vector()
    }

    pub fn new(dirsys: DDS, revsys: RDS) -> NegativeFeedback<DDS, RDS, StateVectorEnum> {
        // NegativeFeedback::<DDS, RDS>::assert_sizes();
        NegativeFeedback { dirsys, revsys, _phantom: PhantomData }
    }

    pub fn dir_ref(&self) -> &DDS {
        &self.dirsys
    }

    pub fn rev_ref(&self) -> &RDS {
        &self.revsys
    }
}

pub type UnityFeedback<DDS: DynamicalSystem, StateVectorEnum> = NegativeFeedback<DDS, UnitySystem<DDS::Input>, StateVectorEnum>;

impl<DDS: DynamicalSystem, StateVectorEnum> NegativeFeedback<DDS, UnitySystem<DDS::Input>, StateVectorEnum> 
where
    StateVectorEnum: EnumArray<f64>,
    EnumMap<StateVectorEnum, f64>: SubMapOps<First = DDS::StateVector, Second = Empty>,
    DVector<f64>: AsMap<DDS::Input>
{
    pub fn new_unity_feedback(dirsys: DDS) -> NegativeFeedback<DDS, UnitySystem<DDS::Input>, StateVectorEnum>
    {
        NegativeFeedback::new(dirsys, UnitySystem{_phantom: PhantomData})
    }
}

impl<DDS, RDS, StateVectorEnum: EnumArray<f64>> DynamicalSystem for 
    NegativeFeedback<DDS, RDS, StateVectorEnum> 
where
    DDS: DynamicalSystem,
    RDS: DynamicalSystem,
    StateVectorEnum: EnumArray<f64>,
    DVector<f64>: AsMap<RDS::Input> + AsMap<DDS::Input> + AsMap<DDS::Output>,
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
            |output| self.residue(t, output, x, u.as_vector())
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

    #[derive(Enum, Clone, Copy, Debug)]
    enum LFInp {
        X
    }

    impl AsMap<LFInp> for DVector<f64> {
        fn as_map(&self) -> EnumMap<LFInp, f64> {
            enum_map! {
                LFInp::X => self[0]
            }
        }
    }

    //terrível
    #[derive(Enum, Clone, Copy, Debug)]
    struct LFSV;

    #[derive(Enum, Clone, Copy, Debug)]
    enum LFOut {
        Y
    }

    impl AsMap<LFOut> for DVector<f64> {
        fn as_map(&self) -> EnumMap<LFOut, f64> {
            enum_map! {
                LFOut::Y => self[0]
            }
        }
    }

    impl ConvertMap<LFInp> for EnumMap<LFOut, f64> {
        fn convert(&self) -> EnumMap<LFInp, f64> {
            enum_map! {
                LFInp::X => self[LFOut::Y]
            }
        }
    }

    impl ConvertMap<LFOut> for EnumMap<LFInp, f64> {
        fn convert(&self) -> EnumMap<LFOut, f64> {
            enum_map! {
                LFOut::Y => self[LFInp::X]
            }
        }
    }

    impl DynamicalSystem for LinearFunc {
        type Input = LFInp;

        type StateVector = LFSV;

        type Output = LFOut;

        fn xdot(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> {
                x
        }

        fn y(&self, t: f64, 
            x: EnumMap<Self::StateVector, f64>, 
            u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> {
            enum_map! {
                LFOut::Y => self.a * u[LFInp::X]
            }
        }
    }

    StateVectorTypes!(Feedback = LFSV, LFSV);

    impl SubMapOps for EnumMap<FeedbackLFSVLFSV, f64> {
        type First = LFSV;

        type Second = LFSV;

        fn first(&self) -> EnumMap<Self::First, f64> {
            enum_map! {
                LFSV => 0.0
            }
        }

        fn second(&self) -> EnumMap<Self::Second, f64> {
            enum_map! {
                LFSV => 0.0
            }
        }

        fn merge(map1: EnumMap<Self::First, f64>, map2: EnumMap<Self::Second, f64>) -> Self {
            enum_map! {
                FeedbackLFSVLFSV::First(LFSV) => map1[LFSV],
                FeedbackLFSVLFSV::Second(LFSV) => map2[LFSV]
            }
        }
    }

    struct ExpSys {
        alpha: f64
    }

    
    #[derive(Enum, Clone, Copy, Debug)]
    enum ESInp {
        Reference
    }
    
    impl AsMap<ESInp> for DVector<f64> {
        fn as_map(&self) -> EnumMap<ESInp, f64> {
            enum_map! {
                ESInp::Reference => self[0]
            }
        }
    }

    impl ConvertMap<ESOut> for EnumMap<ESInp, f64> {
        fn convert(&self) -> EnumMap<ESOut, f64> {
            enum_map! {
                ESOut::Position => self[ESInp::Reference]
            }
        }
    }
    
    #[derive(Enum, Clone, Copy, Debug)]
    enum ESSV {
        Position
    }

    #[derive(Enum, Clone, Copy, Debug)]
    enum ESOut {
        Position
    }

    impl AsMap<ESOut> for DVector<f64> {
        fn as_map(&self) -> EnumMap<ESOut, f64> {
            enum_map! {
                ESOut::Position => self[0]
            }
        }
    }

    impl ConvertMap<ESInp> for EnumMap<ESOut, f64> {
        fn convert(&self) -> EnumMap<ESInp, f64> {
            enum_map! {
                ESInp::Reference => self[ESOut::Position]
            }
        }
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
                ESOut::Position => x[ESSV::Position]
            }
        }
    }

    StateVectorTypes!(Feedback = ESSV, ESSV);

    impl SubMapOps for EnumMap<FeedbackESSVESSV, f64> {
        type First = ESSV;

        type Second = ESSV;

        fn first(&self) -> EnumMap<Self::First, f64> {
            enum_map! {
                ESSV::Position => self[FeedbackESSVESSV::First(ESSV::Position)]
            }
        }

        fn second(&self) -> EnumMap<Self::Second, f64> {
            enum_map! {
                ESSV::Position => self[FeedbackESSVESSV::Second(ESSV::Position)]
            }
        }

        fn merge(map1: EnumMap<Self::First, f64>, map2: EnumMap<Self::Second, f64>) -> Self {
            enum_map! {
                FeedbackESSVESSV::First(ESSV::Position) => map1[ESSV::Position],
                FeedbackESSVESSV::Second(ESSV::Position) => map2[ESSV::Position]
            }
        }
    }

    #[test]
    fn test_feedback_output() {
        let sys = LinearFunc{a: 2.0};
        let feedback_sys: NegativeFeedback<LinearFunc, LinearFunc, FeedbackLFSVLFSV> = NegativeFeedback::new(sys, sys);

        let state = enum_map! {
            FeedbackLFSVLFSV::First(LFSV) => 0.0,
            FeedbackLFSVLFSV::Second(LFSV) => 0.0
        };

        let input = enum_map! {
            LFInp::X => 1.0
        };

        let output = feedback_sys.y(0.0, state, input);

        assert!((output[LFOut::Y] - 0.4).abs() < 1e-15);
    }

    #[test]
    fn test_feedback_xdot() {
        let exp1 = ExpSys{ alpha: 1.0 };
        let exp2 = ExpSys{ alpha: 1.0 };

        let feedback_sys: NegativeFeedback<ExpSys, ExpSys, FeedbackESSVESSV> = NegativeFeedback::new(
            exp1, exp2);

        let state = enum_map! {
            FeedbackESSVESSV::First(ESSV::Position) => 1.0,
            FeedbackESSVESSV::Second(ESSV::Position) => 2.0
        };

        let input = enum_map! {
            ESInp::Reference => 3.0
        };

        let xdot = feedback_sys.xdot(0.0, state, input);

        assert!(xdot.as_vector() == dvector![0.0, -1.0]);
    }
}