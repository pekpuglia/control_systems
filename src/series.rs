
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

macro_rules! StateVectorTypes {
    ($base_name: ident = $sv1:ty , $sv2:ty) => {
        paste!{
            #[derive(Enum, Clone, Copy, Debug)]
            enum [<$base_name $sv1 $sv2>] {
                First($sv1),
                Second($sv2)
            }
        }
    };
}

trait SubMapOps {
    type First: EnumArray<f64>;
    type Second: EnumArray<f64>;
    fn first(&self) -> EnumMap<Self::First, f64>;
    fn second(&self) -> EnumMap<Self::Second, f64>;
    fn merge(map1: EnumMap<Self::First, f64>, map2: EnumMap<Self::Second, f64>) -> Self;
}

trait ConvertMap<DestType: EnumArray<f64>> {
    fn convert(&self) -> EnumMap<DestType, f64>;
}


impl<DS1, DS2, StateVectorEnum: EnumArray<f64>> DynamicalSystem for Series<DS1, DS2, StateVectorEnum> 
where
    DS1: DynamicalSystem,
    DS2: DynamicalSystem,
    EnumMap<DS1::Output, f64>: ConvertMap<DS2::Input>,
    EnumMap<StateVectorEnum, f64>: SubMapOps<First = DS1::StateVector, Second = DS2::StateVector>,
    EnumMap<StateVectorEnum, f64>: Copy,
    EnumMap<DS1::Input, f64>: Copy
{
    type Input = DS1::Input;

    type StateVector = StateVectorEnum;

    type Output = DS2::Output;

    fn xdot(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64> 
    {
        let x1dot = self.dynsys1.xdot(t, x.first(), u);

        let x2dot = self.dynsys2.xdot(
            t, 
            x.second(), 
            self.dynsys1.y(t, x.first(), u).convert());

        SubMapOps::merge(x1dot, x2dot)
    }

    fn y(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64> 
    {
        self.dynsys2.y(t, x.second(), self.dynsys1.y(t, x.first(), u).convert())    
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DynamicalSystem;
    use nalgebra::{DVector, dvector};
    use enum_map::enum_map;
    struct ExpSys {
        alpha: f64
    }

    #[derive(Enum, Clone, Copy, Debug)]
    enum Input {
        Reference
    }

    #[derive(Enum, Clone, Copy, Debug)]
    enum State {
        Position
    }

    #[derive(Enum, Clone, Copy, Debug)]
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

    #[derive(Enum, Clone, Copy, Debug)]
    enum SOInp {
        Force
    }

    impl ConvertMap<SOInp> for EnumMap<Output, f64> {
        fn convert(&self) -> EnumMap<SOInp, f64> {
            enum_map! {
                SOInp::Force => self[Output::Position]
            }
        }
    }

    #[derive(Enum, Clone, Copy, Debug)]
    enum SOSV {
        Position,
        Velocity
    }

    #[derive(Enum, Clone, Copy, Debug)]
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

    }

    StateVectorTypes!(Series = State, SOSV);

    impl SubMapOps for EnumMap<SeriesStateSOSV, f64> {
        type First = State;

        type Second = SOSV;

        fn first(&self) -> EnumMap<Self::First, f64> {
            enum_map! {
                State::Position => self[SeriesStateSOSV::First(State::Position)]
            }
        }

        fn second(&self) -> EnumMap<Self::Second, f64> {
            enum_map! {
                SOSV::Position => self[SeriesStateSOSV::Second(SOSV::Position)],
                SOSV::Velocity => self[SeriesStateSOSV::Second(SOSV::Velocity)]
            }
        }

        fn merge(map1: EnumMap<Self::First, f64>, map2: EnumMap<Self::Second, f64>) -> Self {
            enum_map! {
                SeriesStateSOSV::First(State::Position) => map1[State::Position],
                SeriesStateSOSV::Second(SOSV::Position) => map2[SOSV::Position],
                SeriesStateSOSV::Second(SOSV::Velocity) => map2[SOSV::Velocity]
            }
        }


    }

    #[test]
    fn test_series_xdot() {
        let s: Series<ExpSys, SecondOrder, SeriesStateSOSV> = Series::new(ExpSys{alpha: 0.5}, SecondOrder{ k: 1.0, c: 1.0 });
        
        let state = enum_map! {
            SeriesStateSOSV::First(State::Position) => 0.0,
            SeriesStateSOSV::Second(SOSV::Position) => 0.0,
            SeriesStateSOSV::Second(SOSV::Velocity) => 1.0
        };

        let input = enum_map! {Input::Reference => 1.0};

        let xdot = s.xdot(0.0, state, input);
        dbg!(&xdot);
        assert!(xdot == enum_map! {
            SeriesStateSOSV::First(State::Position) => 0.5,
            SeriesStateSOSV::Second(SOSV::Position) => 1.0,
            SeriesStateSOSV::Second(SOSV::Velocity) => -1.0
        })

    }

    #[test]
    fn test_series_output() {
        // let series = Series::new(ExpSys{alpha: 0.5}, SecondOrder{ k: 1.0, c: 1.0 });

        // let state = enum_map! {
        //     SeriesStateSOSV::First(State::Position) => 0.5,
        //     SeriesStateSOSV::Second(SOSV::Position) => 0.7,
        //     SeriesStateSOSV::Second(SOSV::Velocity) => 1.0
        // };

        // let input = enum_map! {Input::Reference => 1.0};

        // let out = series.y(0.0, state, input);
        // assert!(out[SOOut::Position] == 0.7)
    }
}