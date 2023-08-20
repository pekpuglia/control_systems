use enum_map::EnumMap;
pub use nalgebra::{DVector, dvector};
pub use enum_map::{EnumArray, Enum};


//todo
//remove cloning?
//make it size-safe!
//system lib?

pub trait DynamicalSystem {

    type Input: EnumArray<f64>;
    type StateVector: EnumArray<f64>;
    type Output: EnumArray<f64>;

    fn xdot(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::StateVector, f64>;

    fn y(&self, t: f64, 
        x: EnumMap<Self::StateVector, f64>, 
        u: EnumMap<Self::Input, f64>) -> EnumMap<Self::Output, f64>;
}

#[cfg(test)]
mod test {
    use super::*;
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
        Error
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
                Output::Error => u[Input::Reference] - x[State::Position]
            }
        }



    }

    #[test]
    fn test_trait_impl() {
        let sys = ExpSys{ alpha: 0.1 };

        let state = enum_map! {State::Position => 0.0};

        let input = enum_map! {Input::Reference => 1.0};

        let deriv = sys.xdot(0.0, state, input);

        //garante que essa conversão é possível
        let x: [f64; 1] = deriv.into_array();

        assert!(x == [0.1])
    }
}

mod series;
pub use series::Series;


// mod feedback;
// pub use feedback::{NegativeFeedback, UnitySystem, UnityFeedback};