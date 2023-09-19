use crate::*;

#[derive(Clone, Copy, Debug)]
pub struct Parallel<TDS, BDS> {
    top_sys: TDS,
    bot_sys: BDS
}

impl<TDS: DynamicalSystem, BDS: DynamicalSystem> Parallel<TDS, BDS> {
    pub fn new(top_sys: TDS, bot_sys: BDS) -> Self {
        Self { top_sys, bot_sys }
    }

    pub fn top_sys_ref(&self) -> &TDS {
        &self.top_sys
    }

    pub fn bot_sys_ref(&self) -> &BDS {
        &self.bot_sys
    }

    //inspetores
}

impl<TDS: DynamicalSystem, BDS: DynamicalSystem> DynamicalSystem for Parallel<TDS, BDS> {
    const STATE_VECTOR_SIZE: usize = TDS::STATE_VECTOR_SIZE + BDS::STATE_VECTOR_SIZE;

    const INPUT_SIZE      : usize = TDS::INPUT_SIZE + BDS::INPUT_SIZE;

    const OUTPUT_SIZE     : usize = TDS::INPUT_SIZE + BDS::INPUT_SIZE;

    fn xdot(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64> {
        let mut txdot = self.top_sys.xdot(
            t, 
            StateVector::<Self>::new(x.clone()).topx().data,
            u.rows(0, TDS::INPUT_SIZE).into()
        );

        txdot
            .extend(
                self.bot_sys.xdot(
                    t, 
                    StateVector::<Self>::new(x).botx().data, 
                    u.rows(TDS::INPUT_SIZE, BDS::INPUT_SIZE).into()
                ).iter()
                .copied()
        );

        txdot
    }

    fn y(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64> {
        let mut ty = self.top_sys.y(
            t, 
            StateVector::<Self>::new(x.clone()).topx().data, 
            u.rows(0, TDS::INPUT_SIZE).into());

        ty
            .extend(
                self.bot_sys.y(
                    t, 
                    StateVector::<Self>::new(x).botx().data, 
                    u.rows(TDS::INPUT_SIZE, BDS::INPUT_SIZE).into()
                ).iter()
                .copied()
        );

        ty
    }
}

impl <TDS: DynamicalSystem, BDS: DynamicalSystem> StateVector<Parallel<TDS, BDS>> {
    pub fn topx(&self) -> StateVector<TDS> {
        self.data
        .rows(
            0, 
            TDS::STATE_VECTOR_SIZE)
        .as_slice()
        .into()
    }

    pub fn botx(&self) -> StateVector<BDS> {
        self.data
        .rows(
            TDS::STATE_VECTOR_SIZE, 
            BDS::STATE_VECTOR_SIZE)
        .as_slice()
        .into()
    }
}

impl<TDS: DynamicalSystem> StateVector<TDS>  {
    pub fn parallel<BDS: DynamicalSystem>(&self, bot_sv: StateVector<BDS>) -> StateVector<Parallel<TDS, BDS>> {
        let dataiter = self.data
        .iter()
        .chain(
            bot_sv.data
                .iter()
        ).copied();
        StateVector::new(
            DVector::from_iterator(
                TDS::STATE_VECTOR_SIZE + BDS::STATE_VECTOR_SIZE, 
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

    #[test]
    fn test_parallel_xdot() {
        let par = Parallel {
            top_sys: Exp{alpha: 0.5},
            bot_sys: Exp{alpha: 1.0}
        };

        let xdot = par.xdot(0.0, dvector![1.0, 2.0], dvector![1.0, 0.0]);
        assert!(xdot == dvector![0.0, -2.0])
    }

    #[test]
    fn test_parallel_y() {
        let par = Parallel {
            top_sys: Exp{alpha: 0.5},
            bot_sys: Exp{alpha: 1.0}
        };

        let y = par.y(0.0, dvector![1.0, 2.0], dvector![1.0, 0.0]);
        assert!(y == dvector![1.0, 2.0])
    }

    #[test]
    fn test_parallel_sv_creation() {
        let sv = [1.0].into_sv::<Exp>().parallel([2.0].into_sv::<Exp>());
        assert!(sv.data == dvector![1.0, 2.0])
    }
}