use crate::*;

use self::state_vector::VecConcat;

#[derive(Clone, Copy, Debug)]
pub struct Parallel<TDS, BDS> {
    top_sys: TDS,
    bot_sys: BDS
}

impl<TDS, BDS> Parallel<TDS, BDS> {
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

impl<TIN, BIN, TST, BST, TOUT, BOUT, TDS, BDS> DynamicalSystem<VecConcat<TIN, BIN>, VecConcat<TST, BST>, VecConcat<TOUT, BOUT>> for Parallel<TDS, BDS> 
where
    TIN: ComposableVector,
    BIN: ComposableVector,
    TST: ComposableVector,
    BST: ComposableVector,
    TOUT: ComposableVector,
    BOUT: ComposableVector,
    TDS: DynamicalSystem<TIN, TST, TOUT>,
    BDS: DynamicalSystem<BIN, BST, BOUT>
{
    fn xdot(&self, t: f64, 
        x: &VecConcat<TST, BST>, 
        u: &VecConcat<TIN, BIN>) -> VecConcat<TST, BST> {
        self.top_sys.xdot(t, x.first(), u.first()).concat(
            self.bot_sys.xdot(t, x.second(), u.second())
        )
    }
    
    fn y(&self, t: f64, 
        x: &VecConcat<TST, BST>, 
        u: &VecConcat<TIN, BIN>) -> VecConcat<TOUT, BOUT> {
            self.top_sys.y(t, x.first(), u.first()).concat(
                self.bot_sys.y(t, x.second(), u.second())
            )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{dvector, vector, DVector};
    use crate::systems::Exp;

    const exp05: Exp = Exp::new(0.5);
    const exp1: Exp = Exp::new(1.0);


    #[test]
    fn test_parallel_xdot() {
        let par = Parallel {
            top_sys: exp05,
            bot_sys: exp1
        };

        let xdot = par.xdot(0.0, &vector![1.0].concat(vector![2.0]), &vector![1.0].concat(vector![0.0]));
        assert!(xdot.into_dvector() == dvector![0.0, -2.0])
    }

    #[test]
    fn test_parallel_y() {
        let par = Parallel {
            top_sys: exp05,
            bot_sys: exp1
        };

        let y = par.y(0.0, &vector![1.0].concat(vector![2.0]), &vector![1.0].concat(vector![0.0]));
        assert!(y.into_dvector() == dvector![1.0, 2.0])
    }
}