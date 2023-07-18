use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;

//direct dynamical system, reverse dynamical system
struct NegativeFeedback<DDS, RDS, const DSVS: usize, const RSVS: usize> {
    dirsys: DDS,
    revsys: RDS
}

impl<const DSVS: usize, const RSVS: usize, const IS: usize, const OS: usize, DDS, RDS> 
    DynamicalSystem<{DSVS+RSVS}, IS, OS> for NegativeFeedback<DDS, RDS, DSVS, RSVS> 
where
    DDS: DynamicalSystem<DSVS, IS, OS>,
    RDS: DynamicalSystem<RSVS, OS, IS>    
{
    fn xdot(&self, t: f64, x: SVector<f64, {DSVS+RSVS}>, u: SVector<f64, IS>) -> SVector<f64, {DSVS+RSVS}> {
        todo!()
    }

    //this is most likely only an approximation
    fn y(&self, t: f64, x: SVector<f64, {DSVS + RSVS}>, u: SVector<f64, IS>) -> SVector<f64, OS> {
        // MultiVarNewtonFD::new(
        //     |output| output
        // )
        // .solve(self.dirsys.y(t, x, u))
        // .unwrap()
        todo!()
    }
}