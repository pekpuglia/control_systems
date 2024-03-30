use core::panic;
use std::{marker::PhantomData, process::Output};

use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;
use nalgebra::{dvector, vector, SVector};

use self::state_vector::VecConcat;

#[derive(Clone, Copy, Debug)]
pub struct UnitySystem;

impl<T: ComposableVector> DynamicalSystem<T, SVector<f64, 0>, T> for UnitySystem {
    fn xdot(&self, t: f64, 
        x: &SVector<f64, 0>, 
        u: &T) -> SVector<f64, 0> {
        vector![]
    }

    fn y(&self, t: f64, 
        x: &SVector<f64, 0>, 
        u: &T) -> T {
        *u
    }
}

#[derive(Debug, Clone, Copy)]
//direct dynamical system, reverse dynamical system
pub struct NegativeFeedback<DDS, RDS, DST, RST, DOUTRIN, DINROUT> 
where
    DST: ComposableVector,
    RST: ComposableVector,
    DOUTRIN: ComposableVector,
    DINROUT: ComposableVector,
    DDS: DynamicalSystem<DINROUT, DST, DOUTRIN>,
    RDS: DynamicalSystem<DOUTRIN, RST, DINROUT>
{
    dirsys: DDS,
    revsys: RDS,
    _phantom_dir_st: PhantomData<DST>,
    _phantom_rev_st: PhantomData<RST>,
    _phantom_dout_rin: PhantomData<DOUTRIN>,
    _phantom_din_rout: PhantomData<DINROUT>
}

impl<DDS, RDS, DST, RST, DOUTRIN, DINROUT> 
    NegativeFeedback<DDS, RDS, DST, RST, DOUTRIN, DINROUT> 
where
    DST: ComposableVector,
    RST: ComposableVector,
    DOUTRIN: ComposableVector,
    DINROUT: ComposableVector,
    DDS: DynamicalSystem<DINROUT, DST, DOUTRIN>,
    RDS: DynamicalSystem<DOUTRIN, RST, DINROUT>
{
    fn residue(&self, t: f64, y: &DOUTRIN, x: &VecConcat<DST, RST>, u: &DINROUT) -> DVector<f64> {
        (*y - self.dirsys.y(
            t, 
            &x.first(),
            &(*u - self.revsys.y(
                t, 
                &x.second(), 
                y))
            )).into_dvector()
    }

    pub fn new(dirsys: DDS, revsys: RDS) -> Self {
        NegativeFeedback { dirsys, revsys, _phantom_dir_st: PhantomData, _phantom_rev_st: PhantomData, _phantom_dout_rin: PhantomData, _phantom_din_rout: PhantomData }
    }

    pub fn dir_ref(&self) -> &DDS {
        &self.dirsys
    }

    pub fn rev_ref(&self) -> &RDS {
        &self.revsys
    }

    pub fn revy(&self, t: f64, x: &VecConcat<DST, RST>, u: &DINROUT) -> DINROUT {
        let output = self.y(t, x, u);
        self.revsys.y(t, x.second(), &output)
    }

    pub fn error(&self, t: f64, x: &VecConcat<DST, RST>, u: &DINROUT) -> DINROUT {
        *u - self.revy(t, x, u)
    }
}

pub type UnityFeedback<INOUT: ComposableVector, DST: ComposableVector, DDS> = NegativeFeedback<DDS, UnitySystem, DST, SVector<f64, 0>, INOUT, INOUT>;

impl<INOUT: ComposableVector, DST: ComposableVector, DDS: DynamicalSystem<INOUT, DST, INOUT>> UnityFeedback<INOUT, DST, DDS> {
    pub fn new_unity_feedback(dirsys: DDS) -> Self
    {
        NegativeFeedback::new(dirsys, UnitySystem{})
    }
}

impl<DDS, RDS, DST, RST, DOUTRIN, DINROUT>  DynamicalSystem<DINROUT, VecConcat<DST, RST>, DOUTRIN> for NegativeFeedback<DDS, RDS, DST, RST, DOUTRIN, DINROUT> 
where
    DST: ComposableVector,
    RST: ComposableVector,
    DOUTRIN: ComposableVector,
    DINROUT: ComposableVector,
    DDS: DynamicalSystem<DINROUT, DST, DOUTRIN>,
    RDS: DynamicalSystem<DOUTRIN, RST, DINROUT>
{
    fn xdot(&self, t: f64, 
        x: &VecConcat<DST, RST>, 
        u: &DINROUT) -> VecConcat<DST, RST> {
        let output = self.y(t, x, u);
        let rev_out = self.revsys.y(t, x.second(), &output);

        let dir_xdot = self.dirsys.xdot(t, x.first(), &(*u - rev_out));
        let rev_xdot = self.revsys.xdot(t, x.second(), &output);
        VecConcat::new(dir_xdot, rev_xdot)
    }

    fn y(&self, t: f64, 
        x: &VecConcat<DST, RST>, 
        u: &DINROUT) -> DOUTRIN {
        DOUTRIN::from_dvector(
            MultiVarNewtonFD::new(
            |output| self.residue(t, 
                                        &DOUTRIN::from_dvector(output).expect("output should have size equal to DOUTRIN"), x, u)
            )
            .solve(self.dirsys.y(
                t, 
                &x.first(), 
                u).into_dvector()
            ).expect("solver should find a solution")
        ).expect("result should have appropriate size")
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::vector;

    use super::*;

    use crate::systems::{Exp, LinearFunc};

    const LF: LinearFunc = LinearFunc::new(2.0);

    const EXP: Exp = Exp::new(1.0);

    #[test]
    fn test_feedback_output() {
        let feedback_sys = NegativeFeedback::new(LF, LF);

        let output = feedback_sys.y(0.0, &vector![].concat(vector![]), &vector![1.0]);

        assert!((output[0] - 0.4).abs() < 1e-15);
    }

    #[test]
    fn test_feedback_xdot() {

        let feedback_sys = NegativeFeedback::new(EXP, EXP);

        let xdot = feedback_sys.xdot(0.0, &vector![1.0].concat(vector![2.0]), &vector![3.0]);
        dbg!(xdot);
        assert!(xdot.into_dvector() == dvector![0.0, -1.0]);
    }

    #[test]
    fn test_revy() {
        let feedback_sys = NegativeFeedback::new(LF, LF);

        assert!(feedback_sys.revy(0.0, &vector![].concat(vector![]), &vector![1.0]).into_dvector() == dvector![0.8])
    }

    #[test]
    fn test_error() {
        let feedback_sys = NegativeFeedback::new(LF, LF);
        assert!((feedback_sys.error(0.0, &vector![].concat(vector![]), &vector![1.0]) - vector![0.2]).abs().max() < 1e-4);        
    }

    #[test]
    fn test_unity_feedback() {
        let uf = UnityFeedback::new_unity_feedback(LF);
        assert!(uf.y(0.0, &vector![].concat(vector![]), &vector![1.0]).x == (2.0/3.0))
    }
}