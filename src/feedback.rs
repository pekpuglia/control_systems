use core::panic;

use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;
use nalgebra::dvector;

#[derive(Clone, Copy, Debug)]
pub struct UnitySystem<const SIZE: usize>;

impl<const SIZE: usize> DynamicalSystem for UnitySystem<SIZE> {
    fn xdot(&self, _t: f64, 
        _x: &StateVector<Self>, 
        _u: DVector<f64>) -> DVector<f64> {
        dvector![]
    }

    fn y(&self, _t: f64, 
        _x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64> {
        u
    }

    const STATE_VECTOR_SIZE: usize = 0;

    const INPUT_SIZE      : usize = SIZE;

    const OUTPUT_SIZE     : usize = SIZE;
}

#[derive(Debug, Clone, Copy)]
//direct dynamical system, reverse dynamical system
pub struct NegativeFeedback<DDS, RDS> 
{
    dirsys: DDS,
    revsys: RDS
}

impl<DDS, RDS> 
    NegativeFeedback<DDS, RDS>
where
    DDS: DynamicalSystem,
    RDS: DynamicalSystem 
{
    fn residue(&self, t: f64, y: DVector<f64>, x: &StateVector<Self>, u: DVector<f64>) -> DVector<f64> {
        y.clone() - self.dirsys.y(
            t, 
            &x.dirx(),
            u - self.revsys.y(
                t, 
                &x.revx(), 
                y)
            )
    }

    pub fn new(dirsys: DDS, revsys: RDS) -> NegativeFeedback<DDS, RDS> {
        NegativeFeedback::<DDS, RDS>::assert_sizes();
        NegativeFeedback { dirsys, revsys }
    }

    const fn assert_sizes() {
        if DDS::OUTPUT_SIZE != RDS::INPUT_SIZE || DDS::INPUT_SIZE != RDS::OUTPUT_SIZE {
            panic!("Wrong sizes!")
        }
    }

    pub fn dir_ref(&self) -> &DDS {
        &self.dirsys
    }

    pub fn rev_ref(&self) -> &RDS {
        &self.revsys
    }

    pub fn revy(&self, t: f64, x: StateVector<Self>, u: DVector<f64>) -> DVector<f64> {
        let output = self.y(t, &x, u.clone());
        self.revsys.y(t, &x.revx(), output.clone())
    }

    pub fn error(&self, t: f64, x: StateVector<Self>, u: DVector<f64>) -> DVector<f64> {
        u.clone() - self.revy(t, x, u)
    }
}

pub type UnityFeedback<DDS, const SIZE: usize> = NegativeFeedback<DDS, UnitySystem<SIZE>>;

impl<DDS: DynamicalSystem, const SIZE: usize> NegativeFeedback<DDS, UnitySystem<SIZE>> {
    pub fn new_unity_feedback(dirsys: DDS) -> NegativeFeedback<DDS, UnitySystem<SIZE>>
    {
        NegativeFeedback::new(dirsys, UnitySystem{})
    }
}

impl<DDS, RDS> DynamicalSystem for NegativeFeedback<DDS, RDS> 
where
    DDS: DynamicalSystem,
    RDS: DynamicalSystem,
    
{
    const STATE_VECTOR_SIZE: usize = DDS::STATE_VECTOR_SIZE + RDS::STATE_VECTOR_SIZE;

    const INPUT_SIZE      : usize = DDS::INPUT_SIZE;

    const OUTPUT_SIZE     : usize = DDS::OUTPUT_SIZE;

    fn xdot(&self, t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64> {
        
        let output = self.y(t, x.clone(), u.clone());
        let rev_output = self.revsys.y(t, &x.revx(), output.clone());

        let mut dirxdot = self.dirsys.xdot(
            t, 
            &x.dirx(), 
            u - rev_output
        );

        let revxdot = self.revsys.xdot(
            t, 
            &x.revx(),
            output
        );

        dirxdot
            .extend(revxdot.iter().copied());
        dirxdot
    }

    fn y(&self, t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64> {
        MultiVarNewtonFD::new(
            |output| self.residue(t, output, x, u.clone())
        )
        .solve(self.dirsys.y(
            t, 
            &x.dirx(), 
            u.clone())
        )
        .unwrap()
    }
}

impl<DDS: DynamicalSystem, RDS: DynamicalSystem> StateVector<NegativeFeedback<DDS, RDS>>  {
    pub fn dirx(&self) -> StateVector<DDS> {
        self.data
        .rows(
            0, 
            DDS::STATE_VECTOR_SIZE)
        .as_slice()
        .into_sv::<DDS>()
    }
    
    pub fn revx(&self) -> StateVector<RDS> {
        self.data
        .rows(
            DDS::STATE_VECTOR_SIZE, 
            RDS::STATE_VECTOR_SIZE)
        .as_slice()
        .into_sv::<RDS>() 
    }
}

impl<DDS: DynamicalSystem> StateVector<DDS>  {
    pub fn feedback<RDS: DynamicalSystem>(&self, rev_sv: StateVector<RDS>) -> StateVector<NegativeFeedback<DDS, RDS>> {
        let dataiter = self.data
        .iter()
        .chain(
            rev_sv.data
                .iter()
        ).copied();
        StateVector::new(DVector::from_iterator(
                DDS::STATE_VECTOR_SIZE + RDS::STATE_VECTOR_SIZE, 
                dataiter))
    }
}

#[cfg(test)]
mod tests {

    use crate::{DynamicalSystem, systems::{Exp, LinearFunc}};
    use super::*;

    const LF: LinearFunc = LinearFunc::new(2.0);
    const EXP1: Exp = Exp::new(1.0);

    #[test]
    fn test_feedback_output() {
        let feedback_sys = NegativeFeedback::new(LF, LF);

        let output = feedback_sys.y(0.0, &[].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]);

        assert!((output[0] - 0.4).abs() < 1e-15);
    }

    #[test]
    fn test_feedback_xdot() {

        let feedback_sys = NegativeFeedback::new(EXP1, EXP1);

        let xdot = feedback_sys.xdot(0.0, &[1.0, 2.0].into_sv::<NegativeFeedback<Exp, Exp>>(), dvector![3.0]);

        assert!(xdot == dvector![0.0, -1.0]);
    }

    #[test]
    fn test_ss_slicing() {
        let sv = StateVector::<NegativeFeedback<Exp, Exp>>::new(dvector![1.0, 2.0]);
        assert!(sv.dirx().data == [1.0].into());
        assert!(sv.revx().data == [2.0].into())
    }

    #[test]
    fn test_revy() {
        let feedback_sys = NegativeFeedback::new(LF, LF);

        assert!(feedback_sys.revy(0.0, [].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]) == dvector![0.8])
    }

    #[test]
    fn test_error() {
        let feedback_sys = NegativeFeedback::new(LF, LF);
        dbg!(feedback_sys.error(0.0, [].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]));
        assert!((feedback_sys.error(0.0,  [].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]) - dvector![0.2]).abs().max() < 1e-4);        
    }

    #[test]
    fn test_feedback_builder() {
        let sv_feedback = [1.0].into_sv::<Exp>()
            .feedback([2.0].into_sv::<Exp>());

        assert!(sv_feedback.data == dvector![1.0, 2.0])
    }
}