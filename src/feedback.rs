use core::panic;
use std::{marker::PhantomData, process::Output};

use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;
use nalgebra::dvector;

use self::state_vector::VecConcat;

// #[derive(Clone, Copy, Debug)]
// pub struct UnitySystem<const SIZE: usize>;

// impl<const SIZE: usize> DynamicalSystem for UnitySystem<SIZE> {
//     fn xdot(&self, _t: f64, 
//         _x: &StateVector<Self>, 
//         _u: DVector<f64>) -> DVector<f64> {
//         dvector![]
//     }

//     fn y(&self, _t: f64, 
//         _x: &StateVector<Self>, 
//         u: DVector<f64>) -> DVector<f64> {
//         u
//     }

//     const STATE_VECTOR_SIZE: usize = 0;

//     const INPUT_SIZE      : usize = SIZE;

//     const OUTPUT_SIZE     : usize = SIZE;
// }

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

// pub type UnityFeedback<DDS, const SIZE: usize> = NegativeFeedback<DDS, UnitySystem<SIZE>>;

// impl<DDS: DynamicalSystem, const SIZE: usize> NegativeFeedback<DDS, UnitySystem<SIZE>> {
//     pub fn new_unity_feedback(dirsys: DDS) -> NegativeFeedback<DDS, UnitySystem<SIZE>>
//     {
//         NegativeFeedback::new(dirsys, UnitySystem{})
//     }
// }

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

        let dir_xdot = self.dirsys.xdot(t, x.first(), u);
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

// #[cfg(test)]
// mod tests {

//     use crate::DynamicalSystem;
//     use super::*;

//     #[derive(Clone, Copy)]
//     struct LinearFunc {
//         a: f64
//     }

//     impl DynamicalSystem for LinearFunc {
//         fn xdot(&self, _t: f64, 
//             _x: &StateVector<Self>, 
//             _u: nalgebra::DVector<f64>) -> nalgebra::DVector<f64> {
//             dvector![]
//         }

//         fn y(&self, _t: f64, 
//             _x: &StateVector<Self>, 
//             u: nalgebra::DVector<f64>) -> nalgebra::DVector<f64> {
//             self.a * u
//         }

//         const STATE_VECTOR_SIZE: usize = 0;

//         const INPUT_SIZE      : usize = 1;

//         const OUTPUT_SIZE     : usize = 1;
//     }
//     #[derive(Clone, Copy)]
//     struct Exp {
//         alpha: f64
//     }

//     impl DynamicalSystem for Exp {
//         fn xdot(&self, _t: f64, x: &StateVector<Self>, u: DVector<f64>) -> DVector<f64> {
//             self.alpha * (u - x.data.clone())
//         }

//         fn y(&self, _t: f64, x: &StateVector<Self>, _u: DVector<f64>) -> DVector<f64> {
//             x.data.clone()
//         }

//         const STATE_VECTOR_SIZE: usize = 1;

//         const INPUT_SIZE      : usize = 1;

//         const OUTPUT_SIZE     : usize = 1;
//     }

//     #[test]
//     fn test_feedback_output() {
//         let sys = LinearFunc{a: 2.0};
//         let feedback_sys = NegativeFeedback::new(sys, sys);

//         let output = feedback_sys.y(0.0, &[].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]);

//         assert!((output[0] - 0.4).abs() < 1e-15);
//     }

//     #[test]
//     fn test_feedback_xdot() {
//         let exp1 = Exp{ alpha: 1.0 };
//         let exp2 = Exp{ alpha: 1.0 };

//         let feedback_sys = NegativeFeedback { 
//             dirsys: exp1, revsys: exp2 };

//         let xdot = feedback_sys.xdot(0.0, &[1.0, 2.0].into_sv::<NegativeFeedback<Exp, Exp>>(), dvector![3.0]);

//         assert!(xdot == dvector![0.0, -1.0]);
//     }

//     #[test]
//     fn test_ss_slicing() {
//         let sv = StateVector::<NegativeFeedback<Exp, Exp>>::new(dvector![1.0, 2.0]);
//         assert!(sv.dirx().data == [1.0].into());
//         assert!(sv.revx().data == [2.0].into())
//     }

//     #[test]
//     fn test_revy() {
//         let sys = LinearFunc{a: 2.0};
//         let feedback_sys = NegativeFeedback::new(sys, sys);

//         assert!(feedback_sys.revy(0.0, [].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]) == dvector![0.8])
//     }

//     #[test]
//     fn test_error() {
//         let sys = LinearFunc{a: 2.0};
//         let feedback_sys = NegativeFeedback::new(sys, sys);
//         dbg!(feedback_sys.error(0.0, [].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]));
//         assert!((feedback_sys.error(0.0,  [].into_sv::<NegativeFeedback<LinearFunc, LinearFunc>>(), dvector![1.0]) - dvector![0.2]).abs().max() < 1e-4);        
//     }

//     #[test]
//     fn test_feedback_builder() {
//         let sv_feedback = [1.0].into_sv::<Exp>()
//             .feedback([2.0].into_sv::<Exp>());

//         assert!(sv_feedback.data == dvector![1.0, 2.0])
//     }
// }