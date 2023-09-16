use core::panic;

use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;
use nalgebra::dvector;

#[derive(Clone, Copy, Debug)]
pub struct UnitySystem<const SIZE: usize>;

impl<const SIZE: usize> DynamicalSystem for UnitySystem<SIZE> {
    fn xdot(&self, _t: f64, 
        _x: DVector<f64>, 
        _u: DVector<f64>) -> DVector<f64> {
        dvector![]
    }

    fn y(&self, _t: f64, 
        _x: DVector<f64>, 
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
    fn residue(&self, t: f64, y: DVector<f64>, x: DVector<f64>, u: DVector<f64>) -> DVector<f64> {
        y.clone() - self.dirsys.y(
            t, 
            StateVector::<Self>::new(x.clone()).dirx().data, 
            u - self.revsys.y(
                t, 
                StateVector::<Self>::new(x.clone()).revx().data, 
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
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64> {
        
        let output = self.y(t, x.clone(), u.clone());
        let rev_output = self.revsys.y(t, x.rows(DDS::STATE_VECTOR_SIZE, RDS::STATE_VECTOR_SIZE).into(), output.clone());

        let mut dirxdot = self.dirsys.xdot(
            t, 
            StateVector::<Self>::new(x.clone()).dirx().data, 
            u - rev_output
        );

        let revxdot = self.revsys.xdot(
            t, 
            StateVector::<Self>::new(x.clone()).revx().data, 
            output
        );

        dirxdot
            .extend(revxdot.iter().copied());
        dirxdot
    }

    fn y(&self, t: f64, 
        x: DVector<f64>, 
        u: DVector<f64>) -> DVector<f64> {
        MultiVarNewtonFD::new(
            |output| self.residue(t, output, x.clone(), u.clone())
        )
        .solve(self.dirsys.y(
            t, 
            StateVector::<Self>::new(x.clone()).dirx().data, 
            u.clone())
        )
        .unwrap()
    }
}

impl<DDS: DynamicalSystem, RDS: DynamicalSystem> StateVector<NegativeFeedback<DDS, RDS>>  {
    pub fn dirx(&self) -> StateVector<DDS> {
        StateVector { 
            data: self.data.rows(0, DDS::STATE_VECTOR_SIZE).into(), 
            _phantom: PhantomData }
    }
    pub fn revx(&self) -> StateVector<RDS> {
        StateVector { 
            data: self.data.rows(DDS::STATE_VECTOR_SIZE, RDS::STATE_VECTOR_SIZE).into(), 
            _phantom: PhantomData }
    }
}


#[cfg(test)]
mod tests {

    use crate::DynamicalSystem;
    use super::*;

    #[derive(Clone, Copy)]
    struct LinearFunc {
        a: f64
    }

    impl DynamicalSystem for LinearFunc {
        fn xdot(&self, _t: f64, 
            _x: nalgebra::DVector<f64>, 
            _u: nalgebra::DVector<f64>) -> nalgebra::DVector<f64> {
            dvector![]
        }

        fn y(&self, _t: f64, 
            _x: nalgebra::DVector<f64>, 
            u: nalgebra::DVector<f64>) -> nalgebra::DVector<f64> {
            self.a * u
        }

        const STATE_VECTOR_SIZE: usize = 0;

        const INPUT_SIZE      : usize = 1;

        const OUTPUT_SIZE     : usize = 1;
    }

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
    fn test_feedback_output() {
        let sys = LinearFunc{a: 2.0};
        let feedback_sys = NegativeFeedback::new(sys, sys);

        let output = feedback_sys.y(0.0, dvector![], dvector![1.0]);

        assert!((output[0] - 0.4).abs() < 1e-15);
    }

    #[test]
    fn test_feedback_xdot() {
        let exp1 = Exp{ alpha: 1.0 };
        let exp2 = Exp{ alpha: 1.0 };

        let feedback_sys = NegativeFeedback { 
            dirsys: exp1, revsys: exp2 };

        let xdot = feedback_sys.xdot(0.0, dvector![1.0, 2.0], dvector![3.0]);

        assert!(xdot == dvector![0.0, -1.0]);
    }
}