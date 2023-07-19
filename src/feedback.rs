use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;
use nalgebra::Dyn;

//direct dynamical system, reverse dynamical system
struct NegativeFeedback<DDS, RDS, const DSVS: usize, const RSVS: usize, const IS: usize, const OS: usize> {
    dirsys: DDS,
    revsys: RDS
}

impl<const DSVS: usize, const RSVS: usize, const IS: usize, const OS: usize, DDS, RDS> 
    NegativeFeedback<DDS, RDS, DSVS, RSVS, IS, OS> 
where
    DDS: DynamicalSystem<DSVS, IS, OS>,
    RDS: DynamicalSystem<RSVS, OS, IS> 
{
    fn residue(&self, t: f64, y: SVector<f64, OS>, x: SVector<f64, {DSVS+RSVS}>, u: SVector<f64, IS>) -> SVector<f64, OS> {
        y - self.dirsys.y(
            t, 
            SVector::from_row_slice(x.fixed_rows::<DSVS>(0).data.into_slice()), 
            u - self.revsys.y(
                t, 
                SVector::from_row_slice(x.fixed_rows::<RSVS>(DSVS).data.into_slice()), 
                y)
            )
    }
}



impl<const DSVS: usize, const RSVS: usize, const IS: usize, const OS: usize, DDS, RDS> 
    DynamicalSystem<{DSVS+RSVS}, IS, OS> for NegativeFeedback<DDS, RDS, DSVS, RSVS, IS, OS> 
where
    DDS: DynamicalSystem<DSVS, IS, OS>,
    RDS: DynamicalSystem<RSVS, OS, IS>    
{
    fn xdot(&self, t: f64, x: SVector<f64, {DSVS+RSVS}>, u: SVector<f64, IS>) -> SVector<f64, {DSVS+RSVS}> {
        todo!()
    }

    fn y(&self, t: f64, x: SVector<f64, {DSVS + RSVS}>, u: SVector<f64, IS>) -> SVector<f64, OS> {
        MultiVarNewtonFD::new(
            |output| self.residue(t, output, x, u)
        )
        .solve(self.dirsys.y(
            t, 
            SVector::from_row_slice(x.fixed_rows::<DSVS>(0).data.into_slice()), u))
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::ComplexField;

    use crate::DynamicalSystem;
    use super::*;

    #[derive(Clone, Copy)]
    struct LinearFunc {
        a: f64
    }

    impl DynamicalSystem<0, 1, 1> for LinearFunc {
        fn xdot(&self, t: f64, 
            x: nalgebra::SVector<f64, 0>, 
            u: nalgebra::SVector<f64, 1>) -> nalgebra::SVector<f64, 0> {
            todo!()
        }

        fn y(&self, t: f64, 
            x: nalgebra::SVector<f64, 0>, 
            u: nalgebra::SVector<f64, 1>) -> nalgebra::SVector<f64, 1> {
            self.a * u
        }
    }

    #[test]
    fn test_feedback_output() {
        let sys = LinearFunc{a: 2.0};
        let feedback_sys: NegativeFeedback<LinearFunc, LinearFunc, 0, 0, 1, 1> = NegativeFeedback{ 
            dirsys: sys, revsys: sys };

        let output = feedback_sys.y(0.0, [].into(), [1.0].into());

        assert!((output.x - 0.4).abs() < 1e-15);
    }
}