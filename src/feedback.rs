use crate::*;

use eqsolver::multivariable::MultiVarNewtonFD;

pub struct UnitySystem<const SIZE: usize>;

impl<const SIZE: usize> DynamicalSystem<0, SIZE, SIZE> for UnitySystem<SIZE> {
    fn xdot(&self, _t: f64, 
        _x: SVector<f64, 0>, 
        _u: SVector<f64, SIZE>) -> SVector<f64, 0> {
        [].into()
    }

    fn y(&self, _t: f64, 
        _x: SVector<f64, 0>, 
        u: SVector<f64, SIZE>) -> SVector<f64, SIZE> {
        u
    }
}

//direct dynamical system, reverse dynamical system
pub struct NegativeFeedback<DDS, RDS, const DSVS: usize, const RSVS: usize, const IS: usize, const OS: usize> 
where
    DDS: DynamicalSystem<DSVS, IS, OS>,
    RDS: DynamicalSystem<RSVS, OS, IS> 
{
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

    pub fn new(dirsys: DDS, revsys: RDS) -> NegativeFeedback<DDS, RDS, DSVS, RSVS, IS, OS> {
        NegativeFeedback { dirsys, revsys }
    }

    pub fn new_unity_feedback<SquareDS>(dirsys: SquareDS) -> NegativeFeedback<SquareDS, UnitySystem<IS>, DSVS, 0, IS, IS>
    where
        SquareDS: DynamicalSystem<DSVS, IS, IS>,
    {
        NegativeFeedback{ dirsys, 
            revsys: UnitySystem{} }
    }
}



impl<const DSVS: usize, const RSVS: usize, const IS: usize, const OS: usize, DDS, RDS> 
    DynamicalSystem<{DSVS+RSVS}, IS, OS> for NegativeFeedback<DDS, RDS, DSVS, RSVS, IS, OS> 
where
    DDS: DynamicalSystem<DSVS, IS, OS>,
    RDS: DynamicalSystem<RSVS, OS, IS>    
{
    fn xdot(&self, t: f64, x: SVector<f64, {DSVS+RSVS}>, u: SVector<f64, IS>) -> SVector<f64, {DSVS+RSVS}> {
        let output = self.y(t, x, u);

        let dirxdot = self.dirsys.xdot(
            t, 
            SVector::from_row_slice(x.fixed_rows::<DSVS>(0).data.into_slice()), 
            u - self.revsys.y(
                t, 
                SVector::from_row_slice(x.fixed_rows::<RSVS>(DSVS).data.into_slice()), 
                output
            )
        );

        let revxdot = self.revsys.xdot(
            t, 
            SVector::from_row_slice(x.fixed_rows::<RSVS>(DSVS).data.into_slice()), 
            output
        );

        SVector::from_iterator(
            dirxdot
            .into_iter()
            .chain(revxdot.into_iter())
            .copied()
        )
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

    use crate::DynamicalSystem;
    use super::*;

    #[derive(Clone, Copy)]
    struct LinearFunc {
        a: f64
    }

    impl DynamicalSystem<0, 1, 1> for LinearFunc {
        fn xdot(&self, _t: f64, 
            _x: nalgebra::SVector<f64, 0>, 
            _u: nalgebra::SVector<f64, 1>) -> nalgebra::SVector<f64, 0> {
            [].into()
        }

        fn y(&self, _t: f64, 
            _x: nalgebra::SVector<f64, 0>, 
            u: nalgebra::SVector<f64, 1>) -> nalgebra::SVector<f64, 1> {
            self.a * u
        }
    }

    struct Exp {
        alpha: f64
    }

    impl DynamicalSystem<1, 1, 1> for Exp {
        fn xdot(&self, _t: f64, x: SVector<f64, 1>, u: SVector<f64, 1>) -> SVector<f64, 1> {
            self.alpha * (u - x)
        }

        fn y(&self, _t: f64, x: SVector<f64, 1>, _u: SVector<f64, 1>) -> SVector<f64, 1> {
            x
        }
    }

    #[test]
    fn test_feedback_output() {
        let sys = LinearFunc{a: 2.0};
        let feedback_sys = NegativeFeedback{ 
            dirsys: sys, revsys: sys };

        let output = feedback_sys.y(0.0, [].into(), [1.0].into());

        assert!((output.x - 0.4).abs() < 1e-15);
    }

    #[test]
    fn test_feedback_xdot() {
        let exp1 = Exp{ alpha: 0.0 };
        let exp2 = Exp{ alpha: 1.0 };

        let feedback_sys = NegativeFeedback { 
            dirsys: exp1, revsys: exp2 };

        let xdot = feedback_sys.xdot(0.0, [1.0, 2.0].into(), [3.0].into());

        assert!(xdot == SVector::<f64, 2>::from_row_slice(&[0.0, -1.0]));
    }
}