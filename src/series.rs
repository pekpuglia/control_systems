
use crate::*;
#[derive(Clone, Copy, Debug)]
pub struct Series<DS1, DS2>
{
    dynsys1: DS1,
    dynsys2: DS2
}

impl<DS1: DynamicalSystem, DS2: DynamicalSystem> Series<DS1, DS2>  {
    pub fn new(ds1: DS1, ds2: DS2) -> Self {
        Series { dynsys1: ds1, dynsys2: ds2 }
    }
    pub fn ds1_ref(&self) -> &DS1 {
        &self.dynsys1
    }

    pub fn ds2_ref(&self) -> &DS2 {
        &self.dynsys2
    }

    pub fn y1(&self, t: f64, x: StateVector<Self>, u:DVector<f64>) -> DVector<f64> {
        self.dynsys1.y(
            t, 
            &x.x1(), 
            u
        )
    }
}

impl<DS1, DS2> DynamicalSystem for Series<DS1, DS2> 
where
    DS1: DynamicalSystem,
    DS2: DynamicalSystem
{
    const STATE_VECTOR_SIZE: usize = DS1::STATE_VECTOR_SIZE + DS2::STATE_VECTOR_SIZE;

    const INPUT_SIZE      : usize = DS1::INPUT_SIZE;

    const OUTPUT_SIZE     : usize = DS2::OUTPUT_SIZE;

    fn xdot(&self, t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64> {
        let mut x1dot = self.dynsys1.xdot(
            t, 
            &x.x1(), 
            u.clone()
        );
        let x2dot = self.dynsys2.xdot(
            t, 
            &x.x2(), 
            self.dynsys1.y(
                t, &x.x1(), u)
        );

        x1dot
            .extend(x2dot.iter().copied());
        x1dot
    }

    fn y(&self, t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64> {
        let u1 = self.dynsys1.y(
            t, 
            &x.x1(), 
            u
        );

        self.dynsys2.y(
            t, 
            &x.x2(), 
            u1)
    }
}


impl<DS1: DynamicalSystem, DS2: DynamicalSystem> StateVector<Series<DS1, DS2>> {
    pub fn x1(&self) -> StateVector<DS1> {
        self.data
        .rows(
            0, 
            DS1::STATE_VECTOR_SIZE)
        .as_slice()
        .into_sv::<DS1>()
    }
    pub fn x2(&self) -> StateVector<DS2> {
        self.data
        .rows(
            DS1::STATE_VECTOR_SIZE, 
            DS2::STATE_VECTOR_SIZE)
        .as_slice()
        .into_sv::<DS2>()
    }
}

impl<System: DynamicalSystem> StateVector<System> {
    pub fn series<S2: DynamicalSystem>(&self, sv2: StateVector<S2>) -> StateVector<Series<System, S2>>{
        let dataiter = self.data
        .iter()
        .chain(
            sv2.data
                .iter()
        ).copied();
        StateVector::new(DVector::from_iterator(
                System::STATE_VECTOR_SIZE + S2::STATE_VECTOR_SIZE, 
                dataiter)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::systems::*;
    use nalgebra::{vector};

    const EXP: Exp = Exp::new(0.5);
    const SO: SecondOrder = SecondOrder::new(1.0, 1.0);

    #[test]
    fn test_series_xdot() {
        let series = Series::new(
            EXP, SO
        );

        let xdot = series.xdot(0.0, &[0.0,0.0, 1.0].into_sv::<Series<Exp, SecondOrder>>(), dvector![1.0]);
        dbg!(&xdot);
        assert!(xdot == DVector::from_row_slice(&[0.5, 1.0, -1.0]))
    }

    #[test]
    fn test_series_output() {
        let series = Series::new(EXP, SO);

        let out = series.y(0.0, &[0.5, 0.7, 1.0].into_sv::<Series<Exp, SecondOrder>>(), dvector![1.0]);
        assert!(out == [0.7].into())
    }

    #[test]
    fn test_state_vector_slicing() {
        let x = dvector![1.0, 2.0, 3.0];
        
        let sv = StateVector::<Series<Exp, SecondOrder>>::new(x);
        
        assert!(sv.x1().data == [1.0].into());
        assert!(sv.x2().data == dvector![2.0, 3.0]);
    }

    #[test]
    fn test_y1() {
        let series = Series::new(EXP, SO);

        let sys1 = series.ds1_ref();

        assert!(series.y1(0.0, [0.5, 0.7, 1.0].into_sv::<Series<Exp, SecondOrder>>(), dvector![1.0]) == sys1.y(0.0, &[0.5].into_sv::<Exp>(), dvector![1.0]))
    }

    #[test]
    fn test_series_builder() {
        let svseries = [0.5].into_sv::<Exp>()
            .series(
                [0.7, 1.0].into_sv::<SecondOrder>()
        );
        assert!(svseries.data == dvector![0.5, 0.7, 1.0])
    }
}
