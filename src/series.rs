
use std::marker::PhantomData;

use crate::*;

use self::state_vector::VecConcat;
#[derive(Clone, Copy, Debug)]
pub struct Series<DS1, DS2, OUT1IN2: ComposableVector>
where
    DS1: DynamicalSystem<OUT = OUT1IN2>,
    DS2: DynamicalSystem<IN = OUT1IN2>
{
    dynsys1: DS1,
    dynsys2: DS2
}

impl<DS1, OUT1IN2, DS2> Series<DS1, DS2, OUT1IN2>
where
    OUT1IN2: ComposableVector,
    DS1: DynamicalSystem<OUT = OUT1IN2>,
    DS2: DynamicalSystem<IN = OUT1IN2>
{
    pub fn new(ds1: DS1, ds2: DS2) -> Self {
        Series { dynsys1: ds1, dynsys2: ds2 }
    }
    pub fn ds1_ref(&self) -> &DS1 {
        &self.dynsys1
    }

    pub fn ds2_ref(&self) -> &DS2 {
        &self.dynsys2
    }

    pub fn y1(&self, t: f64, x: &VecConcat<DS1::ST, DS2::ST>, u: &DS1::IN) -> OUT1IN2 {
        self.dynsys1.y(
            t, 
            x.first(), 
            u
        )
    }
}

impl<DS1, OUT1IN2, DS2> DynamicalSystem for Series<DS1, DS2, OUT1IN2> 
where
    OUT1IN2: ComposableVector,
    DS1: DynamicalSystem<OUT = OUT1IN2>,
    DS2: DynamicalSystem<IN = OUT1IN2> 
{
    fn xdot(&self, t: f64, 
        x: &Self::ST, 
        u: &Self::IN) -> Self::ST {
        self.dynsys1.xdot(t, x.first(), u).concat(
            self.dynsys2.xdot(t, x.second(), 
                &self.dynsys1.y(
                    t, x.first(), u
                )
            ))
    }
    
    fn y(&self, t: f64, 
        x: &Self::ST, 
        u: &Self::IN) -> Self::OUT {
            self.dynsys2.y(t, x.second(), &self.dynsys1.y(t, x.first(), u))
    }
    
    type IN = DS1::IN;
    
    type ST = VecConcat<DS1::ST, DS2::ST>;
    
    type OUT = DS2::OUT;
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

        let xdot = series.xdot(0.0, &vector![0.0].concat(vector![0.0, 1.0]), &vector![1.0]);
        assert!(xdot.first().x == 0.5);
        assert!(xdot.second().x == 1.0);
        assert!(xdot.second().y == -1.0);
    }

    #[test]
    fn test_series_output() {
        let series = Series::new(EXP, SO);

        let out = series.y(0.0, &vector![0.5].concat(vector![0.7, 1.0]), &vector![1.0]);
        assert!(out == [0.7].into())
    }

    #[test]
    fn test_y1() {
        let series = Series::new(EXP, SO);

        let sys1 = series.ds1_ref();

        assert!(series.y1(0.0, &vector![0.5].concat(vector![0.7, 1.0]), &vector![1.0]) == sys1.y(0.0, &vector![0.5], &vector![1.0]))
    }

    //TODO BOTAR EM MACRO - testar
    type SESO = Series<Exp, SecondOrder, <Exp as DynamicalSystem>::IN>;
    #[test]
    fn test_type_alias() {
        let seso = SESO::new(EXP, SO);
        
    }
}
