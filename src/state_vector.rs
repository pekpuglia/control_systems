use std::{fmt::format, ops::Sub, process::Output};

use nalgebra::{SVector, Vector};

use crate::*;

//fazer conversão para DVector
//igualdade!!!
//indexação!!!
pub trait ComposableVector: Sub<Output=Self> + Sized + Copy {
    fn concat<CVO: ComposableVector>(self, other: CVO) -> VecConcat<Self, CVO> where Self: Sized {
        VecConcat(self, other)
    }

    fn into_dvector(&self) -> DVector<f64>;

    fn from_dvector(v: DVector<f64>) -> Result<Self, String>;

    fn size() -> usize;
}


impl<const N: usize> ComposableVector for SVector<f64, N> {
    fn into_dvector(&self) -> DVector<f64> {
        DVector::from_row_slice(self.data.as_slice())
    }
    
    fn from_dvector(v: DVector<f64>) -> Result<Self, String> {
        let s = Self::size();
        match v.len() {
            s => Ok(SVector::from_iterator(v.iter().copied())),
            l => Err(format!("wrong size, expected {}, got {}", N, l))
        }
    }
    
    fn size() -> usize {
        N
    }
}

#[derive(Clone, Copy, Debug)]
pub struct VecConcat<CV0: ComposableVector, CV1: ComposableVector>(CV0, CV1);

impl<CV0: ComposableVector, CV1: ComposableVector> VecConcat<CV0, CV1> {
    pub fn first(&self) -> &CV0 {
        &self.0
    }

    pub fn second(&self) -> &CV1 {
        &self.1
    }

    pub fn new(v0: CV0, v1: CV1) -> VecConcat<CV0, CV1> {
        Self(v0, v1)
    }
}

impl<CV0: ComposableVector , CV1: ComposableVector> Sub for VecConcat<CV0, CV1> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        VecConcat::new(self.0 - rhs.0, self.1-rhs.1)
    }
}

impl<CV0: ComposableVector, CV1: ComposableVector> ComposableVector for VecConcat<CV0, CV1> {
    fn into_dvector(&self) -> DVector<f64> {
        let mut ret = self.0.into_dvector();
        ret.extend(self.1.into_dvector().iter().copied());
        ret
    }
    
    fn from_dvector(v: DVector<f64>) -> Result<Self, String> {
        let s = Self::size();
        match v.len() {
            s => {
                let dv0 = v.rows(0, CV0::size()).into();
                let dv1 = v.rows(CV0::size(), CV1::size()).into();
                Ok(Self(CV0::from_dvector(dv0).expect("dv0 should have the appropriate size"), CV1::from_dvector(dv1).expect("dv1 should have appropriate size")))
            },
            l => Err(format!("wrong size, expected {}, got {}", s, l))
        }
    }
    
    fn size() -> usize {
        CV0::size() + CV1::size()
    }
}

#[cfg(test)]
mod test {
    use nalgebra::{dvector, SVector};

    use crate::ComposableVector;

    use super::VecConcat;

    #[test]
    fn test_from_dvector() {
        let dv = dvector![1.0, 2.0, 3.0];
        let vc = VecConcat::<SVector<f64, 1>, SVector<f64, 2>>::from_dvector(dv);
        assert!(vc.is_ok());
        assert!(vc.as_ref().unwrap().0.x == 1.0);
        assert!(vc.as_ref().unwrap().1.x == 2.0);
        assert!(vc.as_ref().unwrap().1.y == 3.0)
    }
}

// #[derive(Clone, Debug)]
// pub struct StateVector<System: DynamicalSystem + ?Sized> {
//     pub data: DVector<f64>,
//     _phantom: PhantomData<System>
// }

// impl<System: DynamicalSystem> StateVector<System>  {
//     pub fn new(x: DVector<f64>) -> Self {
//         assert!(x.len() == System::STATE_VECTOR_SIZE);
//         StateVector { data: x, _phantom: PhantomData }
//     }
// }
// //refazer o from dvector
// pub trait IntoSV {
//     fn into_sv<System: DynamicalSystem>(self) -> StateVector<System>;
// }

// impl<const N: usize> IntoSV for [f64; N] {
//     fn into_sv<System: DynamicalSystem>(self) -> StateVector<System> {
//         assert!(N == System::STATE_VECTOR_SIZE);
//         StateVector { data: DVector::from_iterator(N, self.into_iter()), _phantom: PhantomData }
//     }
// }

// impl IntoSV for &[f64] {
//     fn into_sv<System: DynamicalSystem>(self) -> StateVector<System> {
//         assert!(self.len() == System::STATE_VECTOR_SIZE);
//         StateVector { data: DVector::from_row_slice(self), _phantom: PhantomData }
//     }
// }

// use std::ops::Index;

// impl<System: DynamicalSystem> Index<usize> for StateVector<System> {
//     type Output = f64;

//     fn index(&self, index: usize) -> &Self::Output {
//         &self.data[index]
//     }


// }