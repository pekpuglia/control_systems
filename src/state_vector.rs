use std::{ops::Sub, process::Output};

use nalgebra::SVector;

use crate::*;

//fazer conversão para DVector
//igualdade!!!
//indexação!!!
pub trait ComposableVector: Sub<Output=Self> + Sized + Copy {
    fn concat<CVO: ComposableVector>(self, other: CVO) -> VecConcat<Self, CVO> where Self: Sized {
        VecConcat(self, other)
    }

    fn into_dvector(&self) -> DVector<f64>;
}


impl<const N: usize> ComposableVector for SVector<f64, N> {
    fn into_dvector(&self) -> DVector<f64> {
        DVector::from_row_slice(self.data.as_slice())
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