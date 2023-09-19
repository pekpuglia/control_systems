use crate::*;

#[derive(Clone, Debug)]
pub struct StateVector<System: DynamicalSystem + ?Sized> {
    pub data: DVector<f64>,
    _phantom: PhantomData<System>
}

impl<System: DynamicalSystem> StateVector<System>  {
    pub fn new(x: DVector<f64>) -> Self {
        assert!(x.len() == System::STATE_VECTOR_SIZE);
        StateVector { data: x, _phantom: PhantomData }
    }
}
//refazer o from dvector
pub trait IntoSV {
    fn into_sv<System: DynamicalSystem>(self) -> StateVector<System>;
}

impl<const N: usize> IntoSV for [f64; N] {
    fn into_sv<System: DynamicalSystem>(self) -> StateVector<System> {
        assert!(N == System::STATE_VECTOR_SIZE);
        StateVector { data: DVector::from_iterator(N, self.into_iter()), _phantom: PhantomData }
    }
}

impl IntoSV for &[f64] {
    fn into_sv<System: DynamicalSystem>(self) -> StateVector<System> {
        assert!(self.len() == System::STATE_VECTOR_SIZE);
        StateVector { data: DVector::from_row_slice(self), _phantom: PhantomData }
    }
}

use std::ops::Index;

impl<System: DynamicalSystem> Index<usize> for StateVector<System> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }


}