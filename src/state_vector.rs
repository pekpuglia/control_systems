use crate::*;

#[derive(Clone, Debug)]
pub struct StateVector<System: DynamicalSystem> {
    pub data: DVector<f64>,
    _phantom: PhantomData<System>
}

impl<System: DynamicalSystem> StateVector<System>  {
    pub fn new(x: DVector<f64>) -> Self {
        assert!(x.len() == System::STATE_VECTOR_SIZE);
        StateVector { data: x, _phantom: PhantomData }
    }
}

impl<System: DynamicalSystem> From<DVector<f64>> for StateVector<System> {
    fn from(value: DVector<f64>) -> Self {
        StateVector { data: value, _phantom: PhantomData }
    }
}

impl<System: DynamicalSystem> From<&[f64]> for StateVector<System>  {
    fn from(value: &[f64]) -> Self {
        assert!(value.len() == System::STATE_VECTOR_SIZE);
        StateVector { data: DVector::from_row_slice(value), 
            _phantom: PhantomData }
    }
}

impl<System: DynamicalSystem, const N: usize> From<[f64; N]> for StateVector<System> {
    fn from(value: [f64; N]) -> Self {
        assert!(N == System::STATE_VECTOR_SIZE);
        StateVector { data: DVector::from_iterator(N, value.into_iter()), _phantom: PhantomData }
    }
}

pub trait IntoSV {
    fn into_sv<System: DynamicalSystem>(self) -> StateVector<System>;
}

impl<const N: usize> IntoSV for [f64; N] {
    fn into_sv<System: DynamicalSystem>(self) -> StateVector<System> {
        assert!(N == System::STATE_VECTOR_SIZE);
        self.into()
    }
}