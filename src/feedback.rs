use crate::*;

//direct dynamical system, reverse dynamical system
struct NegativeFeedback<DDS, RDS> {
    dirsys: DDS,
    revsys: RDS
}

impl<InputAndReverseOutput, OutputAndReverseInput, DDS, RDS> 
    DynamicalSystem for NegativeFeedback<DDS, RDS> 
where
    DDS: DynamicalSystem<Input = InputAndReverseOutput, Output = OutputAndReverseInput>,
    RDS: DynamicalSystem<Input = OutputAndReverseInput, Output = InputAndReverseOutput>    
{
    const STATE_VECTOR_SIZE: usize = DDS::STATE_VECTOR_SIZE + RDS::STATE_VECTOR_SIZE;

    type Input = InputAndReverseOutput;

    type Output = OutputAndReverseInput;

    fn xdot(&self, t: f64, x: DVector<f64>, u: Self::Input) -> DVector<f64> {
        todo!()
    }

    fn y(&self, t: f64, x: DVector<f64>, u: Self::Input) -> Self::Output {
        todo!()
    }
}