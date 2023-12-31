use control_systems::{DynamicalSystem, DVector, dvector, UnityFeedback, StateVector, IntoSV, NegativeFeedback};

#[derive(Clone, Copy)]
struct SecondOrder {
    k: f64,
    c: f64
}

impl DynamicalSystem for SecondOrder {
    fn xdot(&self, _t: f64, 
        x: &StateVector<Self>, 
        u: DVector<f64>) -> DVector<f64> {
        dvector![
            x[1],
            -self.k*x[0]-self.c*x[1]+u[0]
        ]
    }

    fn y(&self, _t: f64, 
        x: &StateVector<Self>, 
        _u: DVector<f64>) -> DVector<f64> {
        dvector![x[0]]
    }

    const STATE_VECTOR_SIZE: usize = 2;

    const INPUT_SIZE      : usize = 1;

    const OUTPUT_SIZE     : usize = 1;
}

fn main() {
    let test = SecondOrder{ k: 1.0, c: 1.0 };

    let feedback1 = UnityFeedback::<SecondOrder, 1>::new_unity_feedback( test);

    let output = feedback1.y(0.0, &[0.0, 1.0, 0.0, 1.0].into_sv::<UnityFeedback<SecondOrder, 1>>(), dvector![1.0]);

    println!("{}", output);
}