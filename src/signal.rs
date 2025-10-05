pub mod sine;
pub use sine::Sine;
pub mod example_spline;
pub use example_spline::ExampleSpectral;
pub mod example_gaussian;
use crate::types::{Amplitude, Frequency};
pub use example_gaussian::ExampleGaussianSpectrum;

pub trait Sampler {
    fn sample(&self, sample_index: usize, sample_rate: Frequency) -> Amplitude;
}

/// Signal is the interface that `audio` uses to obtain samples from.
pub enum Signal {
    Sine(Sine),
    ExampleSpectral(ExampleSpectral),
    ExampleGaussianSpectrum(ExampleGaussianSpectrum),
}

// Update the Sampler implementation
impl Sampler for Signal {
    fn sample(&self, sample_index: usize, sample_rate: Frequency) -> Amplitude {
        match self {
            Signal::Sine(sine) => sine.sample(sample_index, sample_rate),
            Signal::ExampleSpectral(example_spectral) => {
                example_spectral.sample(sample_index, sample_rate)
            }
            Signal::ExampleGaussianSpectrum(example_gaussian) => {
                example_gaussian.sample(sample_index, sample_rate)
            }
        }
    }
}
