use crate::{
    envelope::SpectralEnvelopeGaussians,
    signal::Sampler,
    types::{Amplitude, Frequency},
};

pub struct ExampleGaussianSpectrum {
    envelope: SpectralEnvelopeGaussians,
}

impl ExampleGaussianSpectrum {
    pub fn new(base_frequency: Frequency) -> Self {
        let mut envelope = SpectralEnvelopeGaussians::new();

        // Standard deviation for the Gaussian
        let bandwidth_factor = 0.0015;
        let sigma = Frequency::from_hz(base_frequency.as_hz() * bandwidth_factor);

        // Parameters:
        // - phase: 0.0 (zero phase)
        // - decay_rate: 2.0 seconds (exponential decay time constant)
        envelope.add_normalized_magnitude_gaussian(
            base_frequency,
            sigma,
            1.0,
            0.0, // phase
            2.0, // decay_rate (τ in seconds)
        );

        Self { envelope }
    }
}

impl Sampler for ExampleGaussianSpectrum {
    fn sample(&self, sample_index: usize, sample_rate: Frequency) -> Amplitude {
        let t = sample_index as f64 / sample_rate.as_hz();

        // analytical time-domain sampling
        self.envelope.sample_time_domain(t)
    }
}
