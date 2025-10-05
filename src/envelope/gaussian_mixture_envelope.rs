pub mod gaussian;
pub use gaussian::Gaussian;

use crate::types::Frequency;

pub struct SpectralEnvelopeGaussians {
    magnitude_gaussians: Vec<Gaussian>,
    frequency_range: Option<(Frequency, Frequency)>,
}

impl SpectralEnvelopeGaussians {
    pub fn new() -> Self {
        Self {
            magnitude_gaussians: Vec::new(),
            frequency_range: None,
        }
    }

    /// Add a magnitude Gaussian normalized to produce a specific integral
    pub fn add_normalized_magnitude_gaussian(
        &mut self,
        center: Frequency,
        sigma: Frequency,
        target_integral: f64,
        phase: f64,
        decay_rate: f64,
    ) {
        let gaussian = Gaussian::with_integral(
            center.as_hz(),
            sigma.as_hz(),
            target_integral,
            phase,
            decay_rate,
        );
        self.magnitude_gaussians.push(gaussian);
        self.update_range(center, sigma);
    }

    /// Update the frequency range based on a new Gaussian (±3σ coverage)
    fn update_range(&mut self, center: Frequency, sigma: Frequency) {
        let span = 3.0; // ±3σ covers 99.7% of energy
        let min_freq = Frequency::from_hz((center.as_hz() - span * sigma.as_hz()).max(0.0));
        let max_freq = Frequency::from_hz(center.as_hz() + span * sigma.as_hz());

        match self.frequency_range {
            None => {
                self.frequency_range = Some((min_freq, max_freq));
            }
            Some((current_min, current_max)) => {
                let new_min = if min_freq.as_hz() < current_min.as_hz() {
                    min_freq
                } else {
                    current_min
                };
                let new_max = if max_freq.as_hz() > current_max.as_hz() {
                    max_freq
                } else {
                    current_max
                };
                self.frequency_range = Some((new_min, new_max));
            }
        }
    }

    /// **ANALYTICAL TIME-DOMAIN SAMPLING**
    ///
    /// Compute the time-domain signal directly using the closed-form
    /// inverse Fourier transform of Gaussian spectra.
    ///
    /// Each Gaussian component contributes analytically with its own
    /// phase and decay rate.
    pub fn sample_time_domain(&self, t: f64) -> f64 {
        self.magnitude_gaussians
            .iter()
            .map(|gaussian| gaussian.evaluate_time_domain(t))
            .sum()
    }
}
