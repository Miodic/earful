use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct Gaussian {
    pub center: f64,
    pub sigma: f64,
    pub amplitude: f64,
    pub phase: f64,
    pub decay_rate: f64,
}

impl Gaussian {
    pub fn new(center: f64, sigma: f64, amplitude: f64, phase: f64, decay_rate: f64) -> Self {
        Self {
            center,
            sigma,
            amplitude,
            phase,
            decay_rate,
        }
    }

    /// Normalize the Gaussian to have a specified integral
    pub fn with_integral(
        center: f64,
        sigma: f64,
        target_integral: f64,
        phase: f64,
        decay_rate: f64,
    ) -> Self {
        let amplitude = target_integral / (sigma * (2.0 * PI).sqrt());
        Self::new(center, sigma, amplitude, phase, decay_rate)
    }

    /// Compute the time-domain contribution with both envelopes
    ///
    /// For complex spectrum: S(f) = A * exp(-(f - f₀)²/(2σ²)) * exp(iφ)
    /// The time-domain signal includes:
    ///   - Gaussian envelope from spectral width: exp(-2π²σ²t²)
    ///   - Exponential amplitude decay: exp(-t/τ)
    ///   - Carrier with phase: cos(2πf₀t + φ)
    pub fn evaluate_time_domain(&self, t: f64) -> f64 {
        let f0 = self.center;
        let sigma = self.sigma;
        let amplitude = self.amplitude;

        // Gaussian envelope in time domain (from spectral width)
        let gaussian_envelope =
            amplitude * sigma * (2.0 * PI).sqrt() * (-2.0 * PI * PI * sigma * sigma * t * t).exp();

        // Exponential amplitude decay envelope
        let amplitude_envelope = if self.decay_rate > 0.0 {
            (-t / self.decay_rate).exp()
        } else {
            1.0 // No decay if decay_rate is 0 or negative
        };

        // Carrier wave with phase offset
        let carrier = (2.0 * PI * f0 * t + self.phase).cos();

        gaussian_envelope * amplitude_envelope * carrier
    }
}
