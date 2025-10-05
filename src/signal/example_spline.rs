use crate::{
    envelope::{Complex, FilonQuadrature, PHCIPSplineEnvelope},
    signal::Sampler,
    types::{Amplitude, Frequency},
};
use std::f64::consts::PI;

pub struct ExampleSpectral {
    envelope: PHCIPSplineEnvelope,
    integrator: FilonQuadrature,
}

impl ExampleSpectral {
    pub fn new(base_frequency: Frequency) -> Self {
        let mut envelope = PHCIPSplineEnvelope::new();

        let bandwidth_factor = 1e-6;
        let half_bandwidth = base_frequency * bandwidth_factor;

        let bandwidth_hz = (2.0 * half_bandwidth).as_hz();
        let magnitude = 0.75 / bandwidth_hz;

        envelope.add_point(base_frequency - half_bandwidth, 0.0, 0.0);
        envelope.add_point(base_frequency, magnitude, 0.0);
        envelope.add_point(base_frequency + half_bandwidth, 0.0, 0.0);

        envelope.finalize();

        Self {
            envelope,
            integrator: FilonQuadrature::new(),
        }
    }
}

impl Sampler for ExampleSpectral {
    fn sample(&self, sample_index: usize, sample_rate: Frequency) -> Amplitude {
        let t = sample_index as f64 / sample_rate.as_hz();

        let (freq_min, freq_max) = self
            .envelope
            .get_range()
            .unwrap_or((Frequency::from_hz(0.0), Frequency::from_hz(0.0)));

        let spectrum_function = |omega: f64| -> Complex {
            self.envelope.sample(Frequency::from_hz(omega / (2.0 * PI)))
        };

        let spectrum_derivative = |omega: f64| -> Complex {
            // Derivative w.r.t. ω: d/dω = d/df * df/dω = d/df * 1/(2π)
            self.envelope
                .derivative(Frequency::from_hz(omega / (2.0 * PI)))
                .scale(1.0 / (2.0 * PI))
        };

        let omega_min = freq_min.as_hz() * 2.0 * PI;
        let omega_max = freq_max.as_hz() * 2.0 * PI;

        let result = self.integrator.integrate_oscillatory(
            spectrum_function,
            spectrum_derivative,
            omega_min,
            omega_max,
            t,
        );

        result.re / PI
    }
}
