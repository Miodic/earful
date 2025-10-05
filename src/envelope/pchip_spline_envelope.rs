use crate::types::{Frequency, Phase};

pub mod complex;
pub use complex::Complex;

pub mod pchip_spline;
pub use pchip_spline::PchipSpline;

pub mod filon_quadrature;
pub use filon_quadrature::FilonQuadrature;

pub struct PHCIPSplineEnvelope {
    magnitude_spline: PchipSpline,
    phase_spline: PchipSpline,
}

impl PHCIPSplineEnvelope {
    pub fn new() -> Self {
        Self {
            magnitude_spline: PchipSpline::new(),
            phase_spline: PchipSpline::new(),
        }
    }

    /// Add a control point with magnitude and phase
    pub fn add_point(&mut self, frequency: Frequency, magnitude: f64, phase: Phase) {
        self.magnitude_spline.add_point(frequency, magnitude);
        self.phase_spline.add_point(frequency, phase);
    }

    /// Finalize the envelope by computing spline coefficients
    pub fn finalize(&mut self) {
        self.magnitude_spline.compute_coefficients();
        self.phase_spline.compute_coefficients();
    }

    /// Sample the complex spectrum at a given frequency
    pub fn sample(&self, frequency: Frequency) -> Complex {
        let magnitude = self.magnitude_spline.evaluate(frequency);
        let phase = self.phase_spline.evaluate(frequency);
        Complex::from_polar(magnitude, phase)
    }

    /// Compute the derivative of the complex spectrum: d/df[M(f) * e^(i*P(f))]
    /// Returns: e^(i*P) * (dM/df + i * M * dP/df)
    pub fn derivative(&self, frequency: Frequency) -> Complex {
        let magnitude = self.magnitude_spline.evaluate(frequency);
        let phase = self.phase_spline.evaluate(frequency);
        let d_magnitude = self.magnitude_spline.derivative(frequency);
        let d_phase = self.phase_spline.derivative(frequency);

        // d/df[M(f) * e^(i*P(f))] = e^(i*P) * (dM/df + i * M * dP/df)
        let derivative_magnitude = d_magnitude;
        let derivative_phase_contribution = magnitude * d_phase;

        // Convert to complex and multiply by e^(i*P)
        let deriv_polar_form = Complex::new(derivative_magnitude, derivative_phase_contribution);
        let phase_factor = Complex::from_polar(1.0, phase);

        deriv_polar_form.mul(&phase_factor)
    }

    pub fn get_range(&self) -> Option<(Frequency, Frequency)> {
        self.magnitude_spline.get_range()
    }
}
