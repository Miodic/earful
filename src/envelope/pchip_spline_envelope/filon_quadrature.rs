use std::f64::consts::PI;

use crate::envelope::pchip_spline_envelope::Complex;

pub struct FilonQuadrature {}

impl FilonQuadrature {
    pub fn new() -> Self {
        Self {}
    }

    /// Compute ∫_a^b F(ω) * e^(i*ω*t) dω using Filon-type method
    /// This is specialized for inverse Fourier transform
    pub fn integrate_oscillatory<F, DF>(
        &self,
        integrand: F,
        derivative: DF,
        a: f64,
        b: f64,
        oscillation_param: f64, // This is 't' in e^(i*ω*t)
    ) -> Complex
    where
        F: Fn(f64) -> Complex,
        DF: Fn(f64) -> Complex,
    {
        if (b - a).abs() < 1e-12 {
            return Complex::new(0.0, 0.0);
        }

        // For Filon method, we use adaptive quadrature based on oscillation
        // Split into subintervals if highly oscillatory
        let omega_range = (b - a) * oscillation_param.abs();

        if omega_range < 10.0 {
            // Low oscillation: use standard Gauss-Legendre quadrature
            self.gauss_legendre_complex(&integrand, a, b, oscillation_param)
        } else {
            // High oscillation: use Hermite Filon-type method with subdivision
            let num_subdivisions = (omega_range / (2.0 * PI)).ceil() as usize;
            let num_subdivisions = num_subdivisions.max(4).min(100);

            let mut result = Complex::new(0.0, 0.0);
            let h = (b - a) / num_subdivisions as f64;

            for i in 0..num_subdivisions {
                let sub_a = a + i as f64 * h;
                let sub_b = sub_a + h;
                let sub_result = self.hermite_filon_method(
                    &integrand,
                    &derivative,
                    sub_a,
                    sub_b,
                    oscillation_param,
                );
                result = result.add(&sub_result);
            }

            result
        }
    }

    /// Hermite Filon-type integration for ∫_a^b F(ω) * e^(i*ω*t) dω
    /// Uses cubic Hermite interpolation of F(ω) and exact integration of oscillatory part
    fn hermite_filon_method<F, DF>(
        &self,
        integrand: &F,
        derivative: &DF,
        a: f64,
        b: f64,
        t: f64,
    ) -> Complex
    where
        F: Fn(f64) -> Complex,
        DF: Fn(f64) -> Complex,
    {
        // Evaluate function and derivatives at endpoints
        let fa = integrand(a);
        let fb = integrand(b);
        let dfa = derivative(a);
        let dfb = derivative(b);
        let h = b - a;
        let theta = h * t;

        // For small oscillation parameter, use Simpson's rule
        if theta.abs() < 0.1 {
            let mid = (a + b) / 2.0;
            let fm = integrand(mid);

            let exp_a = Complex::from_polar(1.0, a * t);
            let exp_m = Complex::from_polar(1.0, mid * t);
            let exp_b = Complex::from_polar(1.0, b * t);

            let h_sixth = h / 6.0;
            let contrib_a = fa.mul(&exp_a).scale(h_sixth);
            let contrib_m = fm.mul(&exp_m).scale(4.0 * h_sixth);
            let contrib_b = fb.mul(&exp_b).scale(h_sixth);

            return contrib_a.add(&contrib_m).add(&contrib_b);
        }

        // Compute moment integrals for Hermite basis functions
        // M_k = ∫_a^b (ω - a)^k * e^(i*ω*t) dω for k = 0, 1, 2, 3
        let moments = self.compute_hermite_moments(a, b, t);

        // Hermite cubic interpolant:
        // H(ω) = fa * h00(s) + dfa * h * h10(s) + fb * h01(s) + dfb * h * h11(s)
        // where s = (ω - a) / h
        //
        // Hermite basis functions in terms of u = ω - a:
        // h00: coefficient is fa
        // h10: coefficient is dfa
        // h01: coefficient is fb
        // h11: coefficient is dfb
        //
        // Expanded in powers of u = (ω - a):
        // H(ω) = c0 + c1*u + c2*u² + c3*u³

        // Compute coefficients in power basis
        let h2 = h * h;
        let h3 = h2 * h;

        // h00(s) = 2s³ - 3s² + 1 = 2(u/h)³ - 3(u/h)² + 1
        //        = 2u³/h³ - 3u²/h² + 1
        // Contributes: fa * (2/h³ * u³ - 3/h² * u² + 1)

        // h10(s) = s³ - 2s² + s = (u/h)³ - 2(u/h)² + u/h
        //        = u³/h³ - 2u²/h² + u/h
        // Contributes: dfa * h * (1/h³ * u³ - 2/h² * u² + 1/h * u)
        //            = dfa * (1/h² * u³ - 2/h * u² + u)

        // h01(s) = -2s³ + 3s² = -2(u/h)³ + 3(u/h)²
        //        = -2u³/h³ + 3u²/h²
        // Contributes: fb * (-2/h³ * u³ + 3/h² * u²)

        // h11(s) = s³ - s² = (u/h)³ - (u/h)²
        //        = u³/h³ - u²/h²
        // Contributes: dfb * h * (1/h³ * u³ - 1/h² * u²)
        //            = dfb * (1/h² * u³ - 1/h * u²)

        // Coefficients for constant term (u^0)
        let c0 = fa;

        // Coefficients for linear term (u^1)
        let c1 = dfa;

        // Coefficients for quadratic term (u^2)
        let c2_fa = fa.scale(-3.0 / h2);
        let c2_dfa = dfa.scale(-2.0 / h);
        let c2_fb = fb.scale(3.0 / h2);
        let c2_dfb = dfb.scale(-1.0 / h);
        let c2 = c2_fa.add(&c2_dfa).add(&c2_fb).add(&c2_dfb);

        // Coefficients for cubic term (u^3)
        let c3_fa = fa.scale(2.0 / h3);
        let c3_dfa = dfa.scale(1.0 / h2);
        let c3_fb = fb.scale(-2.0 / h3);
        let c3_dfb = dfb.scale(1.0 / h2);
        let c3 = c3_fa.add(&c3_dfa).add(&c3_fb).add(&c3_dfb);

        // Compute integral: ∫_a^b H(ω) * e^(i*ω*t) dω
        // = c0*M0 + c1*M1 + c2*M2 + c3*M3
        let result = c0
            .mul(&moments.0)
            .add(&c1.mul(&moments.1))
            .add(&c2.mul(&moments.2))
            .add(&c3.mul(&moments.3));

        result
    }

    /// Compute moment integrals M_k = ∫_a^b (ω - a)^k * e^(i*ω*t) dω for k = 0, 1, 2, 3
    fn compute_hermite_moments(
        &self,
        a: f64,
        b: f64,
        t: f64,
    ) -> (Complex, Complex, Complex, Complex) {
        let h = b - a;
        let it = Complex::new(0.0, t); // i*t
        let exp_a = Complex::from_polar(1.0, a * t);
        let exp_b = Complex::from_polar(1.0, b * t);

        // M0 = ∫_a^b e^(i*ω*t) dω = (e^(i*b*t) - e^(i*a*t)) / (i*t)
        let m0 = exp_b.add(&exp_a.scale(-1.0)).div(&it);

        // For higher moments, use recurrence relation:
        // M_k = h^k * e^(i*b*t) / (i*t) - k * M_{k-1} / (i*t)
        // This follows from integration by parts

        // M1 = ∫_a^b (ω - a) * e^(i*ω*t) dω
        // Using substitution u = ω - a:
        // M1 = ∫_0^h u * e^(i*(a+u)*t) du = e^(i*a*t) * ∫_0^h u * e^(i*u*t) du
        // Integration by parts: ∫ u * e^(i*u*t) du = [u * e^(i*u*t) / (i*t)]_0^h - ∫ e^(i*u*t) / (i*t) du
        //                                            = h * e^(i*h*t) / (i*t) - [e^(i*u*t) / (i*t)²]_0^h
        //                                            = h * e^(i*h*t) / (i*t) - (e^(i*h*t) - 1) / (i*t)²
        let m1 = exp_b.scale(h).div(&it).add(&m0.div(&it).scale(-1.0));

        // M2 = ∫_a^b (ω - a)² * e^(i*ω*t) dω
        // Using recurrence: M2 = h² * e^(i*b*t) / (i*t) - 2 * M1 / (i*t)
        let m2 = exp_b.scale(h * h).div(&it).add(&m1.div(&it).scale(-2.0));

        // M3 = ∫_a^b (ω - a)³ * e^(i*ω*t) dω
        // Using recurrence: M3 = h³ * e^(i*b*t) / (i*t) - 3 * M2 / (i*t)
        let m3 = exp_b
            .scale(h * h * h)
            .div(&it)
            .add(&m2.div(&it).scale(-3.0));

        (m0, m1, m2, m3)
    }

    /// Standard Gauss-Legendre quadrature for complex-valued oscillatory integrals
    fn gauss_legendre_complex<F>(&self, integrand: &F, a: f64, b: f64, t: f64) -> Complex
    where
        F: Fn(f64) -> Complex,
    {
        // Gauss-Legendre nodes and weights for [-1, 1]
        let (nodes, weights) = self.gauss_legendre_nodes_weights();

        let mut result = Complex::new(0.0, 0.0);
        let half_length = (b - a) / 2.0;
        let mid = (a + b) / 2.0;

        for i in 0..nodes.len() {
            let omega = mid + half_length * nodes[i];
            let f_val = integrand(omega);
            let oscillator = Complex::from_polar(1.0, omega * t);
            let integrand_val = f_val.mul(&oscillator).scale(weights[i] * half_length);
            result = result.add(&integrand_val);
        }

        result
    }

    /// Get Gauss-Legendre nodes and weights (using 8-point rule for accuracy)
    fn gauss_legendre_nodes_weights(&self) -> (Vec<f64>, Vec<f64>) {
        // 8-point Gauss-Legendre quadrature
        let nodes = vec![
            -0.9602898564975362,
            -0.7966664774136267,
            -0.5255324099163290,
            -0.1834346424956498,
            0.1834346424956498,
            0.5255324099163290,
            0.7966664774136267,
            0.9602898564975362,
        ];

        let weights = vec![
            0.1012285362903763,
            0.2223810344533745,
            0.3137066458778873,
            0.3626837833783620,
            0.3626837833783620,
            0.3137066458778873,
            0.2223810344533745,
            0.1012285362903763,
        ];

        (nodes, weights)
    }
}
