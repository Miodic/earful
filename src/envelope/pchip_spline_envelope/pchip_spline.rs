use crate::types::Frequency;

pub struct PchipSpline {
    control_points: Vec<(Frequency, f64)>,
    coefficients: Vec<[f64; 4]>, // [a, b, c, d] for each segment
    computed: bool,
}

impl PchipSpline {
    pub fn new() -> Self {
        Self {
            control_points: Vec::new(),
            coefficients: Vec::new(),
            computed: false,
        }
    }

    pub fn add_point(&mut self, frequency: Frequency, value: f64) {
        let pos = self
            .control_points
            .binary_search_by(|probe| probe.0.partial_cmp(&frequency).unwrap())
            .unwrap_or_else(|e| e);

        self.control_points.insert(pos, (frequency, value));
        self.computed = false;
    }

    /// Compute PCHIP coefficients using Fritsch-Carlson algorithm
    /// with improved endpoint derivative handling
    pub fn compute_coefficients(&mut self) {
        let n = self.control_points.len();
        if n < 2 {
            self.coefficients.clear();
            self.computed = true;
            return;
        }

        let x: Vec<f64> = self
            .control_points
            .iter()
            .map(|(freq, _)| freq.as_hz())
            .collect();
        let y: Vec<f64> = self.control_points.iter().map(|(_, val)| *val).collect();

        // Step 1: Compute secant slopes
        let mut delta = vec![0.0; n - 1];
        let mut h = vec![0.0; n - 1];
        for i in 0..n - 1 {
            h[i] = x[i + 1] - x[i];
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        // Step 2: Compute derivatives at each point
        let mut d = vec![0.0; n];

        // IMPROVED ENDPOINT DERIVATIVES
        // Use shape-preserving formula that prevents overshoot
        if n == 2 {
            // Only two points: use the slope
            d[0] = delta[0];
            d[1] = delta[0];
        } else {
            // Left endpoint: use one-sided three-point formula
            // This prevents overshoot by considering curvature
            d[0] = self.endpoint_derivative(&h, &delta, 0);

            // Right endpoint: use one-sided three-point formula
            d[n - 1] = self.endpoint_derivative_right(&h, &delta, n);
        }

        // Interior point derivatives using weighted harmonic mean
        for i in 1..n - 1 {
            if delta[i - 1] * delta[i] <= 0.0 {
                d[i] = 0.0;
            } else {
                let w1 = 2.0 * h[i] + h[i - 1];
                let w2 = h[i] + 2.0 * h[i - 1];
                d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
            }
        }

        // Step 3: Ensure monotonicity (Fritsch-Butland condition)
        for i in 0..n - 1 {
            if delta[i].abs() < 1e-10 {
                d[i] = 0.0;
                d[i + 1] = 0.0;
            } else {
                let alpha = d[i] / delta[i];
                let beta = d[i + 1] / delta[i];
                let tau_squared = alpha * alpha + beta * beta;

                if tau_squared > 9.0 {
                    let tau = 3.0 / tau_squared.sqrt();
                    d[i] = tau * alpha * delta[i];
                    d[i + 1] = tau * beta * delta[i];
                }
            }
        }

        // Step 4: Compute cubic Hermite coefficients
        self.coefficients.clear();
        for i in 0..n - 1 {
            let a = y[i];
            let b = d[i];
            let c = (3.0 * delta[i] - 2.0 * d[i] - d[i + 1]) / h[i];
            let dd = (d[i] + d[i + 1] - 2.0 * delta[i]) / (h[i] * h[i]);

            self.coefficients.push([a, b, c, dd]);
        }

        self.computed = true;
    }

    /// Compute left endpoint derivative using shape-preserving formula
    fn endpoint_derivative(&self, h: &[f64], delta: &[f64], _index: usize) -> f64 {
        // Three-point endpoint formula that respects monotonicity
        // d[0] = ((2*h[0] + h[1])*delta[0] - h[0]*delta[1]) / (h[0] + h[1])

        if delta.len() < 2 {
            return delta[0];
        }

        let d_candidate = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);

        // Ensure the derivative doesn't cause overshoot
        // If delta[0] and delta[1] have the same sign, d should be between them
        if delta[0] * delta[1] > 0.0 {
            // Clamp to range [0, 3*delta[0]] to prevent overshoot
            let min_val = 0.0_f64.min(3.0 * delta[0]);
            let max_val = 0.0_f64.max(3.0 * delta[0]);
            d_candidate.clamp(min_val, max_val)
        } else {
            // Sign change: be conservative
            if d_candidate.abs() > 3.0 * delta[0].abs() {
                delta[0]
            } else {
                d_candidate
            }
        }
    }

    /// Compute right endpoint derivative using shape-preserving formula
    fn endpoint_derivative_right(&self, h: &[f64], delta: &[f64], n: usize) -> f64 {
        if delta.len() < 2 {
            return delta[delta.len() - 1];
        }

        let n_idx = n - 2; // Index into h and delta arrays
        let d_candidate = ((2.0 * h[n_idx] + h[n_idx - 1]) * delta[n_idx]
            - h[n_idx] * delta[n_idx - 1])
            / (h[n_idx] + h[n_idx - 1]);

        // Ensure the derivative doesn't cause overshoot
        if delta[n_idx] * delta[n_idx - 1] > 0.0 {
            let min_val = 0.0_f64.min(3.0 * delta[n_idx]);
            let max_val = 0.0_f64.max(3.0 * delta[n_idx]);
            d_candidate.clamp(min_val, max_val)
        } else {
            if d_candidate.abs() > 3.0 * delta[n_idx].abs() {
                delta[n_idx]
            } else {
                d_candidate
            }
        }
    }

    pub fn evaluate(&self, frequency: Frequency) -> f64 {
        assert!(
            self.computed,
            "Spline coefficients not computed. Call compute_coefficients() first."
        );

        if self.control_points.is_empty() {
            return 0.0;
        }

        if self.control_points.len() == 1 {
            return self.control_points[0].1;
        }

        let freq_hz = frequency.as_hz();

        // Clamp to boundary values
        if freq_hz <= self.control_points[0].0.as_hz() {
            return self.control_points[0].1;
        }

        if freq_hz >= self.control_points.last().unwrap().0.as_hz() {
            return self.control_points.last().unwrap().1;
        }

        // Find the segment
        let mut segment = 0;
        for i in 0..self.control_points.len() - 1 {
            let x_i = self.control_points[i].0.as_hz();
            let x_i1 = self.control_points[i + 1].0.as_hz();
            if freq_hz >= x_i && freq_hz <= x_i1 {
                segment = i;
                break;
            }
        }

        let x_i = self.control_points[segment].0.as_hz();
        let dx = freq_hz - x_i;
        let coef = &self.coefficients[segment];

        // Evaluate cubic polynomial
        coef[0] + coef[1] * dx + coef[2] * dx * dx + coef[3] * dx * dx * dx
    }

    pub fn derivative(&self, frequency: Frequency) -> f64 {
        assert!(
            self.computed,
            "Spline coefficients not computed. Call compute_coefficients() first."
        );

        if self.control_points.len() < 2 {
            return 0.0;
        }

        let freq_hz = frequency.as_hz();

        // Return zero derivative outside range
        if freq_hz <= self.control_points[0].0.as_hz()
            || freq_hz >= self.control_points.last().unwrap().0.as_hz()
        {
            return 0.0;
        }

        // Find segment
        let mut segment = 0;
        for i in 0..self.control_points.len() - 1 {
            let x_i = self.control_points[i].0.as_hz();
            let x_i1 = self.control_points[i + 1].0.as_hz();
            if freq_hz >= x_i && freq_hz <= x_i1 {
                segment = i;
                break;
            }
        }

        let x_i = self.control_points[segment].0.as_hz();
        let dx = freq_hz - x_i;
        let coef = &self.coefficients[segment];

        // Derivative: b + 2*c*dx + 3*d*dx^2
        coef[1] + 2.0 * coef[2] * dx + 3.0 * coef[3] * dx * dx
    }

    pub fn get_range(&self) -> Option<(Frequency, Frequency)> {
        if self.control_points.is_empty() {
            None
        } else {
            Some((
                self.control_points[0].0,
                self.control_points.last().unwrap().0,
            ))
        }
    }
}
