use std::{
    fmt::{self, Display},
    ops::{Add, Div, Mul, Sub},
    time::Duration,
};

#[derive(Clone, Copy, PartialEq, PartialOrd)]
pub struct Frequency {
    hz: f64,
}

impl Frequency {
    const GIGA: f64 = 1_000_000_000.;
    const MEGA: f64 = 1_000_000.;
    const KILO: f64 = 1_000.;

    #[allow(unused)]
    pub const fn from_hz(hz: f64) -> Self {
        Self { hz }
    }

    #[allow(unused)]
    pub const fn from_kilo_hz(khz: f64) -> Self {
        Self {
            hz: khz * Self::KILO,
        }
    }

    #[allow(unused)]
    pub const fn from_mega_hz(mhz: f64) -> Self {
        Self {
            hz: mhz * Self::MEGA,
        }
    }

    pub const fn from_duration_and_samples(duration: Duration, samples: usize) -> Self {
        let duration_nanos =
            duration.as_secs() * Self::GIGA as u64 + duration.subsec_nanos() as u64;
        if duration_nanos == 0 {
            panic!("Duration cannot be zero");
        }

        let hz = samples as f64 / (duration_nanos as f64 / Self::GIGA);

        Self { hz } // Adjust cast based on your Self type
    }

    pub const fn as_hz(&self) -> f64 {
        self.hz
    }
}

impl Display for Frequency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (value, unit) = match self.hz {
            x if x >= 1e9 => (x / 1e9, "GHz"),
            x if x >= 1e6 => (x / 1e6, "MHz"),
            x if x >= 1e3 => (x / 1e3, "kHz"),
            x => (x, "Hz"),
        };

        if (value.fract()).abs() < 1e-10 {
            // We can ignore the fractional part.
            write!(f, "{:.0} {}", value, unit)
        } else {
            // Remove trailing zeros after decimal point
            let formatted = format!("{:.3}", value);
            let trimmed = formatted.trim_end_matches('0').trim_end_matches('.');
            write!(f, "{} {}", trimmed, unit)
        }
    }
}

/// Adds two frequency values together.
impl Add for Frequency {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            hz: self.hz + rhs.hz,
        }
    }
}

/// Computes the difference between two frequencies.
impl Sub for Frequency {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            hz: self.hz - rhs.hz,
        }
    }
}

/// Frequency divison produces a factor.
impl Div for Frequency {
    type Output = f64;

    fn div(self, rhs: Self) -> Self::Output {
        self.hz / rhs.hz
    }
}

// Divide a factor by a frequency, resulting in time[s]
impl Div<Frequency> for f64 {
    type Output = f64;

    fn div(self, rhs: Frequency) -> Self::Output {
        self / rhs.hz
    }
}

/// Divides Frequency by a factor.
impl Div<f64> for Frequency {
    type Output = Frequency;

    fn div(self, rhs: f64) -> Self::Output {
        Self { hz: self.hz / rhs }
    }
}

/// Scales the frequency by a factor
impl Mul<f64> for Frequency {
    type Output = Self;

    fn mul(self, factor: f64) -> Self::Output {
        Self {
            hz: self.hz * factor,
        }
    }
}

/// Scales the frequency by a factor
impl Mul<Frequency> for f64 {
    type Output = Frequency;

    fn mul(self, rhs: Frequency) -> Self::Output {
        Frequency { hz: self * rhs.hz }
    }
}
