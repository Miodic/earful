use crate::{
    signal::Sampler,
    types::{Amplitude, Frequency},
};
use std::f64::consts::PI;

pub struct Sine {
    frequency: Frequency,
}

impl Sine {
    pub fn new(frequency: Frequency) -> Self {
        Self { frequency }
    }
}

impl Sampler for Sine {
    fn sample(&self, sample_index: usize, sample_rate: Frequency) -> Amplitude {
        let time = sample_index as f64 / sample_rate.as_hz();
        let angular_frequency: f64 = 2.0 * PI * self.frequency.as_hz();
        (angular_frequency * time).sin()
    }
}
