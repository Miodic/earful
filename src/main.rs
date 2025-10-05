use std::time::Duration;

mod types;
use crate::{
    signal::{ExampleGaussianSpectrum, ExampleSpectral, Sine},
    types::Frequency,
};

mod audio;
use audio::Audio;

mod signal;
use signal::Signal;

mod envelope;

fn main() {
    const DURATION: Duration = Duration::from_secs(6);
    const SAMPLES: usize = 48000 * DURATION.as_secs() as usize;
    const SAMPLE_RATE: Frequency = Frequency::from_duration_and_samples(DURATION, SAMPLES);
    println!(
        "Creating {} ms of Audio over {SAMPLES} total samples at {}",
        DURATION.as_millis(),
        SAMPLE_RATE
    );

    let mut audio: Audio<SAMPLES> = Audio::new(SAMPLE_RATE);

    // A4 fundamental (440 Hz) with harmonic overtones
    // Each harmonic has decreasing amplitude for natural timbre
    let a4: Frequency = Frequency::from_hz(440.0);

    // Fundamental (A4 - 440 Hz) at full amplitude with pure Sine signal
    let sine_a4 = Sine::new(a4);
    audio.add(&Signal::Sine(sine_a4), 0, SAMPLES / 3, 1.0);

    // Fundamental (A4 - 440 Hz) at full amplitude with narrow triangular PCHIP spline
    let spectral_a4 = ExampleSpectral::new(a4);
    audio.add(
        &Signal::ExampleSpectral(spectral_a4),
        SAMPLES / 3,
        2 * SAMPLES / 3,
        1.0,
    );

    // Fundamental (A4 - 440 Hz) at full amplitude with Gaussian
    let gaussian_a4 = ExampleGaussianSpectrum::new(a4);
    audio.add(
        &Signal::ExampleGaussianSpectrum(gaussian_a4),
        2 * SAMPLES / 3,
        SAMPLES,
        1.0,
    );

    audio.save("output.wav").expect("Failed to save audio file");
    Audio::<SAMPLES>::generate_visualizations("output.wav", "out_ffmpeg");
}
