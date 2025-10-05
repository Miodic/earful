use crate::signal::{Sampler, Signal};
use crate::types::Frequency;
use std::fs::{self, File};
use std::io::{Result, Write};
use std::path::Path;
use std::process::Command;

pub struct AudioChannel<const SAMPLES: usize> {
    samples: Box<[f64; SAMPLES]>,
}

impl<const SAMPLES: usize> AudioChannel<SAMPLES> {
    pub fn new() -> Self {
        Self {
            samples: vec![0.0; SAMPLES].into_boxed_slice().try_into().unwrap(),
        }
    }
}

pub struct Audio<const SAMPLES: usize> {
    pub left_channel: AudioChannel<SAMPLES>,
    pub right_channel: AudioChannel<SAMPLES>,
    sample_rate: Frequency,
}

impl<const SAMPLES: usize> Audio<SAMPLES> {
    pub fn new(sample_rate: Frequency) -> Self {
        Self {
            left_channel: AudioChannel::new(),
            right_channel: AudioChannel::new(),
            sample_rate,
        }
    }

    /// Add a signal to both channels with specified amplitude
    ///
    /// The signal is sampled with RELATIVE indices starting from 0,
    /// making signals position-independent in the audio buffer.
    pub fn add(&mut self, signal: &Signal, start_sample: usize, end_sample: usize, amplitude: f64) {
        self.add_left(signal, start_sample, end_sample, amplitude);
        self.add_right(signal, start_sample, end_sample, amplitude);
    }

    /// Add a signal with amplitude to the right channel only
    ///
    /// The signal is sampled with RELATIVE indices starting from 0.
    pub fn add_right(
        &mut self,
        signal: &Signal,
        start_sample: usize,
        end_sample: usize,
        amplitude: f64,
    ) {
        let end = end_sample.min(SAMPLES);
        for i in start_sample..end {
            // Pass RELATIVE sample index to the signal
            let relative_index = i - start_sample;
            self.right_channel.samples[i] +=
                signal.sample(relative_index, self.sample_rate) * amplitude;
        }
    }

    /// Add a signal with amplitude to the left channel only
    ///
    /// The signal is sampled with RELATIVE indices starting from 0.
    pub fn add_left(
        &mut self,
        signal: &Signal,
        start_sample: usize,
        end_sample: usize,
        amplitude: f64,
    ) {
        let end = end_sample.min(SAMPLES);
        for i in start_sample..end {
            // Pass RELATIVE sample index to the signal
            let relative_index = i - start_sample;
            self.left_channel.samples[i] +=
                signal.sample(relative_index, self.sample_rate) * amplitude;
        }
    }

    fn to_interleaved_f32_pcm(&self) -> Vec<f32> {
        let mut interleaved = Vec::with_capacity(SAMPLES * 2);
        for i in 0..SAMPLES {
            interleaved.push(self.left_channel.samples[i] as f32);
            interleaved.push(self.right_channel.samples[i] as f32);
        }
        interleaved
    }

    pub fn save(&self, filename: &str) -> Result<()> {
        let mut file = File::create(filename)?;
        let interleaved_data = self.to_interleaved_f32_pcm();

        // WAV header constants for 32-bit float stereo
        let channels: u16 = 2;
        let bits_per_sample: u16 = 32;
        let sample_rate_u32 = match self.sample_rate.as_hz().round() {
            r if r.is_finite() && r >= 0.0 && r <= u32::MAX as f64 => r as u32,
            _ => panic!("FATAL: Invalid sample rate: not finite or out of range for u32"),
        };
        let byte_rate = sample_rate_u32 * channels as u32 * bits_per_sample as u32 / 8;
        let block_align = channels * bits_per_sample / 8;
        let pcm_data_size = interleaved_data.len() * 4;
        let file_size = 36 + pcm_data_size;

        // Write WAV header
        file.write_all(b"RIFF")?;
        file.write_all(&(file_size as u32).to_le_bytes())?;
        file.write_all(b"WAVE")?;
        file.write_all(b"fmt ")?;
        file.write_all(&16u32.to_le_bytes())?;
        file.write_all(&3u16.to_le_bytes())?; // IEEE Float = 3
        file.write_all(&channels.to_le_bytes())?;
        file.write_all(&sample_rate_u32.to_le_bytes())?;
        file.write_all(&byte_rate.to_le_bytes())?;
        file.write_all(&block_align.to_le_bytes())?;
        file.write_all(&bits_per_sample.to_le_bytes())?;
        file.write_all(b"data")?;
        file.write_all(&(pcm_data_size as u32).to_le_bytes())?;

        // Write 32-bit float PCM data
        for &sample in &interleaved_data {
            file.write_all(&sample.to_le_bytes())?;
        }

        Ok(())
    }

    pub fn generate_visualizations(wav_file_path: &str, output_dir: &str) {
        // Configuration flags for each visualization type
        const ENABLE_SHOWWAVES: bool = true;
        const ENABLE_SHOWSPECTRUM: bool = true;
        const ENABLE_SHOWCQT: bool = true;
        const ENABLE_SHOWFREQS: bool = true;
        const ENABLE_SHOWSPECTRUMPIC: bool = true;
        const ENABLE_SHOWWAVESPIC: bool = true;
        const ENABLE_AVECTORSCOPE: bool = true;
        const ENABLE_APHASEMETER: bool = true;
        const ENABLE_AHISTOGRAM: bool = true;
        const ENABLE_SHOWVOLUME: bool = true;

        // Create output directory if it doesn't exist
        fs::create_dir_all(output_dir).expect("Failed to create output directory");
        println!("Output directory: {}", output_dir);

        println!("Generating video visualizations with FFmpeg...");

        // showwaves - Animated waveform
        if ENABLE_SHOWWAVES {
            println!("Creating showwaves visualization...");
            let output_path = Path::new(output_dir).join("output_showwaves.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]showwaves=mode=cline:s=904x904:colors=White,format=yuv420p[v]")
                .arg("-map")
                .arg("[v]")
                .arg("-map")
                .arg("0:a")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for showwaves");
            println!("showwaves exited with status: {}", status);
        }

        // showspectrum - Frequency spectrum over time
        if ENABLE_SHOWSPECTRUM {
            println!("Creating showspectrum visualization...");
            let output_path = Path::new(output_dir).join("output_showspectrum.mp4");
            let status = Command::new("ffmpeg")
            .arg("-i")
            .arg(wav_file_path)
            .arg("-filter_complex")
            .arg("[0:a]showspectrum=s=1920x1080:mode=combined:slide=scroll:saturation=0.2:scale=log,format=yuv420p[v]")
            .arg("-map")
            .arg("[v]")
            .arg("-map")
            .arg("0:a")
            .arg("-y")
            .arg(output_path)
            .status()
            .expect("Failed to execute FFmpeg for showspectrum");
            println!("showspectrum exited with status: {}", status);
        }

        // showcqt - Constant-Q transform visualization
        if ENABLE_SHOWCQT {
            println!("Creating showcqt visualization...");
            let output_path = Path::new(output_dir).join("output_showcqt.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]showcqt=s=1920x1080,format=yuv420p[v]")
                .arg("-map")
                .arg("[v]")
                .arg("-map")
                .arg("0:a")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for showcqt");
            println!("showcqt exited with status: {}", status);
        }

        // showfreqs - Audio power spectrum
        if ENABLE_SHOWFREQS {
            println!("Creating showfreqs visualization...");
            let output_path = Path::new(output_dir).join("output_showfreqs.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]showfreqs=s=1920x1080:mode=line:fscale=log,format=yuv420p[v]")
                .arg("-map")
                .arg("[v]")
                .arg("-map")
                .arg("0:a")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for showfreqs");
            println!("showfreqs exited with status: {}", status);
        }

        // showspectrumpic - Static spectrogram image
        if ENABLE_SHOWSPECTRUMPIC {
            println!("Creating showspectrumpic visualization...");
            let output_path = Path::new(output_dir).join("output_showspectrumpic.png");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-lavfi")
                .arg("showspectrumpic=s=1920x1080:legend=1")
                .arg("-update")
                .arg("1")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for showspectrumpic");
            println!("showspectrumpic exited with status: {}", status);
        }

        // showwavespic - Static waveform image
        if ENABLE_SHOWWAVESPIC {
            println!("Creating showwavespic visualization...");
            let output_path = Path::new(output_dir).join("output_showwavespic.png");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("aformat=channel_layouts=mono,showwavespic=s=1920x1080:colors=white")
                .arg("-frames:v")
                .arg("1")
                .arg("-update")
                .arg("1")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for showwavespic");
            println!("showwavespic exited with status: {}", status);
        }

        // avectorscope - Stereo vectorscope
        if ENABLE_AVECTORSCOPE {
            println!("Creating avectorscope visualization...");
            let output_path = Path::new(output_dir).join("output_avectorscope.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]avectorscope=s=1920x1080:zoom=1.5:rc=0:gc=200:bc=0,format=yuv420p[v]")
                .arg("-map")
                .arg("[v]")
                .arg("-map")
                .arg("0:a")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for avectorscope");
            println!("avectorscope exited with status: {}", status);
        }

        // aphasemeter - Stereo phase correlation meter
        if ENABLE_APHASEMETER {
            println!("Creating aphasemeter visualization...");
            let output_path = Path::new(output_dir).join("output_aphasemeter.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]aphasemeter=s=1920x1080:mpc=cyan:video=1[out0][out1]")
                .arg("-map")
                .arg("[out1]")
                .arg("-map")
                .arg("[out0]")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for aphasemeter");
            println!("aphasemeter exited with status: {}", status);
        }

        // ahistogram - Frequency distribution histogram
        if ENABLE_AHISTOGRAM {
            println!("Creating ahistogram visualization...");
            let output_path = Path::new(output_dir).join("output_ahistogram.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]ahistogram=s=1920x1080:slide=scroll:scale=log,format=yuv420p[v]")
                .arg("-map")
                .arg("[v]")
                .arg("-map")
                .arg("0:a")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for ahistogram");
            println!("ahistogram exited with status: {}", status);
        }

        // showvolume - Volume level meter
        if ENABLE_SHOWVOLUME {
            println!("Creating showvolume visualization...");
            let output_path = Path::new(output_dir).join("output_showvolume.mp4");
            let status = Command::new("ffmpeg")
                .arg("-i")
                .arg(wav_file_path)
                .arg("-filter_complex")
                .arg("[0:a]showvolume=f=0.5:c=VOLUME:b=4:w=1920:h=900,format=yuv420p[v]")
                .arg("-map")
                .arg("[v]")
                .arg("-map")
                .arg("0:a")
                .arg("-y")
                .arg(output_path)
                .status()
                .expect("Failed to execute FFmpeg for showvolume");
            println!("showvolume exited with status: {}", status);
        }

        println!("All enabled visualizations completed!");
    }
}
