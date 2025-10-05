# earful

Learning about sound syntheisis: All roads lead to gaussians.

## Installation:
You need [`ffmpeg`](https://ffmpeg.org/download.html) and some way to play `.wav` audio files on your system.

## Usage:
```bash
cargo run --release
```
This generates a raw `output.wav` file and visual artefacts in `out_ffmpeg/`

## Structure
### main.rs
Entry point, currently creates a single audio-file that concatenates the three example signals.

### audio.rs
This is the central type that you interact with to build your audio file.
It allows you to `add` new signals and `save` the result as a `.wav`.

It also provides a set of useful ffmpeg commands to visualize the generated audio, which is quite useful for debugging.

### signal.rs
Here you can define new signals. Each signal needs to implement the Sampler trait and be registered in the Signal enum.

signal currently implements three examples:
1. a pure sine of a given frequency.
2. an example steep triangle PCHIP-spline with cubic Hermite coefficients.
3. an example gaussian mixture with individual exponential decays.

Each uses a different envelope:
1. Sine is simple enough to implement `sample` directly.
2. the spline is based on the `pchip_spline_envelope.rs`.
3. the gaussian mixture is based on the `gaussian_mixture_envelope.rs`.

When adding a new envelope, always create an example signal that demonstrates how to use it.

### types.rs
If you create re-usable types, put them into `types/`, otherwise create a local module.

### envelopes
TODO