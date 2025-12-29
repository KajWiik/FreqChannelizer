# FreqChannelizer.jl

A Julia package implementing an Arbitrary Frequency Polyphase Channelizer, based on the method described by Fred Harris in [*"Polyphase analysis filter bank down-converts unequal channel bandwidths with arbitrary center frequencies"*](https://www.researchgate.net/publication/257514234_Polyphase_analysis_filter_bank_down-converts_unequal_channel_bandwidths_with_arbitrary_center_frequencies).

## Overview

This package provides an efficient structure for simultaneously down-converting multiple signals with arbitrary center frequencies and equal bandwidths. It uses a polyphase filter bank for efficient filtering and decimation, combined with complex heterodyning (frequency rotation) to correct for the frequency offsets resulting from aliasing.

## Features

-   **Arbitrary Center Frequencies**: Channels can be placed anywhere in the Nyquist zone, they do not need to be aligned with the standard DFT bin centers ($k \cdot f_{\mathrm{s}}/M$).
-   **Polyphase Efficiency**: Uses an $M$-path polyphase filter bank for computational efficiency.
-   **Baseband Correction**: Automatically applies complex frequency rotation to shift the aliased outputs to true DC (0 Hz) baseband.
-   **Configurable**: Supports adjustable number of polyphase branches ($M$), channel bandwidth, and filter tap length.
-   **Lightweight**: Only depends on `DSP.jl`.

## Installation

You can install this package directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/KajWiik/FreqChannelizer.jl.git")
```

## Usage

### Basic Example

```julia
using FreqChannelizer

# System Parameters
fs = 10.0e6           # Input sampling rate (10 MHz)
channel_bw = 200.0e3  # Target channel bandwidth (200 kHz)
M = 32                # Number of polyphase branches (decimation factor)

# Define arbitrary center frequencies for the channels
center_freqs = [0.5e6, 1.25e6, -2.1e6]

# Create the channelizer
channelizer = ArbitraryFrequencyChannelizer(center_freqs, fs, channel_bw, M)

# Create a complex input signal (e.g., from an SDR)
input_signal = rand(ComplexF64, 1024 * M) 

# Process the signal
# Returns a Matrix of size (num_channels x N_output_samples)
output_channels = channelize(channelizer, input_signal)

# output_channels[1, :] contains the baseband signal for 0.5 MHz
# output_channels[2, :] contains the baseband signal for 1.25 MHz
# ...
```

## Algorithm Details

The implementation follows these key steps:
1.  **Prototype Filter**: A Kaiser-windowed sinc lowpass filter is designed based on the target channel bandwidth.
2.  **Modulation**: The prototype filter is modulated to each target center frequency $f_c$.
3.  **Polyphase Decomposition**: The modulated filters are decomposed into $M$ polyphase branches.
4.  **Channelization**: The input signal is fed into a commutator (distributing samples to branches with appropriate delays) and convolved with the polyphase filters.
5.  **Frequency Derotation**: The decimated output, which is aliased to baseband, is multiplied by a complex exponential $e^{-j 2\pi f_{\mathrm{offset}} n}$ to shift the signal center perfectly to DC.

## Example
[`sweep.jl`](examples/sweep.jl) is a demonstration and validation script for the `FreqChannelizer` package. It simulates a real-time wideband digital down-conversion (DDC) environment and visualizes the results.

The example 
* creates a frequency sweep (chirp) test signal with center frequencies from **0 to 4 GHz** at an **8 GS/s** sampling rate.
* configures the polyphase channelizer with **8 randomly placed 128 MHz bands**, demonstrating the ability to extract non-uniformly spaced channels.

The script displays an animation of
*   **Wideband Spectrum**: The input signal moving through the passbands of the fixed channel filters.
*   **Channel Grid**: Individual subplots showing each channel's output, verifying that signals are correctly shifted to baseband (DC) as the sweep passes through their frequency range.

![](/examples/sweep.gif)

## Dependencies

-   **DSP.jl**: For filter window generation.

## License

MIT License
