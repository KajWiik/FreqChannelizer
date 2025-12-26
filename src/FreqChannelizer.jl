module FreqChannelizer

using FFTW
using DSP

export ArbitraryFrequencyChannelizer, channelize, create_channel_filters, design_prototype_filter, generate_multitone_signal

"""
    ArbitraryFrequencyChannelizer

Polyphase channelizer with arbitrary center frequencies.
Based on Fred Harris paper: "Polyphase analysis filter bank down-converts 
unequal channel bandwidths with arbitrary center frequencies"

# Fields
- `M::Int`: Number of polyphase branches
- `num_channels::Int`: Number of output channels
- `center_freqs::Vector{Float64}`: Center frequency for each channel (Hz)
- `fs::Float64`: Input sampling rate (Hz)
- `polyphase_filters::Vector{Vector{Vector{ComplexF64}}}`: [channel][branch][tap]
- `channel_bw::Float64`: Bandwidth per channel (Hz)
- `buffers::Vector{Vector{ComplexF64}}`: M buffers (one per polyphase branch)
- `phase_accumulators::Vector{Float64}`: Phase accumulators for frequency correction
"""
mutable struct ArbitraryFrequencyChannelizer
    M::Int
    num_channels::Int
    center_freqs::Vector{Float64}
    fs::Float64
    polyphase_filters::Vector{Vector{Vector{ComplexF64}}}
    channel_bw::Float64
    buffers::Vector{Vector{ComplexF64}}
    phase_accumulators::Vector{Float64}
end

"""
    ArbitraryFrequencyChannelizer(center_freqs, fs, channel_bw, M=32, taps_per_branch=32)

Create polyphase channelizer with arbitrary center frequencies following Harris method.

# Arguments
- `center_freqs::Vector{Float64}`: Center frequency for each channel (Hz)
- `fs::Float64`: Input sampling rate (Hz)  
- `channel_bw::Float64`: Desired bandwidth per channel (Hz)
- `M::Int`: Number of polyphase branches (default: 32)
- `taps_per_branch::Int`: Filter taps per branch (default: 32)
"""
function ArbitraryFrequencyChannelizer(
    center_freqs::Vector{Float64},
    fs::Float64,
    channel_bw::Float64,
    M::Int = 32,
    taps_per_branch::Int = 32,
)
    num_channels = length(center_freqs)
    L = M * taps_per_branch

    # Design prototype lowpass filter
    prototype = design_prototype_filter(channel_bw, L, fs)

    # Create frequency-shifted polyphase filters for each channel
    # polyphase_filters[ch][br] = complex coefficients for channel ch, branch br
    polyphase_filters = create_channel_filters(prototype, M, center_freqs, fs)

    # Pre-allocate buffers (one per branch)
    buffers = [zeros(ComplexF64, taps_per_branch) for _ = 1:M]

    # Initialize phase accumulators
    phase_accumulators = zeros(Float64, num_channels)

    ArbitraryFrequencyChannelizer(
        M,
        num_channels,
        center_freqs,
        fs,
        polyphase_filters,
        channel_bw,
        buffers,
        phase_accumulators
    )
end

"""
    design_prototype_filter(bw, L, fs)

Design Nyquist prototype lowpass filter.
"""
function design_prototype_filter(bw::Float64, L::Int, fs::Float64)
    # Normalized cutoff: half the channel bandwidth
    fc = bw / (2 * fs)

    # println("  Prototype filter: fc = $fc ($(bw/2/1e3) kHz cutoff)")

    beta = 7.0
    w = kaiser(L, beta)

    # Ideal lowpass using sinc
    n = 0:(L-1)
    center = (L - 1) / 2

    h = zeros(Float64, L)
    for i = 1:L
        if n[i] == center
            h[i] = 2 * fc
        else
            h[i] = sin(2π * fc * (n[i] - center)) / (π * (n[i] - center))
        end
    end

    # Window and normalize
    h_windowed = h .* w
    h_normalized = h_windowed ./ sum(h_windowed)

    # Scale for better numerical performance
    h_scaled = h_normalized .* L

    # println("  Prototype: length=$L, sum=$(sum(h_scaled)), max=$(maximum(abs.(h_scaled)))")

    return h_scaled
end

"""
    create_channel_filters(h, M, center_freqs, fs)

Create modulated polyphase filter banks for each channel.

Following Harris: For channel k at frequency fc_k:
- Modulate prototype: h_mod[n] = h[n] * exp(-j*2π*fc_k*n/fs)
- Decompose to M branches: branch m gets h_mod[m], h_mod[m+M], h_mod[m+2M], ...

Returns: [channel][branch][tap] array
"""
function create_channel_filters(
    h::Vector{Float64},
    M::Int,
    center_freqs::Vector{Float64},
    fs::Float64,
)
    L = length(h)
    taps_per_branch = div(L, M)
    num_channels = length(center_freqs)

    # Initialize: polyphase_filters[channel][branch] = coefficients
    polyphase_filters = [[ComplexF64[] for _ = 1:M] for _ = 1:num_channels]

    for ch = 1:num_channels
        fc = center_freqs[ch]

        # Modulate prototype filter to channel frequency
        # h_ch[n] = h[n] * exp(j*2π*fc*n/fs)
        h_modulated = zeros(ComplexF64, L)
        for n = 0:(L-1)
            h_modulated[n+1] = h[n+1] * exp(2π * im * fc * n / fs)
        end

        # Polyphase decomposition: branch m gets samples m, m+M, m+2M, ...
        for i = 1:L
            branch_idx = mod(i - 1, M) + 1
            push!(polyphase_filters[ch][branch_idx], h_modulated[i])
        end
    end

    return polyphase_filters
end

"""
    channelize(channelizer, input_signal)

Channelize using polyphase structure with arbitrary frequencies.

Following Harris paper structure:
1. Commutator: M input samples → M polyphase branches
2. Each branch: filter with its buffer
3. For each channel: sum all M branch outputs with channel's complex weights

# Arguments
- `channelizer::ArbitraryFrequencyChannelizer`: The channelizer
- `input_signal::Vector{ComplexF64}`: Input at fs

# Returns  
- `Matrix{ComplexF64}`: Channels (num_channels × N_out) at fs/M
"""
function channelize(
    channelizer::ArbitraryFrequencyChannelizer,
    input_signal::Vector{ComplexF64},
)
    M = channelizer.M
    num_channels = channelizer.num_channels
    N = length(input_signal)
    N_out = div(N, M)

    # Output: num_channels × N_out samples at fs/M
    output = zeros(ComplexF64, num_channels, N_out)

    # Process M samples at a time
    for block = 1:N_out
        # Get M input samples
        start_idx = (block - 1) * M + 1
        input_block = input_signal[start_idx:(start_idx+M-1)]

        # Update polyphase branch buffers (commutator distributes samples)
        # Branch m (1-based) corresponds to delay index r = m-1
        # It should receive sample x[now - r]
        # input_block[M] is x[now], input_block[1] is x[now-(M-1)]
        # So Branch m gets input_block[M - (m-1)] = input_block[M - m + 1]
        for m in 1:M
            circshift!(channelizer.buffers[m], 1)
            channelizer.buffers[m][1] = input_block[M - m + 1]
        end

        # For each output channel
        for ch = 1:num_channels
            y = zero(ComplexF64)

            # Sum contributions from all M polyphase branches
            for m = 1:M
                # Convolve: branch m buffer with channel ch's branch m coefficients
                for k = 1:length(channelizer.polyphase_filters[ch][m])
                    y += channelizer.buffers[m][k] * channelizer.polyphase_filters[ch][m][k]
                end
            end

            # Apply frequency derotation to shift alias to baseband (DC)
            # Output sample rate is fs/M
            # We need to shift by -fc
            rot_freq = channelizer.center_freqs[ch]
            output_rate = channelizer.fs / M
            
            # Phase increment per sample
            dphi = -2π * (rot_freq / output_rate)
            
            # Update accumulator and rotate
            channelizer.phase_accumulators[ch] += dphi
             # Keep angle in [-pi, pi] for numerical stability (optional but good practice)
            # channelizer.phase_accumulators[ch] = angle(exp(im * channelizer.phase_accumulators[ch]))
            
            y *= exp(im * channelizer.phase_accumulators[ch])

            output[ch, block] = y
        end
    end

    return output
end

"""
    generate_multitone_signal(fs, duration, frequencies, amplitudes=nothing)

Generate test signal.
"""
function generate_multitone_signal(
    fs::Float64,
    duration::Float64,
    frequencies::Vector{Float64},
    amplitudes::Union{Vector{Float64},Nothing} = nothing,
)
    t = 0:(1/fs):(duration-1/fs)
    signal = zeros(ComplexF64, length(t))

    if amplitudes === nothing
        amplitudes = ones(length(frequencies))
    end

    for (f, a) in zip(frequencies, amplitudes)
        signal .+= a .* exp.(2π * im * f .* t)
    end

    signal ./ sqrt(sum(abs2.(amplitudes)))
end

end # module
