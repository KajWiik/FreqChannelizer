using FreqChannelizer
using DSP
using FFTW
using CairoMakie
using Random

# Parameters
fs = 8.0e9             # 8 GS/s sampling rate
duration = 100.0e-6    # 100 microseconds total simulation time
chunk_duration = 1.0e-6 # 1 microsecond per animation frame
f_start = 0.0          # Sweep start frequency
f_end = 4.0e9          # Sweep end frequency
channel_bw = 128.0e6   # 128 MHz bandwidth
M = 32                 # Polyphase branches (decimation factor)
num_channels = 8

# Random center frequencies in the [0, 4 GHz] range
Random.seed!(42)
center_freqs = sort(rand(num_channels) .* (f_end - f_start) .+ f_start)

println("Initializing Channelizer...")
channelizer = ArbitraryFrequencyChannelizer(center_freqs, fs, channel_bw, M)

# Visualization Setup
fig = Figure(size = (1400, 1200), backgroundcolor = :white)
ax_main = Axis(fig[1, 1:4], 
    title = "Input Spectrum and Channel Filter Responses (8 GS/s)",
    xlabel = "Frequency (GHz)", 
    ylabel = "Magnitude (dB)",
    limits = (-0.1, 4.1, -100, 10))

# Subplots for channel outputs (4x4 grid requested)
axs_channels = Axis[]
for i in 1:num_channels
    row = (i-1) ÷ 4 + 2
    col = (i-1) % 4 + 1
    ax = Axis(fig[row, col], 
                title = "Channel $i ($(round(center_freqs[i]/1e9, digits=3)) GHz)",
                xlabel = "Rel. Freq (MHz)",
                ylabel = "dB",
                limits = (-70, 70, -80, 5))
    push!(axs_channels, ax)
end

# Colors for channels
all_colors = Makie.wong_colors()
colors = [all_colors[mod1(i, length(all_colors))] for i in 1:num_channels]

# Observables for animation
obs_input_freqs = Observable(collect(range(0, 4, length=1024)))
obs_input_mag = Observable(zeros(1024) .- 100)
obs_ch_freqs = [Observable(collect(range(-64, 64, length=128))) for _ in 1:num_channels]
obs_ch_mag = [Observable(zeros(128) .- 80) for _ in 1:num_channels]

# Plot input spectrum
lines!(ax_main, obs_input_freqs, obs_input_mag, color = :black, linewidth = 1.0, label = "Input Sweep")

# Plot channel filter responses
# Compute prototype filter response once
taps_per_branch = length(channelizer.polyphase_filters[1][1])
L_proto = taps_per_branch * M
h_proto = FreqChannelizer.design_prototype_filter(channel_bw, L_proto, fs)
H_f = fftshift(fft(h_proto, 4096))
H_mag = 20log10.(abs.(H_f) ./ maximum(abs.(H_f)) .+ 1e-12)
f_rel = range(-fs/2, fs/2, length=4096)

for i in 1:num_channels
    # Shifted prototype response for Visualization
    h_freq = (f_rel .+ center_freqs[i]) ./ 1e9
    lines!(ax_main, h_freq, H_mag, color = (colors[i], 0.7), linewidth = 2)
    
    # Channel subplots
    lines!(axs_channels[i], obs_ch_freqs[i], obs_ch_mag[i], color = colors[i], linewidth = 1.5)
end

axislegend(ax_main, position=:rt, nbanks=2, labelsize=10)

# Animation logic
frames = Int(duration / chunk_duration)
t_global = 0.0

println("Starting Animation Recording...")
record(fig, "sweep.mp4", 1:frames; framerate = 20) do frame_idx
    global t_global
    
    # Generate chunk of signal (chirp)
    # f(t) = f_start + (f_end - f_start) * t / duration
    # phi(t) = f_start*t + 0.5 * (f_end - f_start)/duration * t^2
    num_samples = Int(chunk_duration * fs)
    t_chunk = range(t_global, t_global + chunk_duration, length = num_samples)
    
    phi = 2π .* (f_start .* t_chunk .+ 0.5 .* (f_end - f_start) / duration .* t_chunk .^ 2)
    chunk = exp.(im .* phi)
    
    # Channelize
    output = channelize(channelizer, chunk)
    
    # Update Input Spectrum Plot
    p = periodogram(chunk, fs = fs)
    obs_input_freqs[] = p.freq ./ 1e9
    obs_input_mag[] = 10log10.(p.power .+ 1e-15)
    
    # Update Channel Plots
    out_fs = fs / M
    for i in 1:num_channels
        ch_p = periodogram(output[i, :], fs = out_fs)
        obs_ch_freqs[i][] = ch_p.freq ./ 1e6
        obs_ch_mag[i][] = 10log10.(ch_p.power .+ 1e-15)
    end
    
    t_global += chunk_duration
    if frame_idx % 20 == 0
        println("Processed frame $frame_idx / $frames")
    end
end

println("Animation saved to sweep.mp4")
