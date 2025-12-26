using FreqChannelizer
using Test

@testset "FreqChannelizer.jl" begin
    @testset "Signal Reconstruction" begin
        # Parameters
        fs = 10e6
        duration = 0.002
        channel_bw = 200e3
        M = 8
        taps_per_branch = 32
        
        # Arbitrary center frequencies
        center_freqs = [
            0.5e6,
            1.2e6,
            2.1e6,
            2.8e6,
            3.5e6,
            4.3e6,
            -1.5e6,
            -2.8e6,
        ]
        
        channelizer = ArbitraryFrequencyChannelizer(center_freqs, fs, channel_bw, M, taps_per_branch)
        
        # Verify structure properties
        @test channelizer.M == M
        @test channelizer.num_channels == length(center_freqs)
        @test length(channelizer.polyphase_filters) == length(center_freqs)
        @test length(channelizer.polyphase_filters[1]) == M
        @test length(channelizer.polyphase_filters[1][1]) == taps_per_branch
        
        # Test tones corresponding to specific channels
        # Channel 1: 0.5 MHz
        # Channel 3: 2.1 MHz
        # Channel 5: 3.5 MHz
        # Channel 7: -1.5 MHz
        test_freqs = [0.5e6, 2.1e6, 3.5e6, -1.5e6]
        test_amps = [1.0, 0.8, 0.9, 0.7]
        
        input_signal = generate_multitone_signal(fs, duration, test_freqs, test_amps)
        
        output = channelize(channelizer, input_signal)
        
        @test size(output, 1) == length(center_freqs)
        # N_out = div(length(input_signal), M)
        @test size(output, 2) == div(length(input_signal), M)
        
        # Calculate channel power
        channel_power = vec(sum(abs2.(output), dims=2) ./ size(output, 2))
        channel_power_db = 10 .* log10.(channel_power ./ maximum(channel_power))
        
        println("Channel Power (dB): ", round.(channel_power_db, digits=1))
        
        # Active channels should be Ch 1, 3, 5, 7
        active_indices = [1, 3, 5, 7]
        inactive_indices = [2, 4, 6, 8]
        
        # Check active channels have high power (> -5 dB relative to max)
        for idx in active_indices
            # If input power is roughly uniform, max channel power should be close to 1
            # But due to scaling and passband ripple, we check relative db
            @test channel_power_db[idx] > -5.0
        end
        
        # Check inactive channels have low power (high rejection, e.g. < -40 dB)
        for idx in inactive_indices
            @test channel_power_db[idx] < -40.0
        end
    end
end
