function h(N::Integer, τ::AbstractFloat)
    # p. 72 eq. 61
    h = zeros(N)
    N2 = div(N, 2)
    for (i, hi) in enumerate(h)
        n = i - N2 - 1
        if mod(n, 2) != 0
            h[i] = -1 / (n * π * τ)^2
        elseif n == 0
            h[i] = 1 / (4 * τ^2)
        end
    end
    h
end

function parallelFilter!(Q::AbstractMatrix, P::AbstractMatrix, τ::AbstractFloat)
    @assert size(Q) == size(P)
    samples, views = size(P)
    N = nextpow(2, 2 * samples - 1)
    ZP(x) = vcat(x, zeros(N - samples))
    w = ZP(h(samples, τ))
    # indices to slice out unpadded projection?
    i = div(samples, 2) + 1
    j = i + samples - 1
    for θ = 1:views
        # p. 74 eq. 68
        @inbounds Q[:, θ] = τ * real.(ifft(fft(ZP(P[:, θ])) .* fft(w)))[i:j]
    end
end

function parallelReconstruction(P::AbstractMatrix, θ::AbstractVector, t::AbstractVector)
    Q = similar(P)
    τ = step(t)
    parallelFilter!(Q, P, τ)

    K = length(θ)
    pixels = 128
    f = zeros(pixels, pixels)

    Q = LinearInterpolation((t, θ), Q, extrapolation_bc = Line())

    for (i, xᵢ) in enumerate(range(-1, 1, length = pixels)),
        (j, yᵢ) in enumerate(range(-1, 1, length = pixels))

        @inbounds f[j, i] = sum(Q(xᵢ * cos(θᵢ) + yᵢ * sin(θᵢ), θᵢ) for θᵢ in θ)
    end
    @. f * π / K
end
