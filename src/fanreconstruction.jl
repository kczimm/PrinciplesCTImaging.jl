fanprojection(e::Ellipse, βᵢ, γᵢ, D) = projection(e, βᵢ + γᵢ, D * sin(γᵢ))

function fanprojection(
    ellipses::AbstractVector{<:Ellipse},
    β::AbstractVector,
    γ::AbstractVector,
    D::AbstractFloat,
)
    [sum(fanprojection(e, βᵢ, γᵢ, D) for e in ellipses) for γᵢ in γ, βᵢ in β]
end

function g(N::Integer, α::AbstractFloat)
    # p. 83 eq. 96
    g = zeros(N)
    N2 = div(N, 2)
    for (i, gi) in enumerate(g)
        n = i - N2 - 1
        if mod(n, 2) != 0
            g[i] = 0.5 * (1 / (π * sin(n * α)))^2
        elseif n == 0
            g[i] = 1 / (8 * α^2)
        end
    end
    g
end

function fanWeight!(R::AbstractMatrix, γ::AbstractVector, D::AbstractFloat)
    views = size(R, 2)
    for βᵢ = 1:views
        @inbounds @. R[:, βᵢ] *= D * cos(γ)
    end
end

function fanFilter!(Q::AbstractMatrix, α::AbstractFloat)
    samples, views = size(Q)
    N = nextpow(2, 2 * samples - 1)
    ZP(x) = vcat(x, zeros(N - samples))
    w = ZP(g(samples, α))
    # indices to slice out unpadded projection?
    i = samples ÷ 2 + 1
    j = i + samples - 1
    for βᵢ = 1:views
        # p. 83 eq. 94
        @inbounds Q[:, βᵢ] .= α * real.(ifft(fft(ZP(Q[:, βᵢ])) .* fftshift(fft(w))))[i:j]
    end
end

function fanReconstruction(
    R::AbstractMatrix,
    β::AbstractVector,
    γ::AbstractVector,
    D::AbstractFloat,
)
    Q = copy(R)
    α = step(γ)
    fanWeight!(Q, γ, D)
    fanFilter!(Q, α)

    Δβ = 2π / length(β)
    pixels = 128
    f = zeros(pixels, pixels)

    Q = LinearInterpolation((γ, β), Q, extrapolation_bc = Line())

    for (i, xᵢ) in enumerate(range(-1, 1, length = pixels)),
        (j, yⱼ) in enumerate(range(-1, 1, length = pixels))

        @inbounds f[j, i] = sum(β) do βᵢ
            r = √(xᵢ^2 + yⱼ^2)
            ϕ = atan(yⱼ, xᵢ)
            γ′ = atan((r * cos(βᵢ - ϕ)) / (D + r * sin(βᵢ - ϕ)))
            L² = (xᵢ - D * cos(βᵢ))^2 + (yⱼ - D * sin(βᵢ))^2
            Q(γ′, βᵢ) / L²
        end
    end
    Δβ .* f
end
