abstract type Shape end

struct Ellipse{T} <: Shape
    x::T
    y::T
    A::T
    B::T
    α::T
    ρ::T
end

function isinside(e::Ellipse, x, y)
    # https://stackoverflow.com/questions/7946187/point-and-ellipse-rotated-position-test-algorithm
    (cosd(e.α) * (x - e.x) + sind(e.α) * (y - e.y))^2 / e.A^2 +
    (sind(e.α) * (x - e.x) - cosd(e.α) * (y - e.y))^2 / e.B^2 ≤ one(x)
end

function projection(e::Ellipse, θᵢ, tᵢ)
    # p. 56
    s = √(e.x^2 + e.y^2)
    γ = atan(e.y, e.x)
    θ = θᵢ - deg2rad(e.α)
    t = tᵢ - s * cos(γ - θᵢ)
    a²θ = e.A^2 * cos(θ)^2 + e.B^2 * sin(θ)^2
    abs(t) ≤ √(a²θ) ? 2 * e.ρ * e.A * e.B * √(a²θ - t^2) / a²θ : zero(a²θ)
end

function projection(
    ellipses::AbstractVector{<:Ellipse},
    θ::AbstractVector,
    t::AbstractVector,
)
    [sum(projection(e, θᵢ, tᵢ) for e in ellipses) for tᵢ in t, θᵢ in θ]
end

# Table 3.1 summary of parameters for tomography simulations
SheppLoganPhantom = [
    Ellipse(0.0, 0.0, 0.92, 0.69, 90.0, 2.0),
    Ellipse(0.0, -0.0184, 0.874, 0.6624, 90.0, -0.98),
    Ellipse(0.22, 0.0, 0.31, 0.11, 72.0, -0.02),
    Ellipse(-0.22, 0.0, 0.41, 0.16, 108.0, -0.02),
    Ellipse(0.0, 0.35, 0.25, 0.21, 90.0, 0.01),
    Ellipse(0.0, 0.1, 0.046, 0.046, 0.0, 0.01),
    Ellipse(0.0, -0.1, 0.046, 0.046, 0.0, 0.01),
    Ellipse(-0.08, -0.605, 0.046, 0.023, 0.0, 0.01),
    Ellipse(0.0, -0.605, 0.023, 0.023, 0.0, 0.01),
    Ellipse(0.06, -0.605, 0.046, 0.023, 90.0, 0.01),
]

function construct(ellipses::AbstractVector{<:Ellipse}, N::Integer)
    f = zeros(N, N)
    for (i, xᵢ) in enumerate(range(-1, 1, length = N)),
        (j, yⱼ) in enumerate(range(-1, 1, length = N))

        @inbounds f[j, i] = sum([isinside(e, xᵢ, yⱼ) ? e.ρ : zero(e.ρ) for e in ellipses])
    end
    f
end
