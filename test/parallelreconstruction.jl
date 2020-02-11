using PrinciplesCTImaging: h

@testset "impulse response" begin
    N = 13
    τ = 0.1
    hn = h(N, τ)
    @test hn[7] ≈ 1/(4*τ^2)
end
