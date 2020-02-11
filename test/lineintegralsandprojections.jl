using PrinciplesCTImaging: Ellipse, isinside

@testset "points inside ellipse" begin
    e = Ellipse(0.0, 0.0, 1.0, 1.0, 0.0, 1.0)
    @test isinside(e, 0.0, 0.0) == true
    @test isinside(e, 1.0, 0.0) == true
    @test isinside(e, 1.0, 1.0) == false
end

using PrinciplesCTImaging: SheppLoganPhantom, construct

@testset "SheppLogan phantom construction" begin
    pixels = 128
    I = construct(SheppLoganPhantom, pixels)
    @test I[1, 1] ≈ 0.0
    @test I[25, 50] ≈ 1.02
    @test I[50, 50] ≈ 1.0
end
