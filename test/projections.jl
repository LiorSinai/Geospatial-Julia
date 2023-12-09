using Test

@testset "WGS-84" begin
    proj = WorldGeodeticSystem84()
    coords_src = (20.0f0, 30.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test coords_dest == coords_src

    proj_inv = inv(proj)
    @test proj_inv == proj
end

@testset "cylindrical equal area - r=1" begin
    proj = CylindricalEqualArea(1.0f0, 0.0f0, 1.0f0)
    coords_src = (20.0f0, 30.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test all(coords_dest .≈ (0.34906584f0, 0.5f0))

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest) 
    @test coords_back isa Tuple{Float32, Float32}
    @test all(coords_back .≈ (20.0, 30.0))
end

@testset "cylindrical equal area" begin
    proj = CylindricalEqualArea()
    coords_src = (20.0, 30.0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test all(coords_dest .≈ (2.2239015f6, 3.1855042f6))

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest) 
    @test coords_back isa Tuple{Float32, Float32}
    @test all(coords_back .≈ (20.0, 30.0))
end

@testset "Mercator" begin
    proj = Mercator(;long0=0.0f0)
    coords_src = (20.0f0, 30.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test all(coords_dest .≈ (2.2239015f6, 3.4996345f6))

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest)
    @test coords_back isa Tuple{Float32, Float32}
    @test all(coords_back .≈ (20.0, 30.00000333930422))
end

@testset "Mercator - pole" begin
    proj = Mercator(1.0f0, 0.0f0, 1.0f0)
    coords_src = (0.0f0, 90.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test all(coords_dest .== (0.0, Inf32))

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest) 
    @test coords_back isa Tuple{Float32, Float32}
    @test all(coords_back .≈ (0.0, 90.00000500895634))
end

@testset "Robinson" begin
    proj = Robinson(1.0f0, 0.0f0, 1.0f0, CubicSpline())
    coords_src = (20.0f0, 30.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest == (0.28440211575222013, 0.5030556)

    proj = Robinson(1.0, 0.0, 1.0, LinearInterpolater())
    coords_dest = proj(coords_src)
    @test coords_dest ==  (0.28440209974417663, 0.5030556)

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest) 
    @test coords_back ==  coords_src
end


@testset "transverse Mercator - r=1" begin
    proj = TransverseMercator(1.0f0, 0.0f0, 1.0f0)
    coords_src = (20.0f0, 30.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test all(coords_dest .≈ (0.30534694f0, 0.5509429f0))

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest) 
    @test coords_back isa Tuple{Float32, Float32}
    @test all(coords_back .≈ (20.0, 30.0))
end

@testset "transverse Mercator - pole" begin
    proj = TransverseMercator(1.0f0, 0.0f0, 1.0f0)
    coords_src = (90.0f0, 0.0f0) # degrees
    coords_dest = proj(coords_src)
    @test coords_dest isa Tuple{Float32, Float32}
    @test all(coords_dest .== (Inf, 0.0))

    proj_inv = inv(proj)
    coords_back = proj_inv(coords_dest) 
    @test coords_back isa Tuple{Float32, Float32}
    @test all(coords_back .≈ (90.0, 0.0))
end

@testset "transverse Mercator - extended" begin
    proj = TransverseMercator(1.0f0, 0.0f0, 1.0f0)
    coords_src = (150.0f0, 15.0f0) # degrees
    coords_dest = project(proj, coords_src; extend=true)
    @test all(coords_dest .== (0.5268411f0, 2.8415337f0))

    proj_inv = inv(proj)
    coords_back = project(proj_inv, coords_dest) 
    @test all(coords_back .≈ (150.0f0, 14.999991f0)) 

    coords_src = (-150.0f0, 15.0f0) # degrees
    coords_dest = project(proj, coords_src; extend=true)
    @test all(coords_dest .== (-0.5268411f0, 2.8415337f0))

    proj_inv = inv(proj)
    coords_back = project(proj_inv, coords_dest) 
    @test all(coords_back .≈ (-150.0f0, 15.0f0))
end