using Test

@testset "linear interpolation" begin
    x1, y1, x2, y2 = 0.6, 10.0, 0.8, 15.5

    # middle
    y_interp = linear_interpolation(x1, y1, x2, y2, 0.7)
    @test y_interp == 12.75

    # boundary
    y_interp = linear_interpolation(x1, y1, x2, y2, x1)
    @test y_interp == y1
    y_interp = linear_interpolation(x1, y1, x2, y2, x2)
    @test y_interp == y2

    # extrapolation
    y_interp = linear_interpolation(x1, y1, x2, y2, 1.5)
    @test y_interp ≈ 34.75
end

@testset "LinearInterpolater" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = [5.0, 7.0, 11.0, 34.0]

    interpolater = LinearInterpolater(nodes, values)

    y_interp = interpolater(0.0)
    @test y_interp == (6 + 1/3)
end

@testset "Lagrange polynomial - degree 3" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = [5.0, 7.0, 11.0, 34.0]

    ys = [lagrange_polynomial(nodes, values, x) for x in nodes]
    @test ys == values

    y_interp = lagrange_polynomial(nodes, values, 0.0)
    @test y_interp ≈ 1087/180
end

@testset "Lagrange polynomial - known" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = nodes.^3 + 2 * (nodes.^2) .+ 2.3

    ys = [lagrange_polynomial(nodes, values, x) for x in nodes]
    @test ys == values

    y_interp = lagrange_polynomial(nodes, values, 1.5)
    expected = 1.5^3 + 2 * 1.5^2 + 2.3
    @test y_interp ≈ expected
end

@testset "Neville's algorithm - degree 3" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = [5.0, 7.0, 11.0, 34.0]

    ys = [nevilles_algorithm(nodes, values, x) for x in nodes]
    @test ys == values

    y_interp = nevilles_algorithm(nodes, values, 0.0)
    @test y_interp ≈ 1087/180
end

@testset "Neville's algorithm - known" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = nodes.^3 + 2 * (nodes.^2) .+ 2.3

    ys = [nevilles_algorithm(nodes, values, x) for x in nodes]
    @test ys ≈ values

    y_interp = nevilles_algorithm(nodes, values, 1.5)
    expected = 1.5^3 + 2 * 1.5^2 + 2.3
    @test y_interp ≈ expected
end

@testset "CubicSpline - coefficients" begin
    nodes = [0.0, 1.0, 2.0, 2.5]
    values = [0.0, 1.0, 8.0, 9.0]

    spline = CubicSpline(nodes, values)
    expected = [
        0    11   88
        -12  57   48
        0    69   -78
        23   -49  52
    ]
    @test spline.coefficients * 11 ≈ expected 

    ys = [spline(x) for x in nodes]
    @test ys ≈ values
end

@testset "CubicSpline" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = [5.0, 7.0, 11.0, 34.0]
    spline = CubicSpline(nodes, values)

    spline_inv = CubicSpline(values, nodes)

    ys = [spline(x) for x in nodes]
    @test ys == values

    y_interp = spline(0.0)
    @test y_interp ≈ 6.0890804597701145
end
