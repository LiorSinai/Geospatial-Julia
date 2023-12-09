import Base: show, similar

function linear_interpolation(x1::AbstractFloat, y1::AbstractFloat, x2::AbstractFloat, y2::AbstractFloat, x::AbstractFloat)
    grad = (y2 - y1) / (x2 - x1)
    inter = (y1 * x2 - y2 * x1) / (x2 - x1)
    grad * x + inter
end

function _get_interval_idx(
    nodes::AbstractVector, x
    )
    idx = 2
    # assume nodes are ordered.
    # could do this with binary search in log(n) time.
    while (idx) < length(nodes) && (x > nodes[idx])
        idx += 1
    end
    idx # idx of first value greater than x
end

"""
    LinearInterpolater(intervals, values)

Linear interpolation between intervals.

```
    y = (y1*(x2-x)-y2*(x1-x))/(x2-x1) where xs[idx] < x < xs[idx+1]
```

Assumes `intervals` are ordered.
"""
struct LinearInterpolater{V1<:AbstractVector, V2<:AbstractVector}
    intervals::V1
    values::V2
end

LinearInterpolater() = LinearInterpolater(Float64[], Float64[])
similar(::LinearInterpolater, xs, ys) = LinearInterpolater(xs, ys)

function interpolate(interpolator::LinearInterpolater, x)
    idx = _get_interval_idx(interpolator.intervals, x)
    xs = interpolator.intervals
    ys = interpolator.values
    linear_interpolation(xs[idx - 1], ys[idx - 1], xs[idx], ys[idx], x)
end

(interpolator::LinearInterpolater)(x) = interpolate(interpolator, x)

function show(io::IO, mime::MIME"text/plain", interpolator::LinearInterpolater)
    print(io, "LinearInterpolater(")
    print(io, interpolator.intervals, ", ")
    print(io, interpolator.values)
    print(io, ")")
end

"""
    lagrange_polynomial(nodes, values, x)

Evaluates the Lagrange polynomial through (node, value) pairs at point `x`.

The Lagrange polynomial is the unique polynomial of lowest degree that intersects all data points.

The values near the boundaries are not guaranteed to be smooth.

Sources:
- https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
- https://en.wikipedia.org/wiki/Lagrange_polynomial
"""
function lagrange_polynomial(
    nodes::AbstractVector{<:AbstractFloat}, values::AbstractVector{<:AbstractFloat}, x::AbstractFloat
    )
    n = length(nodes)
    y = 0.0
    for (j, yj) in enumerate(values)
        ℓ = 1.0
        for k in [1:(j-1) ; (j+1):n]
            ℓ *= (x - nodes[k]) / (nodes[j] - nodes[k])
        end
        y += yj * ℓ
    end
    y
end

"""
    nevilles_algorithm(nodes, values, x)

Evaluates the Lagrange polynomial at point `x`.

The values near the boundaries are not guaranteed to be smooth.

Sources:
- https://mathworld.wolfram.com/NevillesAlgorithm.html
- https://en.wikipedia.org/wiki/Neville%27s_algorithm
"""
function nevilles_algorithm(
    nodes::AbstractVector{<:AbstractFloat}, values::AbstractVector{<:AbstractFloat}, x::AbstractFloat
    )
    n = length(nodes)
    _nevilles_algorithm(nodes, values, x, 1, n)
end

function _nevilles_algorithm(
    nodes::AbstractVector{<:AbstractFloat}, values::AbstractVector{<:AbstractFloat}, x::AbstractFloat, i::Int, j::Int
    )
    if (i == j)
        return values[i]
    end
    p1 = _nevilles_algorithm(nodes, values, x, i + 1, j)
    p2 = _nevilles_algorithm(nodes, values, x, i, j - 1)
    ((x - nodes[i]) * p1 - (x - nodes[j]) * p2) / (nodes[j] - nodes[i])
end

"""
    CubicSpline(xs, ys)
    CubicSpline(coefficients, )

The coefficients for cubic spline interpolation. Boundary condition are that the second derivatives are zero.
    
The coefficients are stored as a matrix. Each column represents a polynomial ``j`` where the coefficient of ``x^i`` is in row ``i``.

That is:
```
    y = sum(A[i, j] * (x - x[j])^i) where xs[j] < x < xs[j+1]
```

Assumes `xs` are ordered.
"""
struct CubicSpline{V<:AbstractVector}
    coefficients::Matrix
    intervals::V
end

function CubicSpline(xs::AbstractVector{T}, ys::AbstractVector{T}) where T
    @assert length(xs) == length(ys)
    n = length(xs)
    hs = xs[2:n] - xs[1:(n - 1)]
    A = zeros(T, n, n)
    b = zeros(T, n)
    A[1, 1] = 1
    A[n, n] = 1
    for i in 1:(n-2)
        A[i + 1, i] = hs[i]
        A[i + 1, i + 1] = 2 * (hs[i] + hs[i + 1])
        A[i + 1, i + 2] = hs[i + 1]
        b[i + 1] = 3 * (ys[i + 2] - ys[i + 1]) / hs[i + 1] - 3 * (ys[i + 1] - ys[i])  / hs[i]
    end
    c = A \ b
    coefficients = zeros(T, 4, n - 1)
    col = 1:(n-1)
    col_next = 2:n 
    coefficients[1, :] = ys[col]
    coefficients[2, :] = (ys[col_next] - ys[col]) ./ hs[col] - hs[col] / 3 .* (2 * c[col] + c[col_next])
    coefficients[3, :] = c[col]
    coefficients[4, :] = (c[col_next] - c[col]) ./ (3 * hs[col])
    CubicSpline{typeof(xs)}(coefficients, xs)
end

CubicSpline() = CubicSpline(zeros(0, 0), zeros(0))
similar(::CubicSpline, xs, ys) = CubicSpline(xs, ys)

function interpolate(spline::CubicSpline, x)
    j = _get_interval_idx(spline.intervals, x) - 1
    coeffs = spline.coefficients[:, j]
    xj = spline.intervals[j]
    polynomial(coeffs, x - xj)
end

(spline::CubicSpline)(x) = interpolate(spline, x)
 
function polynomial(coefficients::AbstractVector, x)
    result = 0.0
    xpow = 1.0
    for coeff in coefficients
        result += coeff * xpow
        xpow *= x
    end
    result
end

function show(io::IO, mime::MIME"text/plain", spline::CubicSpline)
    print(io, "CubicSpline")
    print(io, "\n  intervals = ", spline.intervals)
    print(io, "\n  coefficients = ")
    show(io, mime, spline.coefficients)
end
