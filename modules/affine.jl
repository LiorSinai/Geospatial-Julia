using StaticArrays
import Base: inv, show
using Rasters

"""
    AffineTransform(matrix)

Maps ``(x', y')`` to ``(λ, ϕ)`` where λ is the longitude and ϕ is the latitude. 

Note that for images, one needs to map ``(j, i)`` to ``(x', y')``.
"""
struct AffineTransform{T}
    matrix::Matrix{T}
end

inv(a::AffineTransform) = AffineTransform(inv(a.matrix))

function (a::AffineTransform)(coordinate::Tuple{T, T}) where T
    x, y = coordinate
    xnew = a.matrix[1, 1] * x + a.matrix[1, 2] * y + a.matrix[1, 3]
    ynew = a.matrix[2, 1] * x + a.matrix[2, 2] * y + a.matrix[2, 3]
    (xnew, ynew)
end

function (a::AffineTransform)(coordinate::SVector{2, T}) where T
    x, y = coordinate
    xnew = a.matrix[1, 1] * x + a.matrix[1, 2] * y + a.matrix[1, 3]
    ynew = a.matrix[2, 1] * x + a.matrix[2, 2] * y + a.matrix[2, 3]
    (xnew, ynew)
end

function affine_translation(Δx::T, Δy::T) where T
    A = [
        1 0 Δx;
        0 1 Δy;
        0 0 1
    ]
    AffineTransform{T}(A)
end

function affine_scale(sx::T, sy::T) where T
    A = [
        sx 0  0;
        0  sy 0;
        0  0  1;
    ]
    AffineTransform{T}(A)
end

function show(io::IO, mime::MIME"text/plain", a::AffineTransform)
    println(io, "AffineTransform with matrix:")
    show(io, mime, a.matrix)
end

"""
    affine_from_bounds(west, south, east, north, width, height)

Maps (1, 1) to (west, north) and (width, height) to (east, south).
"""
function affine_from_bounds(
    west::T, south::T, east::T, north::T, width::Int, height::Int
    ) where T <: AbstractFloat
    sx = (east - west) / (width - 1)
    sy = (south - north) / (height - 1)
    A1 = affine_translation(west - sx, north - sy)
    A2 = affine_scale(sx, sy)
    AffineTransform{T}(A1.matrix * A2.matrix)
end

function affine_from_raster(raster::Raster)
    long_min, long_max = extrema(raster.dims[1])
    lat_min, lat_max = extrema(raster.dims[2])
    width, height = size(raster)
    affine_from_bounds(long_min, lat_min, long_max, lat_max, width, height)
end
