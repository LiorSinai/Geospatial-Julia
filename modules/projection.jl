import Base: inv, show
using StaticArrays
using Revise

"""
    AbstractProjection

An abstract projection. It should implement the following methods:

```
project(proj, coordinate::Tuple; options...) → Tuple
inv(proj) → AbstractProjection
```
The `coordinate` is a tuple of `(longitude, latitude)` and the outputs are `(x, y)` pairs.
"""
abstract type AbstractProjection end

function (proj::AbstractProjection)(xy::Tuple{T, T}; options...) where {T <: AbstractFloat}
    project(proj, xy; options...)
end

function (proj::AbstractProjection)(xy::SVector{2, T}; options...) where {T <: AbstractFloat}
    project(proj, (xy[1], xy[2]); options...)
end

degree_to_radian(x::T) where {T<:AbstractFloat} = x * π / 180
degree_to_radian(::Type{T}, x::AbstractFloat) where {T<:AbstractFloat} = convert(T, x * π / 180)
radian_to_degree(x::T) where {T<:AbstractFloat} = x * 180 / π
radian_to_degree(::Type{T}, x::AbstractFloat) where {T<:AbstractFloat} = convert(T, x * 180 / π)

"""
    earth_radius_at_latitude(latitude=35.196f0)

Earth radius at latitude (in degrees) based on the WGS-84 ellipsoid.
The default is set so that it is close to the mean radius: 6,371,008.5 m.

The radius ``r`` and latitude ``ϕ`` satisfy the following equations:
```math
 (x/a)^2 + (y/b)^2 = 1 
    x = r*cos(ϕ)
    y = r*sin(ϕ)
```
"""
function earth_radius_at_latitude(latitude::AbstractFloat=35.196f0)
    a = SEMI_MAJOR_AXIS_WGS_84
    b = SEMI_MINOR_AXIS_WGS_84
    latitude = latitude * π / 180
    radius = a * b / sqrt(b^2 * cos(latitude)^2 + a^2 * sin(latitude)^2)
    radius
end

### constants
SEMI_MAJOR_AXIS_WGS_84 = 6_378_137.0f0
INV_FLATTENING_WGS_84 = 298.257223563f0
SEMI_MINOR_AXIS_WGS_84 = SEMI_MAJOR_AXIS_WGS_84 * (1 - 1/INV_FLATTENING_WGS_84)
## The oblate ellipsoid has 2 of the same major semi-axis (x and y radii) and a slightly flatter minor axis (z radius).
MEAN_RADIUS_WGS_84 = (2 * SEMI_MAJOR_AXIS_WGS_84 + SEMI_MINOR_AXIS_WGS_84) / 3
AUTHALIC_RADIUS_WGS_84 = 6_371_007.2f0 # radius of sphere with same surface area as ellipsoid

include("interpolation.jl")
include("projections\\WGS_84.jl")
include("projections\\cylindrical_equal_area.jl")
include("projections\\mercator.jl")
include("projections\\robinson.jl")
include("projections\\transverse_mercator.jl")