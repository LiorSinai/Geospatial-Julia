"""
    CylindricalEqualArea(radius, long0, k)
    CylindricalEqualArea(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

Projection from a sphere which preserves area.
It has large distortion near the poles.

The equations satisfy:
```math
    x = k*r*(λ-λ_0)
    y = k*r*sin(ϕ)
```

Inputs should be in degrees:
```
proj = CylindricalEqualArea()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (2.22e6, 3.19e6) m
```
"""
struct CylindricalEqualArea{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    k::T ## scale
end

inv(proj::CylindricalEqualArea) = InverseCylindricalEqualArea(
    proj.radius, proj.long0, proj.k)

function show(io::IO, mime::MIME"text/plain", proj::CylindricalEqualArea)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "k=",proj.k)
    print(io, ")")
end

function CylindricalEqualArea(;long0::AbstractFloat=0.0f0, k::AbstractFloat=1.0f0, radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    CylindricalEqualArea(promote(radius, long0, k)...)
end

function project(proj::CylindricalEqualArea{T1}, coordinate::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude -= proj.long0
    longitude = degree_to_radian(T1, longitude)
    latitude  = degree_to_radian(T1, latitude)
    x = proj.k * proj.radius * longitude
    y = proj.k * proj.radius * sin(latitude)
    (x, y)
end

"""
    InverseCylindricalEqualArea(radius, long0, k)

Convert `(x, y)` co-ordinates in a `CylindricalEqualArea` projection to `(longitude, latitude)`. 
"""
struct InverseCylindricalEqualArea{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    k::T ## scale
end

inv(proj::InverseCylindricalEqualArea) = CylindricalEqualArea(proj.radius, proj.long0, proj.k)

function show(io::IO, mime::MIME"text/plain", proj::InverseCylindricalEqualArea)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "k=", proj.k)
    print(io, ")")
end

function project(proj::InverseCylindricalEqualArea{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    x, y = xy
    x = x / (proj.k * proj.radius)
    y = y / (proj.k * proj.radius)
    longitude = x
    latitude = asin(y)
    longitude = radian_to_degree(T1, longitude)
    latitude = radian_to_degree(T1, latitude)
    longitude += proj.long0
    (longitude, latitude)
end