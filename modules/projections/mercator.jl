"""
    Mercator(radius, long0, k)
    Mercator(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

Conformal projection from a sphere which preserves shape but creates large distortion near the poles. 

It is used for navigation where conformality is of the utmost importance.
It should not be used for purely visual applications because of the large distortions.
Instead it is strongly recommended to use the `Robinson` projection for global maps and the `TransverseMercator` projection for small local regions.

The equations satisfy:
```math
    x = k*r*(λ-λ_0)
    y = k*r*log(tan(π/4 + ϕ/2))
```

Inputs should be in degrees:
```
proj = Mercator()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (2.22e6, 3.50e6) #m
```
"""
struct Mercator{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    k::T ## scale
end

function show(io::IO, mime::MIME"text/plain", proj::Mercator)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "k=",proj.k)
    print(io, ")")
end

function Mercator(;long0::AbstractFloat=0.0f0, k::AbstractFloat=1.0f0, radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    Mercator(promote(radius, long0, k)...)
end

inv(proj::Mercator) = InverseMercator(proj.radius, proj.long0, proj.k)

function project(proj::Mercator{T1}, coordinate::Tuple{T2, T2}; atol=1e-6) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude -= proj.long0
    longitude = degree_to_radian(T1, longitude)
    latitude  = degree_to_radian(T1, latitude)
    x = proj.k * proj.radius * longitude
    angle = convert(T1, π/4)  + latitude/2
    at_pole = abs(abs(angle) - π/2) < atol
    y = at_pole ? T1(Inf) : proj.k * proj.radius * log(tan(angle))
    (x, y)
end

"""
    InverseMercator(radius, long0, k)

Convert `(x, y)` co-ordinates in a `Mercator` projection back to `(longitude, latitude)`. 
"""
struct InverseMercator{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    k::T ## scale
end

inv(proj::InverseMercator) = InverseMercator(proj.radius, proj.long0, proj.k)

function project(proj::InverseMercator{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    x, y = xy
    x = x / (proj.k * proj.radius)
    y = y / (proj.k * proj.radius)
    longitude = x
    latitude = 2 * atan(exp(y)) - convert(T1, π/2)
    longitude = radian_to_degree(T1, longitude)
    latitude = radian_to_degree(T1, latitude)
    longitude += proj.long0
    (longitude, latitude)
end