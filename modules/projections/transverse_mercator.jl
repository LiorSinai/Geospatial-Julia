"""
    TransverseMercator(radius, long0, k)
    TransverseMercator(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

Projection from a sphere  which corrects distortion along the longitude `long0` at the expense of large distortion away from it.

The equations satisfy:
```math
    t = sin(λ-λ_0)cos(ϕ)
    x = k*r*0.5*log((1+t) / (1-t))
    y = k*r*atan(tan(ϕ)*sec(λ-λ_0))
```

Inputs should be in degrees:
```
proj = TransverseMercator()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (1.95e6, 3.51e6) #m
```
"""
struct TransverseMercator{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    k::T ## scale
end

function show(io::IO, mime::MIME"text/plain", proj::TransverseMercator)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "k=",proj.k)
    print(io, ")")
end

function TransverseMercator(;long0::AbstractFloat=0.0f0, k::AbstractFloat=1.0f0, radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    TransverseMercator(promote(radius, long0, k)...)
end

inv(proj::TransverseMercator) = InverseTransverseMercator(proj.radius, proj.long0, proj.k)

function project(proj::TransverseMercator{T1}, coordinate::Tuple{T2, T2}; extend::Bool=false) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude -= proj.long0
    longitude = degree_to_radian(T1, longitude)
    latitude  = degree_to_radian(T1, latitude)
    t = sin(longitude) * cos(latitude)
    x = proj.k * proj.radius * convert(T1, 0.5) * log((one(T1) + t) / (one(T1) - t))
    y = proj.k * proj.radius * atan(tan(latitude) * sec(longitude))
    if extend && (abs(longitude) > π/2) # extend && rear of sphere
        y += proj.k * proj.radius * π
    end
    (x, y)
end

"""
    InverseTransverseMercator(radius, long0, k)

Convert `(x, y)` co-ordinates in a `SphericalTransverseMercator` projection back to `(longitude, latitude)`. 
"""
struct InverseTransverseMercator{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    k::T ## scale
end

inv(proj::InverseTransverseMercator) = TransverseMercator(proj.radius, proj.long0, proj.k)

function project(proj::InverseTransverseMercator{T1}, xy::Tuple{T2, T2}; extend=true) where {T1,T2 <: AbstractFloat}
    x, y = xy
    x = x / (proj.k * proj.radius)
    y = y / (proj.k * proj.radius)
    longitude = atan(sinh(x) * sec(y))
    latitude = asin(sech(x) * sin(y))
    longitude = radian_to_degree(T1, longitude)
    latitude = radian_to_degree(T1, latitude)
    longitude += proj.long0
    if abs(y) > π/2
        if longitude < 0
            longitude += 180
        else
            longitude -= 180
        end
    end
    (longitude, latitude)
end
