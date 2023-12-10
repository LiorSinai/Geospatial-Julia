"""
    WorldGeodeticSystem84()

Reference: https://earth-info.nga.mil/?dir=wgs84&action=wgs84.

An ellipsoid with semi major axis ``a`` and flattening ``f=(a-b)/a``.

This projection acts as an identity projection.

Example:
```
proj = WorldGeodeticSystem84()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (20.0f0, 30.0f0)
```
"""
struct WorldGeodeticSystem84 <: AbstractProjection
    semi_major_axis::Float32
    inverse_flattening::Float32
    WorldGeodeticSystem84() = new(SEMI_MAJOR_AXIS_WGS_84, INV_FLATTENING_WGS_84)
end

function show(io::IO, mime::MIME"text/plain", proj::WorldGeodeticSystem84)
    print(io, "WorldGeodeticSystem84(")
    print(io, "semi_major_axis=", proj.semi_major_axis, ", ")
    print(io, "flattening=1/", proj.inverse_flattening)
    print(io, ")")
end

project(proj::WorldGeodeticSystem84, coordinate::Tuple) = coordinate
inv(proj::WorldGeodeticSystem84) = proj
