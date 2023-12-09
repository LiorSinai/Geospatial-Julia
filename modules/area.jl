using PolygonAlgorithms: area_polygon, bounds
using Rasters
using Rasters: crop

function area_geometry(polygon::GeoJSON.Polygon)
    area = 0.0
    for region in polygon.coordinates
        area += area_polygon(region)
    end
    area
end

function area_geometry(multipolygon::GeoJSON.MultiPolygon)
    area = 0.0
    for polygon in multipolygon.coordinates
        for region in polygon
            area += area_polygon(region)
        end
    end
    area
end

function calc_density_per_area(
    population_data::Raster,
    features::Vector{GeoJSON.Feature},
    src_proj
    )
    n = length(features)
    areas = Vector{Float64}(undef, n)
    population_counts = Vector{Float64}(undef, n)
    population_densities = Vector{Float64}(undef, n)
    equal_area_proj = CylindricalEqualArea()
    for (idx, feature) in enumerate(features)
        region = crop(population_data; to=feature.geometry)
        region = mask(region; with=feature.geometry) # slow without cropping first
        population_count = sum(x -> x > 0 ? x : 0, region)
        projected = reproject(feature, src_proj, equal_area_proj)
        area = area_geometry(projected.geometry)
        areas[idx] = area
        population_counts[idx] = population_count
        population_densities[idx] = population_count / area
    end
    areas, population_counts, population_densities
end
