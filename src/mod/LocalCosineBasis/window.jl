"""
Construct a bell of length window_length for local cosine analysis
"""
struct OrthonormalBell
    interior::Vector{Float64}
    exterior::Vector{Float64}
    function OrthonormalBell(m::Int)
        interior = Vector{Float64}(undef,m)
        exterior = Vector{Float64}(undef,m)
        map!(t -> sin(pi/4 * (1 + t/m)),interior,.5:m-.5)
        map!(x-> sqrt(1 - x^2),exterior, interior)
        new(interior,exterior)
    end
end
# Base.getindex(o::OrthonormalBall,idx::Int) = (o.interior[idx],o.exterior[idx])
Base.size(o::OrthonormalBell,idx::Int) = size(o.interior)
Base.length(o::OrthonormalBell) = length(o.interior)
interior(o::OrthonormalBell) = o.interior
exterior(o::OrthonormalBell) = o.exterior
