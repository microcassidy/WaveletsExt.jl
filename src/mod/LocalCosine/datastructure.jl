abstract type AbstractPacketPosition end
abstract type LeftPacket <: AbstractPacketPosition end
abstract type RightPacket <: AbstractPacketPosition end
abstract type CentrePacket <: AbstractPacketPosition  end

mutable struct VariableVector{Float64} <: AbstractVector{Float64}
    data::Vector{Float64}
    length::Int64
    max_length::Int64
    VariableVector(m::Int64) = new{Float64}(zeros(Float64,m),m,m)
end
_update_size!(v::VariableVector{Float64},m::Int64) = 1<=m<=v.max_length ? setfield!(v,:length,m) : error("outside bounds of 1 and $(v.max_length)")
Base.length(v::VariableVector{Float64}) = v.length
Base.size(v::VariableVector{Float64}) = (v.length,)
Base.getindex(v::VariableVector{Float64},i::Int)::Float64 = v.data[i]
function Base.setindex!(v::VariableVector{Float64},val::Float64,ind::Int)
    v.data[ind] = val
    return nothing
end
Base.show(v::VariableVector) = show(v.data)

Base.firstindex(v::VariableVector) = 1
Base.lastindex(v::VariableVector) = v.length

"""
the left,right,centre of a window of data.
A packet is constructed at max size (level 0 decomposition).
This is to avoid reallocation of the array during analysis/synthesis steps.
The size remains constant but the displayed amount will shift. See VariableVector.

The left left/right packet are not varied. The display will be the last m terms of the left and the first m terms of the
right, respectively. This corresponds to the external parts of the length 2m window skirt
"""
mutable struct Packet
    left::Vector{Float64}
    right::Vector{Float64}
    centre::VariableVector{Float64}
    m::Int64
    Packet(max_size::Int64,m::Int64) = new(Vector{Float64}(undef,m), Vector{Float64}(undef,m), VariableVector(max_size),m)
end

function reset!(a::VariableVector{Float64})::Nothing
    for idx=eachindex(a)
        a[idx] = zero(eltype(a))
    end
    return nothing
end
function reset!(a::Vector{Float64})::Nothing
    for idx=eachindex(a)
        a[idx] = zero(eltype(a))
    end
    nothing
end
# reset!(p::Packet,::Type{PacketLeft})

function reset!(p::Packet)::Nothing
    reset!(p.left)
    reset!(p.right)
    reset!(p.centre)
    nothing
end

#the only packet that needs to grow or shrink is actually the centre.

function set!(::Type{LeftPacket},p::Packet, v::AbstractVector{Float64})
    isempty(v) && return
    l = length(v)
    rng = l-p.m+1:l
    !isempty(rng) && copyto!(p.left, v[rng]) #back of v
    return nothing
end
set!(::Type{RightPacket},p::Packet, v::AbstractVector{Float64}) = !isempty(v[begin:p.m]) && copyto!(p.right, v[begin:p.m]) #front of v
set!(::Type{CentrePacket},p::Packet, v::AbstractVector{Float64}) = copyto!(p.centre,v)
_update_size!(p::Packet,m::Int) = _update_size!(p.centre,m)

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
        new(interior,reverse(exterior))
    end
end
# Base.getindex(o::OrthonormalBell,idx::Int) = (o.interior[idx],o.exterior[idx])
Base.size(o::OrthonormalBell,idx::Int) = size(o.interior)
Base.length(o::OrthonormalBell) = length(o.interior)
interior(o::OrthonormalBell) = o.interior
exterior(o::OrthonormalBell) = o.exterior

mutable struct LocalCosineBasis
    bell::OrthonormalBell
    dct_plans::Vector{DCT_T}
    packet::Packet
    J::Int64
    m::Int64
    max_depth::Int
    N::Int

    function LocalCosineBasis(vector_length::Int = DEFAULT_BASIS_TOPLEVEL_LENGTH;max_depth::Int=4)
        J = dyadic_length(vector_length)
        m = div(vector_length,1 << (max_depth + 1)) #half skirt filter size
        bell = OrthonormalBell(m)
        packet = Packet(vector_length,m)
        dct_plans = plans = plan_dct_iv(vector_length,max_depth)
        new(bell,dct_plans,packet,J,m,max_depth,vector_length)
    end
end
set!(t::Type{T},lcb::LocalCosineBasis,v::AbstractVector{Float64}) where T <: AbstractPacketPosition = set!(t,lcb.packet,v)
