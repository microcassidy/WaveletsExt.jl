mutable struct VariableVector <: AbstractVector{Float64}
    data::Vector{Float64}
    length::Int
    max_length::Int
    VariableVector(m::Int) = new(zeros(Float64,m),m,m)
end
_update_size!(v::VariableVector,m::Int) = 1<=m<=v.max_length ? setfield!(v,:length,m) : error("outside bounds of 1 and $(v.max_length)")
Base.length(v::VariableVector) = v.length
Base.size(v::VariableVector) = (v.length,)
Base.getindex(v::VariableVector,i::Int) = Base.getindex(v.data,i)
Base.setindex!(v::VariableVector,val::Float64,ind::Int) = setindex!(v.data,val,ind)
Base.show(v::VariableVector) = show(v.data)

Base.firstindex(v::VariableVector) = 1
Base.lastindex(v::VariableVector) = v.length


mutable struct Packet
    left::VariableVector
    right::VariableVector
    centre::VariableVector
    Packet(m::Int) = new(VariableVector(m), VariableVector(m), VariableVector(m))
end
reset!(a::AbstractArray{<:Number}) = fill!(a,0)
_update_size!(p::Packet,m::Int) = map!(Base.Fix{2}(_update_size!,m), [p.left,p.right,p.centre])
reset!(p::Packet) = map!(reset!, [p.left,p.right,p.centre])
set!(p::Packet,v,which) = setfield!(p,which,v)
#foo
#bea
#ranm
