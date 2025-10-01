struct Schreckenberg_12 <: SchreckenbergModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::SchreckenbergParam
    references::Array{String,1}
end

@registermodel Schreckenberg_12
export Schreckenberg_12

function Schreckenberg_12(solvents,ions; userlocations::Vector{String}=String[], verbose::Bool=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    components = format_components(components)
    icomponents = 1:length(components)

    params = getparams(components, ["Electrolytes/RSP/Schreckenberg.csv","Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose, ignore_missing_singleparams=["d_T","d_V"])
    d_T = params["d_T"]
    d_V = params["d_V"]
    charge = params["charge"]
    packagedparams = SchreckenbergParam(d_T,d_V,charge)

    references = String[]
    
    model = Schreckenberg_12(components,icomponents,packagedparams,references)
    return model
end

function dielectric_constant(model::Schreckenberg_12,V,T,z,_data=nothing)
    d_T = model.params.d_T.values
    d_V = model.params.d_V.values
    Z = model.params.charge.values
    ineutral = model.icomponents[Z.==0]

    if isempty(ineutral)
        return 12.0   # enforce boundary directly here
    end
    
    n_solv = zero(first(z))
    for i in ineutral
        n_solv += z[i]
    end
    ρ_solv = n_solv / V
    d̄ = zero(T+first(z))

    for i in ineutral
        di = d_V[i]*(d_T[i]/T-1)
        dij,zi = di,z[i]
        d̄ += dij*zi*zi
        for j in ineutral[ineutral.!=i]
            dj = d_V[j]*(d_T[j]/T-1)
            dij,zj = 0.5*(di+dj),z[j]
            d̄ += dij*zi*zj
        end
    end

    d̄ = d̄/(n_solv*n_solv)
    result = 1+ρ_solv*d̄

    return result < 12 ? 12.0 : result
end
