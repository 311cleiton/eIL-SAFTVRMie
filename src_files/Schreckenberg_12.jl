struct Schreckenberg_12 <: SchreckenbergModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::SchreckenbergParam
    references::Array{String,1}
end

@registermodel Schreckenberg_12
export Schreckenberg_12

function dielectric_constant_smoothed(model::Schreckenberg_12, V, T, z, _data=nothing)
    eps_floor = 12.0
    # --- compute raw result as before, but avoid early hard return ---
    d_T = model.params.d_T.values
    d_V = model.params.d_V.values
    Z = model.params.charge.values
    ineutral = model.icomponents[Z .== 0]

    n_solv = zero(first(z))
    for i in ineutral
        n_solv += z[i]
    end

    # if no neutrals, set raw_result to 1.0 (physically: no solvent contribution)
    if isempty(ineutral) || n_solv == 0
        raw = 1.0
    else
        ρ_solv = n_solv / V
        dbar = zero(T + first(z))
        for i in ineutral
            di = d_V[i] * (d_T[i]/T - 1)
            for j in ineutral
                # use symmetric d_ij; when i==j this is fine
                dj = d_V[j] * (d_T[j]/T - 1)
                dij = 0.5*(di + dj)
                dbar += dij * z[i] * z[j]
            end
        end
        dbar /= (n_solv * n_solv)
        raw = 1 + ρ_solv * dbar
    end

    # smoothing parameter (in epsilon units). Smaller -> sharper transition.
    delta = 0.5   # try 0.2-1.0 for sensitivity testing
    w = 0.5*(1 + tanh((raw - eps_floor)/delta))
    epsr = eps_floor + (raw - eps_floor) * w

    return epsr
end
