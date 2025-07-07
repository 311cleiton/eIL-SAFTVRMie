# ESTIMATION OF segment(m), sigma(σ), epsilon(ε) OF AN IL CATION AND ANION

# STEP 1 - PACKAGES
using Clapeyron
using Metaheuristics
using Statistics
using Printf
using Dates
print(now())
print("\n")

# STEP 2 - MODEL INPUT
cation = "C4MIM"
anion = "BF4"
ions = [cation,anion]
molar_composition = [0.5, 0.5]
model_PE = SAFTVRMie(ions);
print(model_PE)
print("\n")

@time begin
# STEP 3 - PARAMETERS TO FIT
toestimate = [
    # cation
    Dict(
        :param => :segment,
        :indices => 1,
        :lower => 4.0-3.0,
        :upper => 4.0+3.0,
        :guess => 4.0
    ),
    Dict(
        :param => :sigma, # [Å]
        :factor => 1E-10, # convert [Å] to [m]
        :indices => (1,1), # (1,1)=cation; (2,2)=anion
        :lower => 3.0-2.0,
        :upper => 3.0+2.0,
        :guess => 3.0
    ),
    Dict(
        :param => :epsilon, # [K]
        :indices => (1,1),
        :lower => 300.0-200.0,
        :upper => 300.0+200.0,
        :guess => 300.0
    ),

    # anion
    Dict(
        :param => :segment,
        :indices => 2,
        :lower => 4.0-3.0,
        :upper => 4.0+3.0,
        :guess => 4.0
    ),
    Dict(
        :param => :sigma,
        :factor => 1E-10,
        :indices => (2,2),
        :lower => 3.0-2.0,
        :upper => 3.0+2.0,
        :guess => 3.0
    ),
    Dict(
        :param => :epsilon,
        :indices => (2,2),
        :lower => 300.0-200.0,
        :upper => 300.0+200.0,
        :guess => 300.0
    ),
];

# STEP 4 - PROPERTY USED FOR ESTIMATION
function mass_rho(model_PE::EoSModel,p,T)
    md = mass_density(model_PE,p,T,molar_composition) # p = [Pa]; T = [K]
    return md[1] # md = [kg/m³]
end

# STEP 5 - ESTIMATOR
path_rho = "C:/Users/beral/My Drive/usp-projects/2024-11-IL-Electrolytes/eIL-experimental/csv_PE_rho/"*cation*anion*".csv"
estimator,objective,initial,upper,lower = Estimation(model_PE,toestimate,[path_rho]);
local_method = ECA(;options=Options(iterations=100))
local_params, local_model = optimize(objective, estimator, local_method);

# PRINT PARAMETERS
print(cation,"+",anion," [m,sigma,epsilon] =")
print("\n")
@printf("%.4f",local_params[1])
print(",")
@printf("%.4f",local_params[2])
print(",")
@printf("%.2f",local_params[3])
print("\n")
@printf("%.4f",local_params[4])
print(",")
@printf("%.4f",local_params[5])
print(",")
@printf("%.2f",local_params[6])
print("\n")

print("Runtime:")
end#@time
print("\n")