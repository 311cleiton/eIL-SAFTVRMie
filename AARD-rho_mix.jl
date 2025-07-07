# Property evaluated: mixture mass_density (kg/m³)

################
### PACKAGES ###
################
using Clapeyron
using CSV
using DataFrames
using Statistics
using Printf
using PyCall
import PyPlot
using Dates
print(now())
print("\n")

#############
### INPUT ###
#############
solvent = "water"
# solvent = "methanol"
# solvent = "ethanol"
cation_all = [
    "C2MIM" # 1
    "C4MIM" # 2
    "C6MIM" # 3
    "C8MIM" # 4
]
anion_all = [
    "SCN"   # 1
    "BF4"   # 2
    "PF6"   # 3
    "TFO"   # 4
    "TF2N"  # 5
]
cation = cation_all[2]
anion =   anion_all[2]
model_e = SAFTVREMie([solvent],[cation,anion])
# model_e = eSAFTVRMie([solvent],[cation,anion])
# model_e = SAFTVRMie([solvent,cation,anion]) # neutral model
# model_RSP = Schreckenberg([solvent],[cation,anion])
# model_RSP = LinMixRSP([solvent],[cation,anion])
print(model_e,"\n")
# include("C:/Users/beral/My Drive/coding/CL03/Clapeyron.jl/src/models/Electrolytes/RSP/Schreckenberg.jl")

IL = cation*anion
salts = [(IL,(cation=>1,anion=>1))]
IL_pretty = "["*cation*"]"*"["*anion*"]"
fig_title = solvent*"–"*IL_pretty

####################
### EXPERIMENTAL ###
####################
exp_path = "C:/Users/beral/My Drive/usp-projects/2024-11-IL-Electrolytes/eIL-experimental/csv_rho_mix/"*solvent*"/"*cation*anion*".csv"
exp_data = CSV.read(exp_path,DataFrame)
data_P = exp_data[:,1] # pressure [Pa]
data_T = exp_data[:,2] # temperature [K]
data_x = exp_data[:,3] # composition [mol/mol]
data_rho = exp_data[:,4] # mass density [kg/m³]
N_points = length(data_rho)
P_min = minimum(data_P*1E-6) # pressure [MPa]
P_max = maximum(data_P*1E-6)
T_min = minimum(data_T) # temperature [K]
T_max = maximum(data_T)
x_min = minimum(data_x) # composition [mol/mol]
x_max = maximum(data_x)
rho_min = minimum(data_rho) # mass density [kg/m³]
rho_max = maximum(data_rho)
print("data_P,data_T,data_x,data_rho,datapoints:")
print("\n")
@printf("%.1f",P_min)
print("–")
@printf("%.1f",P_max)
print(",")
@printf("%.2f",T_min)
print("–")
@printf("%.2f",T_max)
print(",")
@printf("%.4f",x_min)
print("–")
@printf("%.4f",x_max)
print(",")
@printf("%.2f",rho_min)
print("–")
@printf("%.2f",rho_max)
print(",")
print(N_points)

@time begin
################
### MODELING ###
################
model_x = [[1 - x, x / 2, x / 2] for x in data_x]
model_rho = mass_density.(model_e, data_P, data_T, model_x) # [kg/m³]

# RELATIVE DEVIATION
model_rho = map(x -> x[1], model_rho)
RD_rho = model_rho - data_rho
AARD_rho = 100*(mean(abs.(RD_rho ./ data_rho)))
print("\n")
print("%AARD:")
print("\n")
@printf("%.3f",AARD_rho)
print("\n")
clipboard(AARD_rho)
print("Runtime:")
end#@time

################
### PLOTTING ###
################
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(figsize=(5,5), dpi = 311)

linestyles = ["solid", "dashed", "dotted", "dashdot"]
markers = ["s", "o", "^", "x"]
colors = ["black", "blue", "green", "red"]

PyPlot.plot(data_x,data_rho,label="Experimental",linestyle="",marker=markers[1],color=colors[1])
PyPlot.plot(data_x,model_rho,label="Model",linestyle=linestyles[1],marker="",color=colors[2])

PyPlot.legend(loc="best",frameon=true,fontsize=16)

PyPlot.xlabel("IL mole fraction",fontsize=16)
PyPlot.ylabel("Specific density (kg/m³)",fontsize=16)

PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)

PyPlot.suptitle(fig_title, y=-0.01, fontsize=16, fontweight="bold")  # Adjust y to control the position

display(PyPlot.gcf())