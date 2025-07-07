# Property evaluated: Pure specific density [kg/m³]

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
exp_path = "C:/Users/beral/My Drive/usp-projects/2024-11-IL-Electrolytes/eIL-experimental/csv_rho_pure/"*cation*anion*".csv"
ions = [cation,anion]
model_pure = SAFTVRMie([cation,anion]);
print(model_pure)
print("\n")

####################
### EXPERIMENTAL ###
####################
exp_data = CSV.read(exp_path,DataFrame)
data_P = exp_data[:,1] # [Pa]
data_T = exp_data[:,2] # [K]
data_rho = exp_data[:,3] # [kg/m³]
N_points = length(data_P)
rho_min = minimum(data_rho)
rho_max = maximum(data_rho)
P_min = minimum(data_P*1E-6) # [MPa]
P_max = maximum(data_P*1E-6)
T_min = minimum(data_T)
T_max = maximum(data_T)
print("data_P,data_T,data_rho,datapoints,AARD%:")
print("\n")
@printf("%.1f",P_min)
print("–")
@printf("%.1f",P_max)
print(",")
@printf("%.2f",T_min)
print("–")
@printf("%.2f",T_max)
print(",")
@printf("%.1f",rho_min)
print("–")
@printf("%.1f",rho_max)
print(",")
print(N_points)

################
### MODELING ###
################
z_molar = fill([0.5, 0.5], N_points) # WE HAVE TO SPECIFY THE COMPOSITION FOR EVERY [PRESSURE,TEMPERATURE] ELEMENT!
model_rho = mass_density.(model_pure,data_P,data_T,z_molar)

# RELATIVE DEVIATION
RD_rho = model_rho - data_rho
print(",")
AARD_rho = 100*(mean(abs.(RD_rho ./ data_rho)))
@printf("%.4f",AARD_rho)
print("\n")
clipboard(AARD_rho)

################
### PLOTTING ###
################
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(figsize=(5,5), dpi = 311)

linestyles = ["solid", "dashed", "dotted", "dashdot"]
markers = ["s", "o", "^", "x"]
colors = ["black", "blue", "green", "red"]
labels = ["Experimental", "Model"]

# EXPERIMENTAL
PyPlot.plot(data_T,data_rho,label=labels[1],linestyle="",marker=markers[1],color=colors[1])
# MODEL
PyPlot.plot(data_T,model_rho,label=labels[2],linestyle=linestyles[1],marker="",color=colors[2])

PyPlot.legend(loc="best",frameon=true,fontsize=13)
PyPlot.xlabel("Temperature (K)",fontsize=16)
PyPlot.ylabel("Density (kg/m³)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)

fig_title = cation*anion
PyPlot.suptitle(fig_title, y=-0.01, fontsize=16, fontweight="bold")  # Adjust y to control the position

display(PyPlot.gcf())