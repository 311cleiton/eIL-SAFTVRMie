# Property evaluated: mixture speed_of_sound (m/s)

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
print(now(),"\n")

#############
### INPUT ###
#############
# solvent = "water"
solvent = "methanol"
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

model_e = SAFTVREMie([solvent],[cation,anion]; idealmodel = JobackIdeal)
# model_e = eSAFTVRMie([solvent],[cation,anion]; idealmodel = JobackIdeal)
model_J = JobackIdeal([solvent,cation,anion])
model_N = SAFTVRMie([solvent,cation,anion]; idealmodel = JobackIdeal)
print(model_e,"\n")
# model_RSP = Schreckenberg([solvent],[cation,anion])

IL = cation*anion
salts = [(IL,(cation=>1,anion=>1))]
IL_pretty = "["*cation*"]"*"["*anion*"]"

####################
### EXPERIMENTAL ###
####################
exp_path = "C:/Users/beral/My Drive/usp-projects/2024-11-IL-Electrolytes/eIL-experimental/csv_sound_mix/"*solvent*"/"*cation*anion*".csv"
exp_data = CSV.read(exp_path,DataFrame)
data_P = exp_data[:,1] # [Pa]
data_T = exp_data[:,2] # [K]
data_x = exp_data[:,3] # [mol/mol]
data_sound = exp_data[:,4] # [m/s]
P_min = minimum(data_P*1E-6) # [MPa]
P_max = maximum(data_P*1E-6)
T_min = minimum(data_T) # [K]
T_max = maximum(data_T)
x_min = minimum(data_x) # [mol/mol]
x_max = maximum(data_x)
sound_min = minimum(data_sound) # [m/s]
sound_max = maximum(data_sound)
N_points = length(data_sound)
print("data_P,data_T,data_x,data_sound,datapoints:")
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
@printf("%.1f",sound_min)
print("–")
@printf("%.1f",sound_max)
print(",")
print(N_points)

@time begin
################
### MODELING ###
################
model_salts = fill(salts, N_points)
model_x = [[1 - x, x / 2, x / 2] for x in data_x]

model_sound = speed_of_sound.(model_e, data_P, data_T, model_x)
Joback_sound = speed_of_sound.(model_J, data_P, data_T, model_x)
neutral_sound = speed_of_sound.(model_N, data_P, data_T, model_x)

# clipboard(data_x)
# clipboard(data_sound)
# clipboard(neutral_sound)
# clipboard(model_sound)

############
### AARD ###
############
# electrolyte
model_sound = map(x -> x[1], model_sound)
RD_sound = model_sound - data_sound
AARD_sound = 100*(mean(abs.(RD_sound ./ data_sound)))
# neutral
neutral_sound = map(x -> x[1], neutral_sound)
RD_neutral = neutral_sound - data_sound
AARD_neutral = 100*(mean(abs.(RD_neutral ./ data_sound)))
print("\n")
print("%AARD electrolyte = ")
@printf("%.4f",AARD_sound)
print("\n")
print("%AARD neutral = ")
@printf("%.4f",AARD_neutral)
print("\n")
# clipboard(AARD_sound)
print("Runtime:")
end#@time
# print("\n")

################
### PLOTTING ###
################
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(figsize=(5,5), dpi = 311)

linestyles = ["", "solid", "dashed", "dotted", "dashdot"]
markers = ["s", "o", "^", "x"]
colors = ["black", "blue", "green", "red"]
labels = ["Experimental", "eModel", "SAFTVRMie {Joback}", "Joback"]

PyPlot.plot(data_x,data_sound,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(data_x,model_sound,label=labels[2],linestyle=linestyles[2],marker="",color=colors[2])
PyPlot.plot(data_x,neutral_sound,label=labels[3],linestyle=linestyles[3],marker="",color=colors[3])
PyPlot.plot(data_x,Joback_sound,label=labels[4],linestyle=linestyles[4],marker="",color=colors[4])

PyPlot.legend(loc="best",frameon=true,fontsize=13)
PyPlot.xlabel("IL mole fraction",fontsize=16)
PyPlot.ylabel("Speed of sound (m/s)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)

fig_title = solvent*"–"*IL_pretty
PyPlot.suptitle(fig_title, y=-0.01, fontsize=16, fontweight="bold")  # Adjust y to control the position

display(PyPlot.gcf())
