# Property evaluated: Speed of sound (pure) [m/s]

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
cation = cation_all[1]
anion =   anion_all[1]
model_neutral = SAFTVRMie([cation,anion]; idealmodel = JobackIdeal)
# model_neutral = SAFTVRMie([cation,anion]; idealmodel = MonomerIdeal)
model_Joback = JobackIdeal([cation,anion])
model_MonomerIdeal = MonomerIdeal([cation,anion])
print(model_neutral,"\n")

IL = cation*anion
salts = [(IL,(cation=>1,anion=>1))]
IL_pretty = "["*cation*"]"*"["*anion*"]"

####################
### EXPERIMENTAL ###
####################
exp_path = "C:/Users/beral/My Drive/usp-projects/2024-11-IL-Electrolytes/eIL-experimental/csv_sound_pure/"*cation*anion*".csv"
exp_data = CSV.read(exp_path,DataFrame)
data_P = exp_data[:,1] # [Pa]
data_T = exp_data[:,2] # [K]
data_sound = exp_data[:,3] # [m/s]
P_min = minimum(data_P*1E-6) # [MPa]
P_max = maximum(data_P*1E-6)
T_min = minimum(data_T) # [K]
T_max = maximum(data_T)
sound_min = minimum(data_sound) # [m/s]
sound_max = maximum(data_sound)
N_points = length(data_sound)
print("data_P,data_T,data_sound,datapoints,%AARD:")
print("\n")
@printf("%.1f",P_min)
print("–")
@printf("%.1f",P_max)
print(",")
@printf("%.2f",T_min)
print("–")
@printf("%.2f",T_max)
print(",")
@printf("%.1f",sound_min)
print("–")
@printf("%.1f",sound_max)
print(",")
print(N_points)

################
### MODELING ###
################
z_molar = fill([0.5, 0.5], N_points) # WE HAVE TO SPECIFY THE COMPOSITION FOR EVERY [PRESSURE,TEMPERATURE] ELEMENT!
neutral_sound = speed_of_sound.(model_neutral, data_P, data_T, z_molar)
Joback_sound = speed_of_sound.(model_Joback, data_P, data_T, z_molar)
MonomerIdeal_sound = speed_of_sound.(model_MonomerIdeal, data_P, data_T, z_molar)

############
### AARD ###
############
neutral_sound = map(x -> x[1], neutral_sound)
RD_sound = neutral_sound - data_sound
AARD_sound = 100*(mean(abs.(RD_sound ./ data_sound)))
print(",")
@printf("%.4f",AARD_sound)
print("\n")
clipboard(AARD_sound)

################
### PLOTTING ###
################
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(figsize=(5,5), dpi = 311)

linestyles = ["", "solid", "dashed", "dotted", "dashdot"]
markers = ["s", "o", "^", "x"]
colors = ["black", "blue", "green", "red"]
labels = ["Experimental", "λᵣ = 16", "λᵣ = 12"]

PyPlot.plot(data_T,data_sound,label="Experimental",linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(data_T,neutral_sound,label="SAFT-VR Mie {Joback}",linestyle=linestyles[2],marker="",color=colors[2])
PyPlot.plot(data_T,Joback_sound,label="Joback (soundⁱᵈ)",linestyle=linestyles[3],marker="",color=colors[3])
PyPlot.plot(data_T,MonomerIdeal_sound,label="MonomerIdeal (soundⁱᵈ)",linestyle=linestyles[4],marker="",color=colors[4])

PyPlot.legend(loc="best",frameon=true,fontsize=13, ncol=1)
PyPlot.xlabel("Temperature (K)",fontsize=16)
PyPlot.ylabel("Speed of sound (m/s)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)
PyPlot.suptitle(IL_pretty, y=-0.01, fontsize=16, fontweight="bold")  
display(PyPlot.gcf())
