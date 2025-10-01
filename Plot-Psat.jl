print("\n")
print(now(),"\n")
# vapor pressure prediction

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

model_1 = SAFTVRMie([cation_all[1],anion_all[2]])
model_2 = SAFTVRMie([cation_all[2],anion_all[2]])
model_3 = SAFTVRMie([cation_all[3],anion_all[2]])
model_4 = SAFTVRMie([cation_all[4],anion_all[2]])

################
### MODELING ###
################
z_molar = [0.5, 0.5]
data_T_1 = 411.89 # [K] [C₂mim][BF₄]
data_T_2 = 411.84 # [K] [C₄mim][BF₄]
data_T_3 = 409.98 # [K] [C₆mim][BF₄]
data_T_4 = 410.44 # [K] [C₈mim][BF₄]
data_Psat_1 = 10 # [μPa] [C₂mim][BF₄]
data_Psat_2 = 14 # [μPa] [C₄mim][BF₄]
data_Psat_3 = 7 # [μPa] [C₆mim][BF₄]
data_Psat_4 = 4 # [μPa] [C₈mim][BF₄]

Pbub_1, V_l_1, V_v_1, y_1 = bubble_pressure(model_1, data_T_1, z_molar)
Pbub_2, V_l_2, V_v_2, y_2 = bubble_pressure(model_2, data_T_2, z_molar)
Pbub_3, V_l_3, V_v_3, y_3 = bubble_pressure(model_3, data_T_3, z_molar)
Pbub_4, V_l_4, V_v_4, y_4 = bubble_pressure(model_4, data_T_1, z_molar)
Psat_1 = Pbub_1*1e-6 # μPa
Psat_2 = Pbub_2*1e-6 # μPa
Psat_3 = Pbub_3*1e-6 # μPa
Psat_4 = Pbub_4*1e-6 # μPa

################
### PLOTTING ###
################
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(figsize=(5,4), dpi = 311)

linestyles = ["", "solid", "dashed", "dotted", "dashdot"]
markers = ["s", "o", "^", "x"]
colors = ["black", "blue", "green", "red", "orange"]
labels = ["", "Experimental", "SAFT-VR Mie", "SAFTVRMie {Joback}", "Joback"]

chain_n = [2, 4, 6, 8]

PyPlot.plot(chain_n[1],data_Psat_1,label=labels[2],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(chain_n[1],Psat_1,label=labels[3],linestyle=linestyles[1],marker=markers[2],color=colors[1])

PyPlot.plot(chain_n[1],data_Psat_1,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(chain_n[1],Psat_1,label=labels[1],linestyle=linestyles[1],marker=markers[2],color=colors[1])

PyPlot.plot(chain_n[2],data_Psat_2,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(chain_n[2],Psat_2,label=labels[1],linestyle=linestyles[1],marker=markers[2],color=colors[1])

PyPlot.plot(chain_n[3],data_Psat_3,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(chain_n[3],Psat_3,label=labels[1],linestyle=linestyles[1],marker=markers[2],color=colors[1])

PyPlot.plot(chain_n[4],data_Psat_4,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(chain_n[4],Psat_4,label=labels[1],linestyle=linestyles[1],marker=markers[2],color=colors[1])

PyPlot.legend(loc="best",frameon=true,fontsize=16)

colnames = ["n in [Cₙmim][BF₄]", "Vapor pressure (μPa)"]
PyPlot.xlabel(colnames[1],fontsize=16)
PyPlot.ylabel(colnames[2],fontsize=16)

PyPlot.xticks(collect(2:2:8))
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)

# IL_pretty = "[C₄mim][BF₄]"
# PyPlot.suptitle(IL_pretty, y=-0.01, fontsize=16, fontweight="bold")  # Adjust y to control the position

display(PyPlot.gcf())

print("\n")
