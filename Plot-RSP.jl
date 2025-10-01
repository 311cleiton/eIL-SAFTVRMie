using Clapeyron

# 1. define your model, T, p, etc.
solvent = "water"
cation = "C4MIM"
anion = "BF4"
IL_pretty = "["*cation*"]"*"["*anion*"]"
model_neutral = SAFTVRMie([solvent,cation,anion])
model_electrolyte = SAFTVREMie([solvent],[cation,anion])
model_LinMixRSP = LinMixRSP([solvent],[cation,anion])
model_ZuoFurst = ZuoFurst([solvent],[cation,anion])
model_Schreckenberg = Schreckenberg([solvent],[cation,anion])
input_temperature = 298.15 # K
input_pressure = 100000.0 # Pa
# ε_experimental = [76.9, 71.9, 66.1, 57.5, 46.8, 42.5, 38.2, 27.5, 18.9, 13.1, 8.1, 7.4]
ε_experimental = [78.3, 76.9, 71.9, 66.1, 57.5, 46.8, 42.5, 38.2, 27.5, 18.9, 13.1, 8.1, 7.4, 6.7]

# 2. your measured volume fractions
# φ_IL = [0.02, 0.09, 0.17, 0.29, 0.44, 0.50, 0.56, 0.71, 0.83, 0.91, 0.98, 0.99] # for IL
# φ_IL = [0.0000001, 0.02, 0.09, 0.17, 0.29, 0.44, 0.50, 0.56, 0.71, 0.83, 0.91, 0.98, 0.99, 0.9999999] # for IL
φ_IL = collect(range(0.0001, 0.9999, step=0.0001))
# φ_IL = [0.0000001, 0.02, 0.09, 0.17, 0.29, 0.44, 0.50, 0.56, 0.71, 0.83, 0.91, 0.98, 0.99, 0.9999999] # for IL
model_φ = [[1 - x, x / 2, x / 2] for x in φ_IL] # volume fraction of [water, IL cation, IL anion]

# 3. compute pure-component molar volumes
x_water_pure = [1.0, 0.0, 0.0]
x_IL_pure = [0.0, 0.5, 0.5]
x_SA = [x_water_pure, x_IL_pure, x_IL_pure]
model_volume = volume.(model_neutral, input_pressure, input_temperature, x_SA; phase=:liquid) # [m³/mol]
model_volume = fill(model_volume, length(φ_IL))

# 4. convert to mole fractions
N_unnorm = []
for i in 1:length(φ_IL)
    append!(N_unnorm,[model_φ[i] ./ model_volume[i]])
end
model_x = []
for i in 1:length(φ_IL)
    append!(model_x,[N_unnorm[i] ./ sum(N_unnorm[i])])
end
x_IL = []
for i in 1:length(φ_IL)
    append!(x_IL,[model_x[i][2]+model_x[i][3]])
end

# 5. calculate RSP
# volume_RSP = volume.(model_neutral, input_pressure, input_temperature, model_x) # [m³]
volume_RSP = volume.(model_electrolyte, input_pressure, input_temperature, model_x) # [m³]
ε_ZuoFurst = Clapeyron.dielectric_constant.(model_ZuoFurst, volume_RSP, input_temperature, model_x)
ε_LinMixRSP = Clapeyron.dielectric_constant.(model_LinMixRSP, volume_RSP, input_temperature, model_x)
ε_Schreckenberg = Clapeyron.dielectric_constant.(model_Schreckenberg, volume_RSP, input_temperature, model_x)

# 6. plot
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
# PyPlot.figure(figsize=(5,5), dpi = 311)
PyPlot.figure(figsize=(9,5), dpi = 311)

labels = ["Koeberg et al.", "Zuo and Furst", "LMR (78.4 to 12)", "Schreckenberg et al."]
# labels = ["Schreckenberg", "ZuoFurst", "LinMixRSP", "conditionalRSP"]
linestyles = ["", "solid", "dashed", "dotted", "dashdot"]
# markers = ["s", "o", "^", "x"]
markers = ["s", "", "", ""]
colors = ["black", "blue", "green", "red"]

x_IL_exp_SAFT = [1.9216277898725135e-8, 0.003906369483846949, 0.018650650832239275, 0.03786819760709666, 0.0727768339719944, 0.1311789681929303, 0.16118836232171038, 0.19651014779943468, 0.31994414677883304, 0.4840590853599209, 0.6602082508711844, 0.9039936289531713, 0.9500602113616644, 0.9999994796080358]

PyPlot.plot(x_IL_exp_SAFT,ε_experimental,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
# PyPlot.plot(x_IL,ε_experimental,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])
PyPlot.plot(x_IL,ε_ZuoFurst,label=labels[2],linestyle=linestyles[2],marker=markers[2],color=colors[2])
PyPlot.plot(x_IL,ε_LinMixRSP,label=labels[3],linestyle=linestyles[3],marker=markers[3],color=colors[3])
PyPlot.plot(x_IL,ε_Schreckenberg,label=labels[4],linestyle=linestyles[4],marker=markers[4],color=colors[4])

PyPlot.legend(loc="center right",bbox_to_anchor=(1, 0.73),frameon=true,fontsize=16,facecolor="white")
# PyPlot.legend(loc="best",frameon=true,fontsize=16)
PyPlot.xlabel("Mole fraction of [C₄mim][BF₄]",fontsize=16)
PyPlot.ylabel("Dielectric constant",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)

xlim_min = 0.0
xlim_max = 1.0
# PyPlot.xticks(collect(xlim_min:0.1:xlim_max))
PyPlot.xlim([xlim_min,xlim_max])
ylim_min = 0
ylim_max = 80
# PyPlot.yticks(collect(ylim_min:25:ylim_max))
PyPlot.ylim([ylim_min,ylim_max])

PyPlot.grid(true, axis="y")  # turns on the grid

IL_cute = "(b)"
PyPlot.suptitle(IL_cute, y=-0.01, fontsize=16, fontweight="bold")  # Adjust y to control the position

display(PyPlot.gcf())

print(ε_Schreckenberg)  # this is your ε_Schreckenberg vector
print("\n")






# 6b. plot
PyPlot.clf()
PyPlot.rc("font", family="times new roman")
PyPlot.figure(figsize=(5,5), dpi = 311)
# PyPlot.figure(figsize=(9,5), dpi = 311)

labels = ["Koeberg et al. (experimental)", "Zuo and Furst", "LMR (78.4 to 12)", "Schreckenberg et al."]
# labels = ["Schreckenberg", "ZuoFurst", "LinMixRSP", "conditionalRSP"]
linestyles = ["", "solid", "dashed", "dotted", "dashdot"]
# markers = ["s", "o", "^", "x"]
markers = ["s", "", "", ""]
colors = ["black", "blue", "green", "red"]

vv_IL = [0.0000001, 0.02, 0.09, 0.17, 0.29, 0.44, 0.50, 0.56, 0.71, 0.83, 0.91, 0.98, 0.99, 0.9999999] # for IL

PyPlot.plot(vv_IL,ε_experimental,label=labels[1],linestyle=linestyles[1],marker=markers[1],color=colors[1])

# PyPlot.legend(loc="center right",bbox_to_anchor=(1, 0.73),frameon=true,fontsize=16,facecolor="white")
PyPlot.legend(loc="best",frameon=true,fontsize=16)
# PyPlot.legend(loc="upper right",frameon=true,fontsize=16)
PyPlot.xlabel("Volume fraction of [C₄mim][BF₄]",fontsize=16)
PyPlot.ylabel("Dielectric constant",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)

xlim_min = 0.0
xlim_max = 1.0
# PyPlot.xticks(collect(xlim_min:0.1:xlim_max))
PyPlot.xlim([xlim_min,xlim_max])
ylim_min = 0
ylim_max = 80
# PyPlot.yticks(collect(ylim_min:25:ylim_max))
PyPlot.ylim([ylim_min,ylim_max])

PyPlot.grid(true, axis="y")  # turns on the grid

IL_cute = "(a)"
PyPlot.suptitle(IL_cute, y=-0.01, fontsize=16, fontweight="bold")  # Adjust y to control the position

display(PyPlot.gcf())