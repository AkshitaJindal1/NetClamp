    using Plots, MAT, Statistics

    a = matopen("latest_experiment_4.mat")
    b = read(a)

    ##
    str = "2026_05_27_12_23"
        time = (0:(length(b["response_opsin_light_only_1_"*str*"_12__data2"])-1)).*0.15
        timet = time[1:length([zeros(1,399) -b["target_currents_4_27_signal"]']')]
        timer = time[1:size([zeros(300,399) b["resampling_trajectories"]],2)]

    off = mean(b["response_opsin_light_only_2_"*str*"_14__data2"][1:350]./1000.0)
    plot(timer, [zeros(300,399) b["resampling_trajectories"]]'.+off, linewidth=4, legend=false, color="grey",alpha=0.1, xlabel = "time [ms]")

    plot!(time,b["response_opsin_light_only_1_"*str*"_12__data2"]./1000.0)
    plot!(time,b["response_opsin_light_only_2_"*str*"_14__data2"]./1000.0)
    plot!(time,b["response_opsin_light_only_3_"*str*"_16__data2"]./1000.0)
    plot!(time,b["response_opsin_light_only_4_"*str*"_19__data2"]./1000.0)
    plot!(time,b["response_opsin_light_only_5_"*str*"_21__data2"]./1000.0)
    plot!(timet, [zeros(1,399) -b["target_currents_4_27_signal"]']'.+off, linewidth=4, legend=false, color="black", xlabel = "time [ms]", ylabel="PC [nA]")
    plot(b["Light_signal_4_27"][:,1], b["Light_signal_4_27"][:,2], linewidth=4, legend=false, color="black", xlabel = "time [ms]", ylabel="PC [nA]")



    

##
    str = "2026_05_27_12_26"
    off = mean(b["response_opsin_light_only_3_"*str*"_05__data2"][1:350]./1000.0)
    plot(b["response_opsin_light_only_1_"*str*"_00__data2"]./1000.0)
    plot!(b["response_opsin_light_only_2_"*str*"_03__data2"]./1000.0)
    plot!(b["response_opsin_light_only_3_"*str*"_05__data2"]./1000.0)
    plot!(b["response_opsin_light_only_4_"*str*"_07__data2"]./1000.0)
    plot!(b["response_opsin_light_only_5_"*str*"_09__data2"]./1000.0)
    plot!([zeros(1,399) -b["target_currents_4_28_signal"]']'.+off, linewidth=4, legend=false, color="black")


##
    str = "2026_05_27_12_27"
    off = mean(b["response_opsin_light_only_3_"*str*"_08__data2"][1:350]./1000.0)
    plot(b["response_opsin_light_only_1_"*str*"_04__data2"]./1000.0)
    plot!(b["response_opsin_light_only_2_"*str*"_06__data2"]./1000.0)
    plot!(b["response_opsin_light_only_3_"*str*"_08__data2"]./1000.0)
    plot!(b["response_opsin_light_only_4_"*str*"_11__data2"]./1000.0)
    plot!(b["response_opsin_light_only_5_"*str*"_13__data2"]./1000.0)
    plot!([zeros(1,399) -b["target_currents_4_29_signal"]']'.+off, linewidth=4, legend=false, color="black")

##
    str = "2026_05_27_12_28"
    off = mean(b["response_opsin_light_only_3_"*str*"_23__data2"][1:350]./1000.0)
    plot(b["response_opsin_light_only_1_"*str*"_19__data2"]./1000.0)
    plot!(b["response_opsin_light_only_2_"*str*"_21__data2"]./1000.0)
    plot!(b["response_opsin_light_only_3_"*str*"_23__data2"]./1000.0)
    plot!(b["response_opsin_light_only_4_"*str*"_25__data2"]./1000.0)
    plot!(b["response_opsin_light_only_5_"*str*"_28__data2"]./1000.0)
    plot!([zeros(1,399) -b["target_currents_4_30_signal"]']'.+off, linewidth=4, legend=false, color="black")





##
    str = "2026_05_27_12_48"
    off = mean(b["response_opsin_light_only_3_"*str*"_50__data2"][1:350]./1000.0)
    plot(b["response_opsin_light_only_1_"*str*"_45__data2"]./1000.0)
    plot!(b["response_opsin_light_only_2_"*str*"_47__data2"]./1000.0)
    plot!(b["response_opsin_light_only_3_"*str*"_50__data2"]./1000.0)
    plot!(b["response_opsin_light_only_4_"*str*"_52__data2"]./1000.0)
    plot!(b["response_opsin_light_only_5_"*str*"_54__data2"]./1000.0)
    plot!([zeros(1,399) -b["target_currents_4_31_signal"]']'.+off, linewidth=4, legend=false, color="black")



##################
using DelimitedFiles
str = "Exp4_6_27_May/"
list = readdir(str)
responses = filter(contains("response"), list)
#target = filter(contains("Current_target"),list)
target = filter(contains("target_current"),list)
tar = - read(matopen(str*target[1]))["signal"]

resamp = filter(contains("resampling_trajectories"),list)
res = readdlm(str*resamp[1],',')
time = 0.15.*(0:(size([zeros(300,399) res],2)-1))
p = plot(time, [zeros(300,399) res]'.+off, linewidth=4, legend=false, color="grey",alpha=0.1)

#tar = -readdlm(str*target[1],',')

off = mean(read(matopen(str*responses[1]))["data2"][1:100]./1000.0)
    time = 0.15.*(0:(size([zeros(1,399) tar'],2)-1))

plot!(time, [zeros(1,399) tar']'.+off, linewidth=4, legend=false,color="black")

for i in responses
    a = read(matopen(str*i))
        time = 0.15.*(0:(length(a["data2"]'./1000.0)-1))

    plot!(time, a["data2"]'./1000.0, legend=false, xlabel = "time [ms]", ylabel="PC [nA]")
end

##################
using DelimitedFiles
str = "Exp4_8_27_May/"
list = readdir(str)
responses = filter(contains("response"), list)
#target = filter(contains("Current_target"),list)
target = filter(contains("target_current"),list)
tar = - read(matopen(str*target[1]))["signal"]
resamp = filter(contains("resampling_trajectories"),list)
res = readdlm(str*resamp[1],',')
#tar = -readdlm(str*target[1],',')

off = mean(read(matopen(str*responses[1]))["data2"][1:100]./1000.0)
time = 0.15.*(0:(size([zeros(300,399) res],2)-1))

p = plot(time, [zeros(300,399) res]'.+off, linewidth=4, legend=false, color="grey",alpha=0.1)
    time = 0.15.*(0:(size([zeros(1,399) tar'],2)-1))

plot!(time, [zeros(1,399) tar']'.+off, linewidth=4, legend=false,color="black")

for i in responses
   a = read(matopen(str*i))
        time = 0.15.*(0:(length(a["data2"]'./1000.0)-1))

    plot!(time, a["data2"]'./1000.0, legend=false, xlabel = "time [ms]", ylabel="PC [nA]")

end

##################
using DelimitedFiles
str = "Exp4_9_27_May/"
list = readdir(str)
responses = filter(contains("response"), list)
#target = filter(contains("Current_target"),list)
target = filter(contains("target_current"),list)
tar = - read(matopen(str*target[1]))["signal"]
#tar = -readdlm(str*target[1],',')

off = mean(read(matopen(str*responses[1]))["data2"][1:100]./1000.0)
p = plot([zeros(1,399) tar']'.+off, linewidth=4, legend=false, color="black")

for i in responses
    a = read(matopen(str*i))
    plot!(a["data2"]'./1000.0, legend=false)
end


##################
using DelimitedFiles
str = "Exp1_2_27_May/"
list = readdir(str)
responses = filter(contains("response"), list)
#target = filter(contains("Current_target"),list)
target = filter(contains("target_current"),list)
tar = - read(matopen(str*target[1]))["signal"]
tim =  read(matopen(str*target[1]))["time"]

#tar = -readdlm(str*target[1],',')

off = mean(read(matopen(str*responses[1]))["data2"][1:100]./1000.0)
p = plot([zeros(1,399) tar']'.+off, linewidth=4, legend=false, color="black")

for i in responses
    a = read(matopen(str*i))
    plot!(a["data2"]'./1000.0, legend=false)
end



##################
using DelimitedFiles
str = "Exp3_0_28_May/"
list = readdir(str)
responses = filter(contains("response"), list)
#target = filter(contains("Current_target"),list)
target = filter(contains("target_current"),list)
tar =  read(matopen(str*target[1]))["signal"]
tim =  read(matopen(str*target[1]))["time"]

#tar = -readdlm(str*target[1],',')

off = mean(read(matopen(str*responses[1]))["data2"][1:100]./1000.0)
p = plot(tim.+400.0, tar.+off, linewidth=4, legend=false, color="black")

for i in responses
    a = read(matopen(str*i))
    plot!(a["data2"]'./1000.0, legend=false)
end
