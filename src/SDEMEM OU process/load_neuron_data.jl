# load package
using PyPlot

# define the ISI type
struct ISI
    time::Array{Float64}
    mV::Array{Float64}
end

# function to load the neuron data
function load_neuron_data()

    # load file
    file = open("data/neuron_data_312ISI.dat")

    # read file
    data_file = readlines(file)

    # pre-allocate vectors
    nbr_obs = length(data_file)
    mV_all = zeros(nbr_obs)
    t_all = zeros(nbr_obs)
    ISI_count = zeros(nbr_obs)

    # read each file line
    for i = 1:nbr_obs
        file_line = split(data_file[i])
        t_all[i] = parse(Float64, file_line[1])
        mV_all[i] = parse(Float64, file_line[2])
        ISI_count[i] = parse(Float64, file_line[3])
    end

    # set nbr ISIs
    nbr_ISE = ISI_count[end]

    scaling = 1000 # fix scaling (we should have the data in milliV)

    # find all individual ISIs
    mV_temp = mV_all[findall(x->x == 1, ISI_count)]*scaling
    t_temp = t_all[findall(x->x == 1, ISI_count)]
    ISI_all = [ISI(t_temp,mV_temp)]

    for i in 2:nbr_ISE
        mV_temp = mV_all[findall(x->x == i, ISI_count)]*scaling
        t_temp = t_all[findall(x->x == i, ISI_count)]
        append!(ISI_all, [ISI(t_temp,mV_temp)])
    end

    # return ISI and nbr_ISE
    return ISI_all, Int(nbr_ISE)

end

# load data
ISI_all, nbr_ISE = load_neuron_data()

# Plot data (this plot should be similar to Fig 3 in the paper)
PyPlot.figure()
for i in 1:nbr_ISE
    PyPlot.plot(ISI_all[i].time, ISI_all[i].mV, "k", alpha = 0.2)
end
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("depolarization mV")


# How to deal with the different time lengths for the different subjects?
# - I think this migth cause some annoying problems in the code...
# - Actually since the data lenght is different for different subjects it would make sense to have subject-dependet number of particels

# Which model should we use? Which parameters should be random effects and shared among subjects?

# Eval inference results?
#  - to check the results we can look at the posterior predictive dist
