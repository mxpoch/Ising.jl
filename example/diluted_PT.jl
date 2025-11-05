# parallel_tempering.jl

using Carlo
using Carlo.JobTools
using Ising

tm = TaskMaker()

tm.sweeps = 20000
tm.thermalization = 2000
tm.binsize = 100

Ts = range(1.0, 3, 30)

tm.parallel_tempering = (mc = Ising.MC, parameter = :T, values = Ts, interval = 1)

# initializing the system
lattice_param_list = []
for dp in DP
    # setting iterations per grid size
    for itnum in [(20, 10), (10, 20), (4,40)]
    # for itnum in [(20, 6)]
        for i in 1:itnum[1] 
            # grid size 
            tm.Lx = itnum[2]
            tm.Ly = itnum[2]
            # dilution percent
            tm.DP = dp 
            task(tm)
        end
    end
end

job = JobInfo(
    splitext(@__FILE__)[1],
    ParallelTemperingMC; # the underlying model MC is set in tm.parallel_tempering.mc
    run_time = "24:00:00",
    checkpoint_time = "30:00",
    tasks = make_tasks(tm),
    ranks_per_run = length(tm.parallel_tempering.values), # needs to match!
)

start(job, ARGS)
