# Import libraries
import moments, moments.LD
import demes, demesdraw
import numpy as np
import msprime as ms

# https://tskit.dev/msprime/docs/stable/demography.html#sec-demography-numerical-trajectories

def get_iirc(demes_model,pop,T=None):
    """
    Returns two arrays: the Coalescence rate and the Inferred Inverse Coalescence Rate (Popsize)
    """
    m = ms.Demography.from_demes(demes_model)
    debug = m.debug()
    if T == None:
        T = np.concatenate([
            np.linspace(0, 1000, 2001),
            np.linspace(1000, 1.5e4, 401)[1:]
        ])
    R, _ = debug.coalescence_rate_trajectory(T, {pop: 2})
    inversed_R = 1/(2*R)

    return R,inversed_R,T

# Population N change model
def get_N_times_from_IIRC(IIRC,T):
    """
    Takes an the inverse inferred rate of coalescence and the time points and returns the times at which size changes.
    """
    previous_N = IIRC[0]
    Ns = []
    times = []
    for N,time in zip(np.flip(IIRC),np.flip(T)):
        if int(N) != previous_N:
            Ns.append(previous_N)
            times.append(time)
            previous_N = int(N)
    return Ns,times

def size_change_from_IIRC(IIRC,T,yaml_filename=None,plot=True,plot_filename=None):
    """
    Takes a vector of Ns in form of IIRC and T times to create a model
    """

    Ns,times = get_N_times_from_IIRC(IIRC,T)
    
    m = demes.Builder()

    #current_time=np.array(time_period).sum()

    epochs = []
    for N,time in zip(Ns,times):
        epoch = dict(start_size=N,end_time=time)
        epochs.append(epoch)

    m.add_deme("Deme",epochs=epochs)

    # Resolve model
    graph = m.resolve()

    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        p = demesdraw.tubes(graph, log_time=False, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph 