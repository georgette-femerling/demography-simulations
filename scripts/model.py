# Import libraries
import numpy as np
import demes, demesdraw
import msprime as ms

# Population N change model
def size_change(Ns,time_period,yaml_filename=None,plot=True,plot_filename=None,time_units="years",generation_time=29):
    
    if time_units == "generations":
        generation_time=1

    m = demes.Builder(time_units=time_units,generation_time=generation_time)
    current_time=np.array(time_period).sum()

    epochs = []
    for N,time in zip(Ns,time_period):
        epoch = dict(start_size=N,end_time=current_time-time)
        current_time=current_time-time
        epochs.append(epoch)

    m.add_deme("Deme",epochs=epochs)

    # Resolve model
    graph = m.resolve()

    # Check demographic events
    print(epochs)
    
    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        p = demesdraw.tubes(graph, log_time=True, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph 

# Cake model function
def cake_model(Ns,splits,alpha1,alpha2,time_period_merge,time_period_splits,migration_rate=1e-4,yaml_filename=None,plot=True,plot_filename=None,time_units="generations",generation_time=1):
    """
    alpha1 = proportions of previous N for split
    alpha2 = propottions of contribution to mergerger per splitted deme
    """
    # Check arguments
    assert len(splits) == len(alpha1) == len(time_period_splits), "Proportions and time period list must be the same length as number of split events."
    assert len(splits)+1 == len(time_period_merge), "Time period merge list must be the same length as number of split events + 1."
    if time_units == "generations":
        generation_time=1
        
    merge_events = len(splits)+1
    assert len(Ns) == merge_events, "Length of Ns list must be equal to number of split events + 1"

    total_time = np.sum(np.array(time_period_merge).sum() + np.array(time_period_splits).sum())

    m = demes.Builder(time_units=time_units,generation_time=generation_time)

    # Add first Ancestor
    m.add_deme("Ancestral",epochs=[dict(start_size=Ns[0], end_time=total_time-time_period_merge[0])])
    current_time = total_time-time_period_merge[0]
    split_b = 1
    previous = ["Ancestral"]

    #
    event = 0 
    split_i = 0
    while current_time > 0:
        if split_b:
            pops = []
            assert splits[split_i] == len(alpha1[split_i]), "Proportions list must have the same length as the number of splits"
            for pop_i,proportion in zip(np.arange(splits[split_i]),alpha1[split_i]):
                name="Split_" + str(event) + str(pop_i)
                m.add_deme(name,ancestors=previous,start_time=current_time,epochs=[dict(start_size=Ns[event]*proportion,end_time=current_time-time_period_splits[event])])
                pops.append(name)
            previous = pops
            if migration_rate > 0:
                m.add_migration(demes = pops, rate = migration_rate)
            current_time = current_time-time_period_splits[event]
            split_b = 0
            event = event + 1
        else: 
            assert len(previous) == len(alpha1[split_i]),"Length of ancestors is not equal to proportions"
            proportion = alpha2[split_i] if len(previous) > 1 else [1]
            name="Merge_" + str(event)
            m.add_deme(name,ancestors=previous,proportions=proportion,start_time=current_time,epochs=[dict(start_size=Ns[event],end_time=current_time-time_period_merge[event])])
            previous = [name]
            current_time = current_time-time_period_merge[event]
            split_b = 1
            split_i = split_i + 1

    # Resolve model
    graph = m.resolve()

    # Check demographic events
    print(graph.discrete_demographic_events()['splits'])

    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        p = demesdraw.tubes(graph, log_time=True, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph

# Load model from yaml
def load_yaml(yaml_file):
    return demes.load(yaml_file)

# https://tskit.dev/msprime/docs/stable/demography.html#sec-demography-numerical-trajectories

def get_iicr(demes_model,pop,T=None):
    """
    Returns two arrays: the Coalescence rate and the Inferred Inverse Coalescence Rate (Popsize)
    """
    m = ms.Demography.from_demes(demes_model)
    debug = m.debug()
    if np.sum(T) == None:
        T = np.concatenate([
            np.linspace(0, 1000, 2001),
            np.linspace(1000, 1.5e4, 401)[1:]
        ])
    R, _ = debug.coalescence_rate_trajectory(T, {pop: 2})
    inversed_R = 1/(2*R)

    return R,inversed_R,T

# Population N change model
def get_N_times_from_iicr(iicr,T):
    """
    Takes an the inverse inferred rate of coalescence and the time points and returns the times at which size changes.
    """
    previous_N = np.flip(iicr)[0]
    Ns = []
    times = []
    for N,time in zip(np.flip(iicr),np.flip(T)):
        if int(N) != previous_N:
            Ns.append(previous_N)
            times.append(time)
            previous_N = int(N)
    Ns.append(N)        
    times.append(0)
    return Ns,times

def size_change_from_iicr(iicr,T,yaml_filename=None,plot=True,plot_filename=None,time_units="years",generation_time=29):
    """
    Takes a vector of Ns in form of iicr and T times to create a model
    """
    if time_units == "generations":
        generation_time=1

    Ns,times = get_N_times_from_iicr(iicr,T)
    
    m = demes.Builder(time_units=time_units,generation_time=generation_time)

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
        p = demesdraw.tubes(graph, log_time=True, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph  
    