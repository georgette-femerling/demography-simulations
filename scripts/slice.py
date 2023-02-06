# Import libraries
import moments, moments.LD
from collections import defaultdict
import numpy as np
import demes

def slice_timepoints(g,time_points,yaml_filename=None):
    sliced_dict = defaultdict()
    for time in time_points:
        slicedg = moments.Demes.DemesUtil.slice(g, time)
        sliced_dict[time] = slicedg
        # option to save to yaml
        if yaml_filename != None:
            demes.dump(slicedg, yaml_filename+"_"+str(int(time))+".tmp", format='yaml', simplified=True)
    
    return sliced_dict

def sample_timepoints(graph, time_points_perdeme = 3, anc_end_time=100000):
    times = []
    times_dic = defaultdict()
    for deme in graph.demes:
        deme = deme.name
        time = np.linspace(graph[deme].end_time+1,graph[deme].start_time-1,time_points_perdeme,dtype=int) if graph[deme].start_time != float('inf') else np.linspace(graph[deme].end_time+1,graph[deme].end_time+anc_end_time,time_points_perdeme,dtype=int)
        times = np.concatenate([times,time])
        times_dic[deme] = time
    return np.unique(np.flip(times)),times_dic