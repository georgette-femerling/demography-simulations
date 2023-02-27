# Import libraries
import moments, moments.LD
import numpy as np
from collections import defaultdict


def Dstat(graph,sampled_demes,rhos,normalize = True, norm_pop_idx = 0,sample_times = None):
     if sample_times is None:
          sample_times = np.zeros(len(sampled_demes))
     y = moments.Demes.LD(graph, sampled_demes=sampled_demes, sample_times=sample_times,rho = rhos)
     if normalize:
          sigma = moments.LD.Inference.sigmaD2(y,normalization=norm_pop_idx)
     else:
          sigma = y
     return sigma


def get_LD_from_sliced_demes(sliced_dict,rhos = np.logspace(-2, 2, 21),normalize = True):
    LD_dictionary = defaultdict(list)
    for sliced in sliced_dict:
        alive = [deme.name for deme in sliced_dict[sliced].demes if deme.end_time == 0 ]
        sizes = [deme.epochs[0].start_size for deme in sliced_dict[sliced].demes if deme.end_time == 0 ]
        norm_idx = sizes.index(max(sizes))
        sigma = Dstat(sliced_dict[sliced], sampled_demes=alive, rhos = rhos, normalize = normalize, norm_pop_idx=norm_idx)
        for deme,i in zip(alive,range(len(alive))):
            DD = 'DD_'+str(i)+"_"+str(i)
            Dz = 'Dz_'+str(i)+"_"+str(i)+"_"+str(i)
            pi = 'pi2_'+str(i)+"_"+str(i)+"_"+str(i)+"_"+str(i)
            sigmapop = sigma.LD()[:,[sigma.names()[0].index(stat) for stat in [DD,Dz,pi]]]
            LD_dictionary[deme].append(sigmapop)
    return LD_dictionary

def calculate_signal_Dz(LDpop1,LDpop2):
    '''Calculates the Mean log difference between the two decays of Dz'''
    x = LDpop1[:,1]
    y = LDpop2[:,1]
    log_diff = np.diff(np.log(x) - np.log(y)) # calculate difference
    mld= np.mean(log_diff)
    return mld

def calculate_signal_D2(LDpop1,LDpop2):
    '''Calculates the Mean log difference between the two decays of D2'''
    x = LDpop1[:,0]
    y = LDpop2[:,0]
    log_diff = np.diff(np.log(x) - np.log(y)) # calculate difference
    mld= np.mean(log_diff)
    return mld