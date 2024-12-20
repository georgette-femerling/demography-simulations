# Import libraries
import matplotlib.pylab as plt
import numpy as np

def plot_iicr(iirc,T): 
    plt.plot(T, iirc)
    plt.xticks(fontsize= 12)
    plt.yticks(fontsize= 12)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("time ago (years)",fontsize = 14)
    plt.ylabel(r"IICR",fontsize = 14)
    
def plot_comparison(LDpop1,LDpop2,labels=["Original","Size Change"]):
    rhos = np.logspace(-2, 2, 21)
    
    # plot D2
    f = plt.figure(figsize=(10,3))
    ax = f.add_subplot(121)
    ax2 = f.add_subplot(122)

    ax.plot(rhos,LDpop1[:,0],label=labels[0])
    ax.plot(rhos,LDpop2[:,0],label=labels[1])
    ax.legend()
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel(r"$\sigma_d^2$")
    ax.set_xlabel(r"$\rho$")

    # plot DZ
    ax2.plot(rhos,LDpop1[:,1],label=labels[0])
    ax2.plot(rhos,LDpop2[:,1],label=labels[1])
    ax2.legend()
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylabel(r"$\sigma_{Dz}$")
    ax2.set_xlabel(r"$\rho$")

    plt.tight_layout()
    plt.show()
    
def plot_LD(LD_sigma, times_dic, ancestral, rhos = np.logspace(-2, 2, 21), plot_file = None, figsize = (10,20)):
    fig = plt.figure(constrained_layout=True, figsize = figsize)
    subfigs = fig.subfigures(nrows=len(LD_sigma.keys()), ncols=1)

    for pop,subfig in zip(LD_sigma,subfigs):
        subfig.suptitle(pop)
        (ax1, ax2, ax3) = subfig.subplots(nrows=1, ncols=3)
        for time_point in range(len(LD_sigma[pop])):
            ax1.plot(rhos, LD_sigma[pop][time_point][:, 0],label=str("tp_"+str(times_dic[pop][time_point])))
            ax2.plot(rhos, LD_sigma[pop][time_point][:, 1],label=str("tp_"+str(times_dic[pop][time_point])))
            ax3.plot(rhos, LD_sigma[pop][time_point][:, 2],label=str("tp_"+str(times_dic[pop][time_point])))
        
        ax1.plot(rhos, ancestral[:, 0],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax2.plot(rhos, ancestral[:, 1],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax3.plot(rhos, ancestral[:, 2],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
    
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax3.set_yscale("log")
        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")
        ax1.set_xlabel(r"$\rho$")
        ax2.set_xlabel(r"$\rho$")
        ax3.set_xlabel(r"$\rho$")
        ax1.set_ylabel(r"$\sigma_d^2$")
        ax2.set_ylabel(r"$\sigma_{Dz}$")
        ax3.set_ylabel(r"Pi2")

        ax1.legend()
        #ax2.legend()
        #ax3.legend()
    if plot_file != None :
        plt.savefig(plot_file,format='pdf',transparent = False)

def plot_LD_Linear(LD_sigma, times_dic, ancestral, rhos = np.logspace(-2, 2, 21), plot_file = None, figsize = (10,20)):
    fig = plt.figure(constrained_layout=True, figsize = figsize)
    subfigs = fig.subfigures(nrows=len(LD_sigma.keys()), ncols=1)

    for pop,subfig in zip(LD_sigma,subfigs):
        subfig.suptitle(pop)
        (ax1, ax2, ax3) = subfig.subplots(nrows=1, ncols=3)
        for time_point in range(len(LD_sigma[pop])):
            ax1.plot(rhos, LD_sigma[pop][time_point][:, 0]/ancestral[:, 0],label=str("tp_"+str(times_dic[pop][time_point])))
            ax2.plot(rhos, LD_sigma[pop][time_point][:, 1]/ancestral[:, 1],label=str("tp_"+str(times_dic[pop][time_point])))
            ax3.plot(rhos, LD_sigma[pop][time_point][:, 2]/ancestral[:, 2],label=str("tp_"+str(times_dic[pop][time_point])))
        
        ax1.plot(rhos, ancestral[:, 0],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax2.plot(rhos, ancestral[:, 1],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax3.plot(rhos, ancestral[:, 2],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
    
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax3.set_yscale("log")
        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")
        ax1.set_xlabel(r"$\rho$")
        ax2.set_xlabel(r"$\rho$")
        ax3.set_xlabel(r"$\rho$")
        ax1.set_ylabel(r"$\sigma_d^2/Ancestral$")
        ax2.set_ylabel(r"$\sigma_{Dz}/Ancestral$")
        ax3.set_ylabel(r"Pi2/Ancestral")

        ax1.legend()
        #ax2.legend()
        #ax3.legend()
    if plot_file != None :
        plt.savefig(plot_file,format='pdf',transparent = False)
