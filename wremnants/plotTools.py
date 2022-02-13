import mplhep as hep
import matplotlib.pyplot as plt
from wremnants import boostHistHelpers as hh
from wremnants import histselections as sel
hep.style.use(hep.style.ROOT)

def makePlotWithRatio(obs, histInfo, stackedProcs, unstacked=None, xlabel="", ylabel="Events/bin", 
                rrange=[0.9, 1.1], scale=8.5e6, action=None):
    width=3 if "unrolled" in obs else 1
    fig = plt.figure(figsize=(8*width,8))
    ax1 = fig.add_subplot(4, 1, (1, 3)) 
    ax2 = fig.add_subplot(4, 1, 4) 
    
    op = lambda x: x.project(obs) 
    if obs == "unrolled":
        op = sel.unrolledHist

    stack = [op(histInfo[k]["hist"]) for k in stackedProcs if histInfo[k]["hist"]]
    colors = [histInfo[k]["color"] for k in stackedProcs if histInfo[k]["hist"]]
    labels = [histInfo[k]["label"] for k in stackedProcs if histInfo[k]["hist"]]
            
    hep.histplot(
        stack,
        histtype="fill",
        color=colors,
        label=labels,
        stack=True,
        ax=ax1
    )
    
    if unstacked:
        unstack = op(histInfo[unstacked]["hist"])
        hep.histplot(
            unstack,
            yerr=True, 
            histtype="errorbar",
            color=histInfo[unstacked]["color"],
            label=histInfo[unstacked]["label"],
            ax=ax1,
        )
        hep.histplot(
            hh.divideHists(unstack, sum(stack)),
            histtype="errorbar",
            color=histInfo[unstacked]["color"],
            label=histInfo[unstacked]["label"],
            yerr=True,
            ax=ax2
        )
        
    ax1.set_xlabel("")
    ax2.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_xticklabels([])
    xrange = [stack[0].axes[0].edges[0], stack[0].axes[0].edges[stack[0].axes[0].size-1]]
    ax1.set_xlim(xrange)
    ax2.set_xlim(xrange)
    ax2.set_ylabel("data/pred.", fontsize=22)
    ax2.set_ylim(rrange)
    ax1.set_ylim([0, scale])
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(reversed(handles), reversed(labels), prop={'size' : 20*(0.7 if width == 1 else 1.3)}, ncol=2, loc='upper right')
    return fig
