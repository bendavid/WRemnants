import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import ROOT
import sys
import math

def plotPullsAndImpacts(df, fOut, sort="impact", ascending = True, nmax=-1):

    df = df.sort_values(by=sort, ascending=ascending)

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.1, shared_yaxes=True)

    max_pull = np.max(df["pull"])
    if max_pull == 0:
        pullrange = [-1.1, 1.1]
    else:
        r = np.max([1.1, max_pull])
        pullrange = [-1*r, r]

    ndisplay = len(df["impact"])
    fig.add_trace(
        go.Scatter(
            x=np.zeros(ndisplay),
            y=df['label'],
            mode="markers",
            marker=dict(color='black', size=8,),
            error_x=dict(
                array=df['impact'],
                color="black",
                thickness=1.5,
                width=5,
            ),
            name="impacts",
        ),
        row=1,col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=df['pull'],
            y=df['label'],
            mode="markers",
            marker=dict(color='black', size=8,),
            error_x=dict(
                array=df['constraint'],
                color="black",
                thickness=1.5,
                width=5,
            ),
            name="impacts",
        ),
        row=1,col=2,
    )
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis_title="Impact on average signal strength (%)",
        title={
            'text': "",
            'y':.999 if ndisplay > 100 else 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        margin=dict(l=20,r=20,t=25,b=20, pad=0),
        #xaxis=dict(range=[-25, 25],
        xaxis=dict(range=[-20, 20], # -3, 3
                showgrid=True, gridwidth=2,gridcolor='LightPink',
                zeroline=True, zerolinewidth=4, zerolinecolor='Gray',
                tickmode='linear',
                tick0=0.,
                dtick=2
            ),
        yaxis=dict(range=[-1, ndisplay]),
        showlegend=False,
        height=25+ndisplay*25,width=1000,
    )

    fig.update_layout(
        xaxis2=dict(range=pullrange,
            showgrid=True, gridwidth=2,gridcolor='LightBlue',
            zeroline=True, zerolinewidth=4, zerolinecolor='Gray',
            tickmode='linear',
            tick0=0.,
            dtick=1 if pullrange[1]-pullrange[0] > 2.5 else 0.5,
        ),
        xaxis2_title="Pull + constraints",
        yaxis2=dict(range=[-1, ndisplay]),
        yaxis2_visible=False,
    )

    fig.update_xaxes(side="top")
    fig.write_html(fOut)


def impacts_avg(filename, fOut):

    fIn = ROOT.TFile(filename)
    tree = fIn.Get("fitresults")
    tree.GetEntry(0)

    name = "nuisance_group_impact_mu"
    name = "nuisance_group_impact_pmaskedexpnorm"
    impact_hist = fIn.Get(name)
    xLabels = np.array([impact_hist.GetXaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsX()+1)]) # POI
    yLabels = np.array([impact_hist.GetYaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsY()+1)]) # nuisances
    yLabels_ = np.append(yLabels, "Total")

    scale = 100
    rounded = 2
    mus = [getattr(tree, "%s" % xLabel) for xLabel in xLabels]
    print(mus)
    impacts = np.array([round(np.mean(np.array([impact_hist.GetBinContent(k+1, i+1)/mus[k] for k,xLabel in enumerate(xLabels)]))*scale, rounded) for i,yLabel in enumerate(yLabels)])

    tot = 0.
    for k,xLabel in enumerate(xLabels): tot += getattr(tree, "%s_err" % xLabel)/mus[k]
    tot /= len(xLabels)
    impacts = np.append(impacts, round(tot*scale, rounded))



    df = pd.DataFrame(np.array((np.abs(impacts)), dtype=np.float64).T, columns=["impact"])
    df.insert(0, "label", yLabels_)
    df = df.sort_values(by="impact", ascending=True)

    fig = make_subplots(rows=1,cols=1, vertical_spacing = 0.1)
    fig.add_trace(
        go.Scatter(
            x=np.zeros(len(yLabels_)),
            y=df['label'],
            mode="markers",
            marker=dict(color='black', size=8,),
            error_x=dict(
                array=df['impact'],
                color="black",
                thickness=1.5,
                width=5,
            ),
            name="impacts",
        ),
        row=1,col=1,
    )


    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis_title = "Average impact on signal strength (%)",
        title={
            'text': "",
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        margin=dict(l=20,r=20,t=25,b=20, pad=0),
        xaxis=dict(range=[-20, 20],
                showgrid=True, gridwidth=2,gridcolor='LightPink',
                zeroline=True, zerolinewidth=4, zerolinecolor='Gray',
                tickmode='linear',
                tick0=0.,
                dtick=2
        ),
        yaxis=dict(range=[-1, len(yLabels_)]),
        showlegend=False,
        height=25+len(yLabels_)*35,width=600,
    )

    fig.update_xaxes(side="top")
    fig.write_html(fOut)


def impacts_mu(filename, fOut):

    fIn = ROOT.TFile(filename)
    tree = fIn.Get("fitresults")
    tree.GetEntry(0)

    name = "nuisance_group_impact_nois"
    impact_hist = fIn.Get(name) # nuisance_group_impact_nois
    xLabels = np.array([impact_hist.GetXaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsX()+1)]) # POI
    yLabels = np.array([impact_hist.GetYaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsY()+1)]) # nuisances
    yLabels_ = np.append(yLabels, "Total")
    
    xLabels_unsorted = xLabels
    xLabels = np.sort(xLabels)

    scale = 100
    rounded = 2
    dfs = []
    for k,xLabel in enumerate(xLabels):

        idx = int(np.where(xLabels_unsorted == xLabel)[0][0])
        impacts_ = np.array([round(impact_hist.GetBinContent(idx+1, i+1)*scale, rounded) for i,yLabel in enumerate(yLabels)])

        # add total impact
        impacts_ = np.append(impacts_, round(getattr(tree, "%s_err" % xLabel)*scale, rounded))

        df = pd.DataFrame(np.array((np.abs(impacts_)), dtype=np.float64).T, columns=["impact"])
        df.insert(0, "label", yLabels_)
        df = df.sort_values(by="impact", ascending=True)
        dfs.append(df)

    cols = 2
    rows = math.ceil(1.*xLabels.size/cols)
    fig = make_subplots(rows=rows,cols=cols, horizontal_spacing=0.2, vertical_spacing = 0.1, subplot_titles=xLabels)

    for k,xLabel in enumerate(xLabels):

        row = int(k/cols)+1
        col = k%2+1

        fig.add_trace(
            go.Scatter(
                x=np.zeros(len(yLabels_)),
                y=dfs[k]['label'],
                mode="markers",
                marker=dict(color='black', size=8,),
                error_x=dict(
                    array=dfs[k]['impact'],
                    color="black",
                    thickness=1.5,
                    width=5,
                ),
                name="",
            ),
            row=row,col=col,
        )

        fig.update_xaxes(
            range=[-20, 20],
            showgrid=True, gridwidth=2,gridcolor='LightPink',
            zeroline=True, zerolinewidth=4, zerolinecolor='Gray',
            tickmode='linear',
            tick0=0.,
            dtick=2,
            row=row,col=col
        )

        fig.update_yaxes(
            range=[-1, len(yLabels_)],
            row=row,col=col
        )

    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis_title="",
        title={
            'text': "Impact on signal strength (%)",
            'y': 0.98,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},
        margin=dict(l=20,r=20,t=75,b=20, pad=0),
        showlegend=False,
        height=rows*len(yLabels_)*30,width=1200,
    )

    #fig.update_xaxes(side="top")
    fig.write_html(fOut)


def getPullsAndImpacts(filename):

    fIn = ROOT.TFile(filename)
    tree = fIn.Get("fitresults")
    tree.GetEntry(0)

    name = "nuisance_impact_nois"
    impact_hist = fIn.Get(name)
    xLabels = np.array([impact_hist.GetXaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsX()+1)]) # POI
    yLabels = np.array([impact_hist.GetYaxis().GetBinLabel(i) for i in range(1, impact_hist.GetNbinsY()+1)]) # nuisances

    # get impacts
    scale = 100
    rounded = 2
    impacts = np.array([round(np.mean(np.array([impact_hist.GetBinContent(k+1, i+1) for k,xLabel in enumerate(xLabels)]))*scale, rounded) for i,yLabel in enumerate(yLabels)])

    # get pulls and constraints
    pulls = np.zeros_like(impacts)
    constraints = np.zeros_like(impacts)
    for i, yLabel in enumerate(yLabels):
        pulls[i] = getattr(tree, yLabel)
        constraints[i] = getattr(tree, yLabel+"_err")

    df = pd.DataFrame(np.array((pulls, np.abs(impacts), constraints), dtype=np.float64).T, columns=["pull", "impact", "constraint"])
    df.insert(0, "label", yLabels)

    return df



if __name__ == '__main__':

    baseName = "ZMassWLike" # ZMassWLike_plus ZMassWLike_minus ZMassWLike
    fIn = f"ZMassWLike/{baseName}_output.root"
    
    
    baseName = "lowPU_highPU"
    fIn = f"CombineStudies/lowPU_wmass//{baseName}_output.root"
    
    
    baseDir = "/eos/user/j/jaeyserm/www/wmass/combine/lowPU_highPU/"

    impacts_mu(fIn, f"{baseDir}/{baseName}_impacts.html")

    df = getPullsAndImpacts(fIn)
    plotPullsAndImpacts(df, f"{baseDir}/{baseName}_pulls_impacts_sortImpact.html", "impact")
    plotPullsAndImpacts(df, f"{baseDir}/{baseName}_pulls_impacts_sortConstraint.html", "constraint", ascending=False)


    ###source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
    ###source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh