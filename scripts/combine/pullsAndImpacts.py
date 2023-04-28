import uproot
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import argparse
import dash
from dash import dcc
import dash_daq as daq
from dash import html
from dash.dependencies import Input, Output
from utilities import input_tools, logging
import os

def writeOutput(fig, outfile, extensions=[], postfix=None):
    name, ext = os.path.splitext(outfile)
    if ext not in extensions:
        extensions.append(ext)

    if postfix:
        name += f"_{postfix}"

    for ext in extensions:
        if ext[0] != ".":
            ext = "."+ext
        func = "write_html" if ext == ".html" else "write_image"

        output = name+ext
        logger.debug(f"Write output file {output}")

        getattr(fig, func)(output)

def plotImpacts(df, title, pulls=False, pullrange=[-5,5], POI='Wmass', normalize=False, oneSidedImpacts=False):
    poi_type = POI.split("_")[-1]

    if poi_type == "Wmass":
        impact_title="Impact on mass (MeV)"
    elif poi_type == "mu":
        impact_title=r"$\Delta \mu$"
    elif poi_type == "pmaskedexp":
        impact_title=r"$\delta N$" if normalize else r"$\Delta N$"
    elif poi_type == "pmaskedexpnorm":
        impact_title=r"$\delta (\mathrm{d} \sigma / \sigma)$" if normalize else r"$\Delta(\mathrm{d} \sigma / \sigma)$"

    fig = make_subplots(rows=1,cols=2 if pulls else 1,
            horizontal_spacing=0.1, shared_yaxes=True)

    max_pull = np.max(df["pull"])
    default_pr = pullrange == [-5, 5]
    if max_pull == 0 and default_pr:
        pullrange = [-1.5, 1.5]
    elif default_pr:
        r = np.max([1.5, max_pull])
        pullrange = [-1*r, r]
    
    ndisplay = len(df)
    fig.add_trace(
        go.Bar(
            x=df['impact' if not oneSidedImpacts else 'absimpact'],
            y=df['label'],
            marker_color=df['impact_color'] if oneSidedImpacts else '#377eb8',
            texttemplate="%{x:0.2f}",
			textposition="outside",
            textfont_size=12,
			textangle=0,
            orientation='h',
            name="impacts_down",
        ),
        row=1,col=1,
    )
    if not oneSidedImpacts:
        fig.add_trace(
            go.Bar(
                x=-1*df['impact'],
                y=df['label'],
                marker_color='#e41a1c',
                name="impacts_up",
                orientation='h',
            ),
            row=1,col=1,
    )
    if pulls:
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
                name="pulls",
            ),
            row=1,col=2,
    )
    impact_range = np.ceil(df['impact'].max())
    impact_spacing = 2 if pulls else 3
    if impact_range % impact_spacing:
        impact_range += impact_spacing - (impact_range % impact_spacing)
    tick_spacing = impact_range/impact_spacing
    if pulls and oneSidedImpacts:
        tick_spacing /= 2.
    fig.update_layout(barmode='relative')
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis_title=impact_title,
        title={
            'text': POI,
            'y': 1.-1/float(ndisplay),
            'x': 0.7,
            'xanchor': 'center',
            'yanchor': 'top'},
        margin=dict(l=20,r=20,t=50,b=20),
        xaxis=dict(range=[-impact_range if not oneSidedImpacts else -impact_range/20, impact_range],
                showgrid=True, gridwidth=1,gridcolor='Gray', griddash='dash',
                zeroline=True, zerolinewidth=2, zerolinecolor='Gray',
                tickmode='linear',
                tickangle=0,
                tick0=0.,
                side='top',
                dtick=tick_spacing,
            ),
        yaxis=dict(range=[-1, ndisplay]),
        showlegend=False,
        height=100*(ndisplay<100)+ndisplay*20.5,width=800,
    )
    if pulls:
        fig.update_layout(
            xaxis2=dict(range=pullrange,
                    showgrid=True, gridwidth=2,gridcolor='LightBlue',
                    zeroline=True, zerolinewidth=4, zerolinecolor='Gray',
                    tickmode='linear',
                    tick0=0.,
                    dtick=0.5 if pullrange[1]-pullrange[0] > 2.5 else 0.25,
                    side='top',
                ),
            xaxis2_title="pull+constraint",
            yaxis2=dict(range=[-1, ndisplay]),
            yaxis2_visible=False,
        )
    return fig

def readFitInfoFromFile(rf,filename, group=False, sort=None, ascending=True, stat=0.0, POI='Wmass', normalize=False):    
    treename = "fitresults"
    tree = rf[treename]
    # TODO: Make add_total configurable
    add_total = group
    impacts, labels, _ = input_tools.readImpacts(rf, group, sort=sort, add_total=add_total, stat=stat, POI=POI, normalize=normalize)
    # TODO: Make configurable
    if True:
        impacts = impacts*100
    
    pulls = np.zeros_like(impacts)
    constraints = np.zeros_like(impacts)
    if not group:
        import ROOT
        rtfile = ROOT.TFile.Open(filename)
        rtree = rtfile.Get(treename)
        rtree.GetEntry(0)
        for i, label in enumerate(labels):
            if not hasattr(rtree, label):
                logger.warning(f"Failed to find syst {label} in tree")
                continue
                
            pulls[i] = getattr(rtree, label)
            constraints[i] = getattr(rtree, label+"_err")
    
    df = pd.DataFrame(np.array((pulls, impacts, constraints), dtype=np.float64).T, columns=["pull", "impact", "constraint"])
    df['label'] = labels
    df['absimpact'] = np.abs(df['impact'])
    if not group:
        df.drop(df.loc[df['label']=='massShift100MeV'].index, inplace=True)
    colors = np.full(len(df), '#377eb8')
    if not group:
        colors[df['impact'] > 0.] = '#e41a1c'
    df['impact_color'] = colors

    if sort:
        df = df.sort_values(by=sort, ascending=ascending)

    return df

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--inputFile", type=str, required=True,
        help="fitresults output ROOT file from combinetf")
    parser.add_argument("-s", "--sort", default="absimpact", type=str, choices=["label", "pull", "constraint", "absimpact"], help="Sort mode for nuisances")
    parser.add_argument("--stat", default=0.0, type=float, help="Overwrite stat. uncertainty with this value")
    parser.add_argument("-d", "--sortDescending", dest='ascending', action='store_false', help="Sort mode for nuisances")
    parser.add_argument("-g", "--group", action='store_true', help="Show impacts of groups")
    parser.add_argument("--absolute", action='store_true', help="Not normalize impacts on cross sections and event numbers.")
    parser.add_argument("--debug", action='store_true', help="Print debug output")
    parser.add_argument("--oneSidedImpacts", action='store_true', help="Make impacts one-sided")
    parsers = parser.add_subparsers(dest='mode')
    interactive = parsers.add_parser("interactive", help="Launch and interactive dash session")
    interactive.add_argument("-i", "--interface", default="localhost", help="The network interface to bind to.")
    output = parsers.add_parser("output", help="Produce plots as output (not interactive)")
    output.add_argument("-o", "--outputFile", default="test.html", type=str, help="Output file (extension specifies if html or pdf/png)")
    output.add_argument("--otherExtensions", default=[], type=str, nargs="*", help="Additional output file types to write")
    output.add_argument("-n", "--num", type=int, help="Number of nuisances to plot")
    output.add_argument("--noPulls", action='store_true', help="Don't show pulls (not defined for groups)")
    output.add_argument("-t", "--title", type=str, default="", help="Title of output plot")
    
    return parser.parse_args()

app = dash.Dash(__name__)

@app.callback(
    Output("scatter-plot", "figure"), 
    [Input("maxShow", "value")],
    [Input("sortBy", "value")],
    [Input("sortDescending", "on")],
    [Input("filterLabels", "value")],
    [Input("groups", "on")],
)
def draw_figure(maxShow, sortBy, sortDescending, filterLabels, groups, oneSidedImpacts=False):
    df = groupsdataframe if groups else dataframe
    df = df[0 if not maxShow else -1*maxShow:].copy()
    df = df.sort_values(by=sortBy, ascending=sortDescending)
    if filterLabels:
        filt = np.full(len(df["label"]), False)
        for label in filterLabels.split(","):
            filt = filt | (df["label"].str.find(label.strip()) >= 0)
        df = df[filt]
    return plotImpacts(df, title="Contribution to uncertainty in mW", pulls=not groups, pullrange=[-5,5], oneSidedImpacts=oneSidedImpacts)

dataframe = pd.DataFrame()
groupsdataframe = pd.DataFrame()

        
def producePlots(rtfile,args, POI='Wmass', normalize=False):

    groupsdataframe = readFitInfoFromFile(rtfile,args.inputFile, True, sort=args.sort, ascending=args.ascending, stat=args.stat/100., POI=POI, normalize=normalize)
    if not (args.group and args.mode == 'output'):
        dataframe = readFitInfoFromFile(rtfile,args.inputFile, False, sort=args.sort, ascending=args.ascending, stat=args.stat/100., POI=POI, normalize=normalize)
    if args.mode == "interactive":
        app.layout = html.Div([
                dcc.Input(
                    id="maxShow", type="number", placeholder="maxShow",
                    min=0, max=10000, step=1,
                ),
                dcc.Input(
                    id="filterLabels", type="text", placeholder="filter labels (comma-separated list)",
                    style={
                        'width': '25%'
                    },
                ),
                html.Br(),
                html.Label('Sort by'),
                dcc.Dropdown(
                    id='sortBy',
                    options=[{'label' : v, "value" : v.lower()} for v in ["Impact", "Pull", "Constraint", "Label"]
                    ],
                    placeholder='select sort criteria...',
                    style={
                        'width': '50%'
                    },
                    value=args.sort
                ),
                daq.BooleanSwitch('sortDescending', label="Decreasing order", labelPosition="top", on=True),
                daq.BooleanSwitch('groups', label="Show nuisance groups", labelPosition="top", on=False),
                dcc.Graph(id="scatter-plot", style={'width': '100%', 'height': '100%'}),
            ],
            style={
                "width": "100%",
                "height": "100%",
                "display": "inline-block",
                "padding-top": "10px",
                "padding-left": "1px",
                "overflow": "hidden"
            },
        )

        app.run_server(debug=True, port=3389, host=args.interface)
    elif args.mode == 'output':
        df = dataframe if not args.group else groupsdataframe
        df = df.sort_values(by=args.sort, ascending=args.ascending)
        if args.num and args.num < df.size:
            df = df[-args.num:].sort_values(by=args.sort, ascending=args.ascending)
        fig = plotImpacts(df, title=args.title, pulls=not args.noPulls and not args.group, pullrange=[-5,5], POI=POI, normalize=not args.absolute, oneSidedImpacts=args.oneSidedImpacts)

        outputfilename = args.outputFile
        outputfilename
        writeOutput(fig, args.outputFile, args.otherExtensions, postfix=POI)
    else:
        raise ValueError("Must select mode 'interactive' or 'output'")


if __name__ == '__main__':
    args = parseArgs()

    logger = logging.setup_logger("pullsAndImpacts", 4 if args.debug else 3)

    rtfile = uproot.open(args.inputFile)
    POIs = input_tools.getPOInames(rtfile)
    for POI in POIs:
        producePlots(rtfile,args,POI)
    
    if not (POIs[0]=="Wmass" and len(POIs) == 1):
        # masked channel
        for POI in input_tools.getPOInames(rtfile, poi_type="pmaskedexp"):
            producePlots(rtfile,args,POI, normalize=not args.absolute)

        # masked channel normalized
        for POI in input_tools.getPOInames(rtfile, poi_type="pmaskedexpnorm"):
            producePlots(rtfile,args,POI, normalize=not args.absolute)