# import ROOT
# from math import pi

# can = ROOT.TCanvas('can', 'can')

# f2d = ROOT.TFile('matchingScheme2/skim_all_Matched.root')
# f3d = ROOT.TFile('matchingScheme2/skim_all_Matched_with_deltaE.root')

# treename = 'MatchedTree'

# t2d = f2d.Get(treename)
# t3d = f3d.Get(treename)

# h2d_chi2 = ROOT.TH1D('h2d_chi2', '#chi^2', 10**6, 0, 10**6)
# h3d_chi2 = ROOT.TH1D('h3d_chi2', '#chi^2', 10**6, 0, 10**6)
# t2d.Draw('Metric >> h2d_chi2')
# t3d.Draw('Metric >> h3d_chi2')

# can.SetLogx()
# can.SetLogy()

# h2d_chi2.Draw()
# h3d_chi2.Draw('same')


# def draw_comparison(var, name, nbins=400, xmin=-pi, xmax=-pi, vmin=-0.5, vmax=0.5):
#     h2d = ROOT.TH1D(f'h2d_{name}', '#chi^{2}(#theta,#phi)', nbins, xmin, xmax)
#     h3d = ROOT.TH1D(f'h3d_{name}', '#chi^{2}(#theta,#phi,E)', nbins, xmin, xmax)
#     h2d.SetLineColor(ROOT.kBlack)
#     h3d.SetLineColor(ROOT.kRed)
#     h2d.SetLineWidth(1)
#     h3d.SetLineWidth(1)
#     can.SetLogx(0)
#     can.SetLogy(1)
#     t2d.Draw(f'{var} >> h2d_{name}', 'IsGen && IsReco')
#     t3d.Draw(f'{var} >> h3d_{name}', 'IsGen && IsReco')
#     h2d.GetXaxis().SetRangeUser(vmin, vmax)
#     h3d.GetXaxis().SetRangeUser(vmin, vmax)
#     h2d.Scale(1/h2d.Integral())
#     h3d.Scale(1/h3d.Integral())
#     h3d.Draw('same')
#     h2d.Draw()
#     h3d.Draw('same')
#     can.BuildLegend()
#     h2d.SetTitle(f';{name};probability')
#     can.SaveAs(f'{var}.png')

# draw_comparison('DeltaPhi', '#Delta#phi', 50, -0.2, 0.2, -0.2, 0.2)
# draw_comparison('DeltaTheta', '#Delta#theta', 50, -0.2, 0.2, -0.2, 0.2)

# # draw_comparison('Metric', '#chi^{2}', 100, 0, 100, 0, 100)

import ROOT
from math import pi

ROOT.EnableImplicitMT(40)
ROOT.gStyle.SetOptStat(0)
can = ROOT.TCanvas("can", "can")

# Use RDataFrame
treename = "MatchedTree"
df2d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched.root")
df3d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched_with_deltaE.root")

# # Plot the chi^2 distributions
# h2d_chi2 = df2d.Histo1D(("h2d_chi2", "#chi^{2}", 10**6, 0, 10**6), "Metric")
# h3d_chi2 = df3d.Histo1D(("h3d_chi2", "#chi^{2}", 10**6, 0, 10**6), "Metric")

# can.SetLogx()
# can.SetLogy()
# h2d_chi2.Draw()
# h3d_chi2.Draw("same")

# def rdefine(df, column_name, expression):
#     """
#     Define a new column in the RDataFrame if it doesn't exist.
#     If it already exists, redefine it.
#     """
#     if column_name in df.GetColumnNames():
#         return df
#         # Drop the column if it exists, then define again
#         # df = df.Redefine(column_name, expression)
#     else:
#         df = df.Define(column_name, expression)
#     return df


# def rdefine(self, column_name, expression):
#     """
#     Define a new column in the RDataFrame if it doesn't exist.
#     If it already exists, redefine it.
#     """
#     if column_name in self.GetColumnNames():
#         # return self.Redefine(column_name, expression)
#         return self
#     else:
#         return self.Define(column_name, expression)

# # Attach the method to ROOT.RDF.RInterface
# ROOT.RDF.RInterface.__rdefine__ = rdefine
# ROOT.RDF.RInterface.rdefine = rdefine
# ROOT.RDataFrame.rdefine = rdefine



# rdf_iface = type(ROOT.RDataFrame("dummy", ROOT.TTree()))
# rdf_iface.rdefine = rdefine


def patch_rdefine(iface=None):
    if not iface:
        dummy_df = ROOT.RDataFrame(1)  # minimal dummy frame with 1 entry
        iface = type(dummy_df)

    def rdefine(self, column_name, expression):
        if column_name in self.GetColumnNames():
            # return self.Redefine(column_name, expression)
            return self
        else:
            return self.Define(column_name, expression)

    iface.rdefine = rdefine

patch_rdefine()


# def draw_comparison(var, name, nbins=400, xmin=-pi, xmax=pi, vmin=-0.5, vmax=0.5):
#     can.SetLogx(0)
#     can.SetLogy(1)
#     mvar = f'matched_{var}'

#     df2d_sel = df2d.Define(mvar, f"{var}[IsGen && IsReco]")
#     df3d_sel = df3d.Define(mvar, f"{var}[IsGen && IsReco]")

#     h2d = df2d_sel.Histo1D((f"h2d_{name}", "#chi^{2}(#theta,#phi)", nbins, xmin, xmax), mvar)
#     h3d = df3d_sel.Histo1D((f"h3d_{name}", "#chi^{2}(#theta,#phi,E)", nbins, xmin, xmax), mvar)

#     h2d.SetLineColor(ROOT.kBlack)
#     h3d.SetLineColor(ROOT.kRed)
#     h2d.SetLineWidth(1)
#     h3d.SetLineWidth(1)

#     # Normalize histograms
#     h2d.Scale(1.0 / h2d.Integral())
#     h3d.Scale(1.0 / h3d.Integral())

#     h2d.GetXaxis().SetRangeUser(vmin, vmax)
#     h2d.Draw()
#     h3d.Draw("same")
#     can.BuildLegend()
#     h2d.SetTitle(f";{name};probability")
#     can.SaveAs(f"{var}.png")

# import ROOT
# from math import pi

import ROOT
from math import pi

def draw_sel(var,
             name,
             df,
             selections,
             nbins=400,
             xmin=-pi,
             xmax=pi,
             vmin=-0.5,
             vmax=0.5,
             logy=True,
             normalize=True,
             outname=None):
    """
    var         : branch name to plot (e.g. "theta")
    name        : axis title (e.g. "θ rec-gen")
    df          : ROOT.RDataFrame to draw from
    selections  : list of tuples (label, cutString),
                  cutString is e.g. "IsGen && IsReco"
    nbins,xmin,xmax : histogram binning
    vmin,vmax   : x-range to zoom into
    logy        : if True, set log-scale on y
    outname     : filename to save canvas (defaults to f"{var}.png")
    """
    can.cd()
    can.Clear()
    if logy:
        can.SetLogy(1)
    else:
        can.SetLogy(0)

    # prepare a color palette
    palette = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
               ROOT.kMagenta, ROOT.kOrange, ROOT.kCyan+1]

    mvar = f"matched_{var}"

    histos = []
    for i, (label, cut) in enumerate(selections):
        # define filtered variable array
        df_sel = df.Define(mvar, f"{var}[{cut}]")

        # make histogram
        h = df_sel.Histo1D(
            (f"h_{label}", f"{name}; {name}; probability", nbins, xmin, xmax),
            mvar
        ).GetValue()

        # normalize
        integral = h.Integral()
        if normalize and integral > 0:
            h.Scale(1.0 / integral)

        # style
        col = palette[i % len(palette)]
        # h.SetLineColor(col)
        # h.SetLineWidth(2)
        h.SetFillStyle(1002)
        h.SetFillColorAlpha(col, 0.3)
        histos.append((label, h))

    maxh = max([h[1].GetMaximum() for h in histos])
    # draw overlays
    first = True
    for label, h in histos:
        h.SetMaximum(maxh * 1.2)
        h.SetMinimum(100)
        h.GetXaxis().SetRangeUser(vmin, vmax)
        if first:
            h.Draw("hist")
            first = False
        else:
            h.Draw("hist same")

    # legend
    leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
    leg.SetBorderSize(0)
    for label, h in histos:
        leg.AddEntry(h, label, "f")
    leg.Draw()

    # save
    # outname = outname or f"{var}.png"
    outname = outname or f"{var}.pdf"
    can.SaveAs(outname)

def draw_comparison(var, name, nbins=400, xmin=-pi, xmax=pi, vmin=-0.5, vmax=0.5,
                    sel=""):
    # Create canvas and pads
    can.cd()
    pad1 = ROOT.TPad("pad1", "Main", 0.0, 0.3, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "Ratio", 0.0, 0.0, 1.0, 0.3)
    pad1.SetBottomMargin(0.02)  # little overlap
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.35)
    pad1.Draw()
    pad1.SetLogy(1)
    pad2.Draw()

    # Define and fill histos
    mvar = f"matched_{var}"
    # selection = "(!IsGen && IsReco)"
    selection = "(IsGen && IsReco)"
    if sel != "":
        selection += f" && {sel}"
    df2d_sel = df2d.Define(mvar, f"{var}[{selection}]")
    df3d_sel = df3d.Define(mvar, f"{var}[{selection}]")

    h2d = df2d_sel.Histo1D((f"h2d_{name}", name, nbins, xmin, xmax), mvar).GetValue()
    h3d = df3d_sel.Histo1D((f"h3d_{name}", name, nbins, xmin, xmax), mvar).GetValue()

    # Normalize
    for h in (h2d, h3d):
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)

    # Upper pad: main distributions
    pad1.cd()
    h2d.SetLineColor(ROOT.kBlack)
    h3d.SetLineColor(ROOT.kRed)
    h2d.SetLineWidth(1)
    h3d.SetLineWidth(1)
    h2d.GetXaxis().SetRangeUser(vmin, vmax)
    h2d.SetTitle(f"{sel};{name};Probability")
    h2d.Draw("hist")
    h3d.Draw("hist same")
    leg = ROOT.TLegend(0.65, 0.7, 0.9, 0.85)
    leg.AddEntry(h2d, "Match w/ p_{T}", "l")
    leg.AddEntry(h3d, "Match w/ E", "l")
    leg.SetBorderSize(0)
    leg.Draw()

    # Lower pad: ratio
    pad2.cd()
    h_ratio = h3d.Clone(f"ratio_{name}")
    h_ratio.Divide(h2d)
    h_ratio.SetTitle(f";{name};3D/2D")
    h_ratio.GetYaxis().SetRangeUser(0.7, 1.3)
    h_ratio.GetXaxis().SetLabelSize(0.12)
    h_ratio.GetYaxis().SetLabelSize(0.12)
    h_ratio.GetYaxis().SetTitleSize(0.14)
    h_ratio.GetYaxis().SetTitleOffset(0.5)
    h_ratio.GetXaxis().SetTitleSize(0.14)
    # draw a line at ratio=1
    line = ROOT.TLine(vmin, 1.0, vmax, 1.0)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineWidth(1)

    h_ratio.Draw("ep")
    line.Draw("same")

    # Save
    can.SaveAs(f"{var}_comparison.png")



# draw_comparison("DeltaPhi", "#Delta#phi")
# draw_comparison("DeltaPhi", "#Delta#phi", 100, -0.5, 0.5, -0.5, 0.5)
# draw_comparison("DeltaPhi", "#Delta#phi", 100, -0.1, 0.1, -0.1, 0.1)
# draw_comparison("DeltaPhi", "#Delta#phi", 100, -0.05, 0.05, -0.05, 0.05)
# draw_comparison("DeltaTheta", "#Delta#theta", 100, -0.2, 0.2, -0.2, 0.2)
# draw_comparison("DeltaE", "#DeltaE", 100, -1, 5, -1, 5)


def extend(df):
    scaleFactorTheta = 0.5 * 2.8  # sigma(rz)   = 28 µm
    scaleFactorPhi   = 0.5 * 2.3  # sigma(rphi) = 23 µm

    # sigmaDelta = 25e-6 + 95e-6 / meanP
    # sigmaTheta = sigmaDelta / 0.06 * scaleFactorTheta
    # sigmaPhi   = sigmaDelta / 0.06 * scaleFactorPhi

    # df = df.Define("genvec", "ROOT::VecOps::Construct<ROOT::Math::PxPyPzEVector>(\
    # GenX,\
    # GenY,\
    # GenZ,\
    # GenE)").Define("genp",
    #                "return ROOT::VecOps::Map(genvec, [](const auto& v){ return v.mag(); })").Define(
    #     'sigma_phi', f'(25e-6 + 95e-6 / genp) / 0.06 * {scaleFactorPhi}').Define(
    #     'sigma_theta', f'(25e-6 + 95e-6 / genp) / 0.06 * {scaleFactorTheta}').Define(
    #         'chi2_phi', 'DeltaPhi / sigma_phi').Define(
    #         'chi2_theta', 'DeltaTheta / sigma_theta')
    df = df.\
        Define("genp", """
        return ROOT::VecOps::Map(GenX, GenY, GenZ, []
        (const auto& x, const auto& y, const auto& z){
        return TMath::Sqrt(x*x + y*y + z*z); })
        """).\
        Define("genpt","""
        return ROOT::VecOps::Map(GenX, GenY, []
        (const auto& x, const auto& y){
        return TMath::Sqrt(x*x + y*y); })
        """).\
        Define('sigma_phi', f'(25e-6 + 95e-6 / genp) / 0.06 * {scaleFactorPhi}').\
        Define('sigma_theta', f'(25e-6 + 95e-6 / genp) / 0.06 * {scaleFactorTheta}').\
        Define('chi_phi', 'DeltaPhi / sigma_phi').\
        Define('chi_theta', 'DeltaTheta / sigma_theta')

    df = df.\
        Define('sigma_pt', 'sqrt(pow((6e-4 * genpt), 2) + pow(5e-3, 2)) * genpt').\
        Define('pt_over_40', 'genpt > 40').\
        Define('chi_phi_high', 'chi_phi[pt_over_40]')
    return df

df2d = extend(df2d)
df3d = extend(df3d)


# draw_comparison("genp", "p", 100, 0, 100, 0, 100)
# draw_comparison("sigma_theta", "#sigma#theta", 100, 0, 1e-2, 0, 1e-2)

# draw_comparison("sigma_phi", "#sigma#phi", 100, 0, 1e-2, 0, 1e-2)
# draw_comparison("sigma_pt", "#sigma p_{T}", 100, 0, 1e-2, 0, 1e-2)
# draw_comparison("sigma_pt", "#sigma p_{T}", 100, 0, 1e-1, 0, 1e-1, 'genpt > 10')
# draw_comparison("sigma_pt", "#sigma p_{T}", 100, 0, 10, 0, 10, 'genpt > 50')

# draw_comparison("chi2_phi", "#chi^{2}#phi", 100, -5e-1, 5e-1, -5e-1, 5e-1)

# draw_comparison("chi2_phi", "#chi^{2}#phi", 100, -3, 3, -3, 3)
# draw_comparison("chi2_phi", "#chi^{2}#phi", 100, -1e2, 1e2, -1e2, 1e2)
# draw_comparison("chi2_phi", "#chi^{2}#phi", 100, -3, 3, -3, 3, sel="genp < 5")
# draw_comparison("chi2_theta", "#chi^{2}#theta", 100, -5, 5, -5, 5)
# draw_comparison("chi2_pt", "#chi^{2}p_{T}", 100, -3, 3, -3, 3, sel="genp > 20")
# draw_comparison("DeltaPhi", "#Delta#phi", 100, -1e-2, 1e-2, -1e-2, 1e-2, sel="genp > 60")


# draw_comparison("chi_pt_high", "#chi p_{T}", 100, -1e-0, 1e-0, -1e-0, 1e-0)

# draw_comparison("genpt", "gen pt", 50, 0, 20, 0, 20)
# draw_comparison("chi_pt", "#chi p_{T}", 100, -3, 3, -3, 3, sel="genpt > 20")
# draw_comparison("chi_pt", "#chi p_{T}", 100, -3, 3, -3, 3, sel="")
# draw_comparison("genp", "gen p", 50, 0, 2, 0, 2)
# draw_comparison("GenX", "gen X", 50, 0, 2, 0, 2)
# draw_comparison("GenY", "gen y", 50, 0, 2, 0, 2)
# draw_comparison("GenZ", "gen z", 50, 0, 2, 0, 2)

# df2d.Filter('(Sum(chi_phi_high) > 0)').Display(['genpt', 'chi_phi_high', 'chi_phi', 'GenX', "RecoX", "RecoEta", "GenEta"], 5, 100).Print()

# draw_comparison("chi_phi", "#chi #phi", 100, -3, 3, -3, 3)
# draw_comparison("chi_phi", "#chi #phi", 100, -3, 3, -3, 3, 'genp > 1 && genp < 2')
# draw_comparison("chi_phi", "#chi #phi", 100, -3, 3, -3, 3, "pt_over_40")
# draw_comparison("chi_phi", "#chi #phi", 100, -10, 10, -10, 10, "pt_over_40")
# draw_comparison("chi_phi", "#chi #phi", 100, -10, 10, -10, 10, "pt_over_40")

# draw_comparison("DeltaPhi", "#Delta#phi", 100, -0.5, 0.5, -0.5, 0.5)
# draw_comparison("DeltaPhi", "#Delta#phi", 100, -0.05, 0.05, -0.05, 0.05)

ptlist = [0.2, 2, 5, 10, 50, 90]
# ptlist = [2, 5, 10, 50, 90]
ptlist = [0.2, 2, 5, 10, 50]
ptsel = [(f'Gen pT #in ({ptlow, pthigh})', f'IsGen && IsReco && genpt > {ptlow} && genpt < {pthigh}') for (ptlow, pthigh) in zip(ptlist, ptlist[1:])]

# plist = [1, 2, 3, 4, 8]
# plist = [2, 4, 10]
plist = [10, 20, 30, 40, 50]
psel = [(f'Gen p #in ({plow, phigh})', f'IsGen && IsReco && genp > {plow} && genp < {phigh}') for (plow, phigh) in zip(plist, plist[1:])]
limit = 20e-3
draw_sel("DeltaPhi", "#Delta #phi", df2d, psel, 101,  -limit, limit, -limit, limit)

draw_sel("MinDeltaPhi", "#Delta #phi", df2d, psel, 50,  -limit, limit, -limit, limit)
# draw_sel("sigma_pt", "#sigma p_{T}", df2d, ptsel, 100, 0, 0.2, 0, 0.2)
draw_sel("chi_phi", "#chi #phi", df2d, psel, 100,  -5, 5, -5, 5, normalize=False)
draw_sel("sigma_phi", "#sigma #phi", df2d, psel, 100, 5e-4, 1e-3, 5e-4, 1e-3)



# # 2D plot
# dfsel = df2d.Define('dphi_high', f'DeltaPhi[{psel[-1][1]}]')
# dfsel = dfsel.Define('chiphi_high', f'chi_phi[{psel[-1][1]}]')
# h2d = dfsel.Histo2D(('f2d', 'f2d;dphi;chiphi', 100, -1e-2, 1e-2, 100, -5, 5), 'dphi_high', 'chiphi_high').GetValue()
# can.SetLogx(0)
# can.SetLogy(0)
# can.SetLogz(1)
# h2d.Draw('col'); can.Update()
# can.SaveAs('precision.png')

df2d.
