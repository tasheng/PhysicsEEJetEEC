import ROOT
from math import pi

# f2d = ROOT.TFile('matchingScheme2/skim_all_Matched.root')
f2d = ROOT.TFile('matchingScheme2/skim_all_Matched_no_cutoff.root')
# f3d = ROOT.TFile('matchingScheme2/skim_all_Matched_with_deltaE.root')

treename = 'MatchedTree'

t2d = f2d.Get(treename)
# t3d = f3d.Get(treename)

ROOT.EnableImplicitMT(8)
ROOT.gStyle.SetOptStat(0)
can = ROOT.TCanvas("can", "can", 1200, 1000)

# Use RDataFrame
treename = "MatchedTree"
# df2d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched.root")
# df3d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched_with_deltaE.root")
# df2d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched.root")
df2d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched_no_cutoff.root")
df3d = ROOT.RDataFrame(treename, "matchingScheme2/skim_all_Matched_no_cutoff.root")

# # # Plot the chi^2 distributions
# h2d_chi2 = df2d.Histo1D(("h2d_chi2", "#chi^{2}", 10**6, 0, 10**6), "Metric")
# h3d_chi2 = df3d.Histo1D(("h3d_chi2", "#chi^{2}", 10**6, 0, 10**6), "Metric")

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
             outname=None,
             array=True):
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
        if array:
            df_sel = df.Define(mvar, f"{var}[{cut}]")
        else:
            df_sel = df.Filter(f'Any({cut})').Define(mvar, var)

        # make histogram
        ytitle = "probability" if normalize else "count"
        h = df_sel.Histo1D(
            (f"h_{label}", f"{name}; {name}; {ytitle}", nbins, xmin, xmax),
            mvar
        ).GetValue()

        # normalize
        integral = h.Integral()
        if normalize and integral > 0:
            h.Scale(1.0 / integral)

        # style
        col = palette[i % len(palette)]
        h.SetLineColor(col)
        h.SetLineWidth(2)
        h.SetFillStyle(1002)
        h.SetFillColorAlpha(col, 0.3)
        histos.append((label, h))

    maxh = max([h[1].GetMaximum() for h in histos])
    # draw overlays
    first = True
    for label, h in histos:
        h.SetMaximum(maxh * 1.2)
        # h.SetMinimum(100)
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
    h_ratio = h2d.Clone(f"ratio_{name}")
    h_ratio.Divide(h3d)
    h_ratio.SetTitle(f";{name};3D/2D")
    h_ratio.GetYaxis().SetRangeUser(0.95, 1.02)
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
# dphi_limit = 5e-2
# draw_comparison("DeltaPhi", "#Delta#phi", 61, -dphi_limit, dphi_limit, -dphi_limit, dphi_limit)
# draw_comparison("DeltaTheta", "#Delta#theta", 41, -dphi_limit, dphi_limit, -dphi_limit, dphi_limit)
# draw_comparison("DeltaE", "#DeltaE", 100, -1, 5, -1, 5)
# draw_comparison("DeltaPhi", "#Delta#phi", 101, -dphi_limit, dphi_limit, -dphi_limit, dphi_limit)
# draw_comparison("DeltaTheta", "#Delta#theta", 101, -dphi_limit, dphi_limit, -dphi_limit, dphi_limit)


def extend(df):
    scaleFactorTheta = 2.8  # sigma(rz)   = 28 µm
    scaleFactorPhi   = 2.3  # sigma(rphi) = 23 µm

    # sigmaDelta = 25e-6 + 95e-6 / meanP
    # sigmaTheta = sigmaDelta / 0.06 * scaleFactorTheta
    # sigmaPhi   = sigmaDelta / 0.06 * scaleFactorPhi

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
        Define('pt_over_40', 'genpt > 40')

    df = df.\
        Define('recopt_reci', '1 / sqrt(RecoX * RecoX + RecoY * RecoY)').\
        Define('genpt_reci', '1/genpt').\
        Define('DeltaPt_Reci', 'recopt_reci - genpt_reci').\
        Define('sigma_pt_reci', 'sqrt(pow((6e-4 * genpt), 2) + pow(5e-3, 2))')

    return df

df2d = extend(df2d)
df3d = extend(df3d)

draw_comparison("Metric", "#chi", 101, 0, 10000, 0, 10000, sel='genp > 30')




# plist = [1, 10, 20, 30, 40, 50]
# plist = [0.5, 1, 5, 10, 20, 30, 100]
plist = [.2, .3, .5, .7, 1.0]
# plist = [.5, .7, 1.0]
# plist = [40, 50, 60, 70]
# plist = [1, 2, 3, 5, 10, 40]
# plist = [40, 80]
psel = [(f'Gen p #in ({plow, phigh})', f'IsGen && IsReco && genp > {plow} && genp < {phigh}') for (plow, phigh) in zip(plist, plist[1:])]
limit = 20e-3
draw_sel("DeltaPhi", "#Delta#phi", df2d, psel, 41,  -limit, limit, -limit, limit)
# draw_sel("MinDeltaPhi", "Min. #Delta#phi", df2d, psel, 41,  -limit, limit, -limit, limit)
# limit = 5e-3
# draw_sel("DeltaTheta", "#Delta#theta", df2d, psel, 201,  -limit, limit, -limit, limit)

# ptlist = [0.2, 2, 5, 10, 20, 30]
ptlist = [2, 5, 10, 20, 30, 40, 50]
ptsel = [(f'Gen pT #in {ptlow, pthigh}', f'IsGen && IsReco && genpt > {ptlow} && genpt < {pthigh}') for (ptlow, pthigh) in zip(ptlist, ptlist[1:])]
# ptsel = [(f'Gen 1/pT #in 1/{pthigh}, 1/{ptlow}', f'IsGen && IsReco && genpt_reci > {1/pthigh} && genpt_reci < {1/ptlow}') for (ptlow, pthigh) in zip(ptlist, ptlist[1:])]
ptlimit = 1e1
draw_sel("DeltaPt", "#Deltap_{T}", df2d, ptsel[::-1], 101,  -ptlimit, ptlimit, -ptlimit, ptlimit)
limit = 1e-2
draw_sel("DeltaPt_Reci", "#Delta#frac{1}{p_{T}}", df2d, ptsel, 101,  -limit, limit, -limit, limit)

# draw_sel("sigma_phi", "#sigma #phi", df2d, psel, 100, 5e-4, 2e-3, 5e-4, 2e-3, normalize=False)
# draw_sel("sigma_theta", "#sigma #theta", df2d, psel, 100, 5e-4, 2e-3, 5e-4, 2e-3, normalize=False)


####################
psel = [('outside', '(IsGen && IsReco && genp > 40 && MinDeltaPhi > 0.002)'),
        ('inside', '(IsGen && IsReco && genp > 40 && abs(MinDeltaPhi) < 0.001)'),
        ]
nlimit = 60
draw_sel("NReco", "# of reco particles", df2d, psel, nlimit,  0, nlimit, 0, nlimit, array=False)
draw_sel("NReco", "#Delta #phi", df2d, psel, nlimit,  0, nlimit, 0, nlimit, array=False)


# draw_sel("MinDeltaPhi", "#Delta #phi", df2d, psel, 50,  -limit, limit, -limit, limit)
draw_sel("sigma_pt", "#sigma p_{T}", df2d, ptsel, 100, 0, 0.1, 0, 0.1, normalize=False)

draw_sel("chi_phi", "#chi #phi", df2d, psel, 101,  -5, 5, -5, 5, normalize=False)
draw_sel("chi_phi", "#chi #phi", df2d, psel, 101,  -5, 5, -5, 5)
# draw_sel("sigma_phi", "#sigma #phi", df2d, psel, 100, 5e-4, 2e-3, 5e-4, 2e-3, normalize=False)


# df2d.Define('genp_high', 'genp > 40').Filter('Sum(genp_high) > 0').Filter('DeltaPhi[genp_high] >0.01').Count()
dfsel = df2d.Define('genp_high', 'genp > 40')
# dfsel.Filter('Any(genp_high && IsGen && IsReco) & Any(DeltaPhi[genp_high && IsGen && IsReco] > 0.01)').Display(['EventID']).Print()
# h2d = dfsel.Histo2D(('f2d', 'f2d;dphi;chiphi', 100, -1e-2, 1e-2, 100, -5, 5), 'dphi_high', 'chiphi_high').GetValue()

# # 2D plot
# dfsel = df2d.Define('dphi_high', f'DeltaPhi[{psel[-1][1]}]')
# dfsel = dfsel.Define('chiphi_high', f'chi_phi[{psel[-1][1]}]')
# h2d = dfsel.Histo2D(('f2d', 'f2d;dphi;chiphi', 100, -1e-2, 1e-2, 100, -5, 5), 'dphi_high', 'chiphi_high').GetValue()
# can.SetLogx(0)
# can.SetLogy(0)
# can.SetLogz(1)
# h2d.Draw('col'); can.Update()
# can.SaveAs('precision.png')


# df2d.Filter('(Sum(chi_phi_high) > 0)').Display(['genpt', 'chi_phi_high', 'chi_phi', 'GenX', "RecoX", "RecoEta", "GenEta"], 5, 100).Print()

def fit_gaussians_to_selections(var,
                                 name,
                                 df,
                                 selections,
                                 nbins,
                                 xmin=-pi,
                                 xmax=pi,
                                 fit_frac=0.5,
                                 array=True):
    """
    Fit Gaussians to histograms created from selections.

    Parameters:
    - var        : branch name to plot
    - name       : axis title
    - df         : ROOT.RDataFrame
    - selections : list of (label, cutString)
    - nbins      : list of nbins (one per selection)
    - xmin, xmax : histogram range
    - fit_frac   : fraction of histogram entries to use for fit (centered)
    - array      : whether var is an array or not

    Returns:
    - List of RooFitResult* (one per selection)
    """
    from ROOT import RooRealVar, RooDataHist, RooGaussian, RooFit, RooPlot, TCanvas
    import re

    def sanitize_label(label):
        # Replace special characters for use in RooFit object names
        return re.sub(r"[^a-zA-Z0-9_\-]", "_", label)

    mvar = f"matched_{var}"
    x = RooRealVar("x", name, xmin, xmax)

    results = []
    canvas_height = 300 * len(selections)
    can = TCanvas("can", "fit canvas", 800, canvas_height)
    can.Divide(1, len(selections))

    for i, (label, cut) in enumerate(selections):
        sanitized_label = sanitize_label(label)
        this_nbins = nbins[i] if isinstance(nbins, list) else nbins

        # define variable
        if array:
            df_sel = df.Define(mvar, f"{var}[{cut}]")
        else:
            df_sel = df.Filter(f"Any({cut})").Define(mvar, var)

        h = df_sel.Histo1D(
            (f"h_{sanitized_label}", f"{name}; {name}; count", this_nbins, xmin, xmax),
            mvar
        ).GetValue()

        # determine fit range based on histogram content
        total_entries = h.Integral()
        cumulative = 0
        fit_bins = []
        center_bin = h.GetXaxis().FindBin(0)

        for delta in range(this_nbins // 2):
            left = center_bin - delta
            right = center_bin + delta
            if left < 1 or right > this_nbins:
                break
            cumulative = h.Integral(left, right)
            if cumulative / total_entries >= fit_frac:
                fit_low = h.GetXaxis().GetBinLowEdge(left)
                fit_high = h.GetXaxis().GetBinUpEdge(right)
                break
        else:
            fit_low = xmin
            fit_high = xmax

        x.setRange("fitRange", fit_low, fit_high)

        # convert to RooDataHist
        dh = RooDataHist(f"dh_{sanitized_label}", f"dh_{label}", [x], h)

        # define gaussian
        mean = RooRealVar(f"mean_{sanitized_label}", "mean", 0, xmin, xmax)
        sigma = RooRealVar(f"sigma_{sanitized_label}", "sigma", (xmax - xmin)/10, 0.001, xmax - xmin)
        gauss = RooGaussian(f"gauss_{sanitized_label}", "gaussian PDF", x, mean, sigma)

        # fit
        result = gauss.fitTo(dh, RooFit.Save(True), RooFit.Range("fitRange"))
        results.append((label, result))

        # plot
        frame = x.frame()
        dh.plotOn(frame)
        gauss.plotOn(frame)
        gauss.paramOn(frame, RooFit.Parameters((mean, sigma)), RooFit.Format("NELU", RooFit.AutoPrecision(2)), RooFit.Layout(0.6, 0.9, 0.85))
        frame.SetTitle(label)
        frame.GetXaxis().SetTitleSize(0.06)
        frame.GetYaxis().SetTitleSize(0.06)
        frame.GetXaxis().SetLabelSize(0.05)
        frame.GetYaxis().SetLabelSize(0.05)

        can.cd(i+1)
        frame.Draw()

    can.Update()
    can.SaveAs(f"fit_{var}.pdf")

    return results


def fit_gaussians_to_selections(var,
                                 name,
                                 df,
                                 selections,
                                 nbins=400,
                                 xmin=-pi,
                                 xmax=pi,
                                 fit_frac=0.5,
                                 array=True):
    """
    Fit Gaussians to histograms created from selections.

    Parameters:
    - var        : branch name to plot
    - name       : axis title
    - df         : ROOT.RDataFrame
    - selections : list of (label, cutString)
    - nbins, xmin, xmax : histogram binning
    - fit_frac   : fraction of histogram entries to use for fit (centered)
    - array      : whether var is an array or not

    Returns:
    - List of RooFitResult* (one per selection)
    """
    from ROOT import RooRealVar, RooDataHist, RooGaussian, RooFit, RooPlot, TCanvas
    import re

    def sanitize_label(label):
        # Replace special characters for use in RooFit object names
        return re.sub(r"[^a-zA-Z0-9_]", "_", label)

    mvar = f"matched_{var}"
    x = RooRealVar("x", name, xmin, xmax)

    results = []
    canvas_height = 300 * len(selections)
    can = TCanvas("fit_can", "fit canvas", 800, canvas_height)
    can.Divide(1, len(selections))

    for i, (label, cut) in enumerate(selections):
        sanitized_label = sanitize_label(label)

        # define variable
        if array:
            df_sel = df.Define(mvar, f"{var}[{cut}]")
        else:
            df_sel = df.Filter(f"Any({cut})").Define(mvar, var)

        h = df_sel.Histo1D(
            (f"h_{sanitized_label}", f"{name}; {name}; count", nbins, xmin, xmax),
            mvar
        ).GetValue()

        # determine fit range based on histogram content
        total_entries = h.Integral()
        cumulative = 0
        fit_bins = []
        center_bin = h.GetXaxis().FindBin(0)

        for delta in range(nbins // 2):
            left = center_bin - delta
            right = center_bin + delta
            if left < 1 or right > nbins:
                break
            cumulative = h.Integral(left, right)
            if cumulative / total_entries >= fit_frac:
                fit_low = h.GetXaxis().GetBinLowEdge(left)
                fit_high = h.GetXaxis().GetBinUpEdge(right)
                break
        else:
            fit_low = xmin
            fit_high = xmax

        x.setRange("fitRange", fit_low, fit_high)

        # convert to RooDataHist
        dh = RooDataHist(f"dh_{sanitized_label}", f"dh_{label}", [x], h)

        # define gaussian
        mean = RooRealVar(f"mean_{sanitized_label}", "mean", 0, xmin, xmax)
        sigma = RooRealVar(f"sigma_{sanitized_label}", "sigma", (xmax - xmin)/10, 1e-4, xmax - xmin)
        gauss = RooGaussian(f"gauss_{sanitized_label}", "gaussian PDF", x, mean, sigma)

        # fit
        result = gauss.fitTo(dh, RooFit.Save(True), RooFit.Range("fitRange"))
        results.append(result)

        # plot
        frame = x.frame()
        dh.plotOn(frame)
        gauss.plotOn(frame)
        gauss.paramOn(frame, RooFit.Parameters((mean, sigma)),
                      RooFit.Format("NELU", RooFit.AutoPrecision(2)),
                      RooFit.Layout(0.6, 0.9, 0.85))
        frame.SetTitle(label)
        frame.GetXaxis().SetTitleSize(0.06)
        frame.GetYaxis().SetTitleSize(0.06)
        frame.GetXaxis().SetLabelSize(0.05)
        frame.GetYaxis().SetLabelSize(0.05)

        can.cd(i+1)
        frame.Draw()

    can.Update()
    can.SaveAs(f"fit_{var}.pdf")

    return results

results = fit_gaussians_to_selections(
    var="DeltaPt_Reci",
    name="#Delta#frac{1}{p_{T}}",
    df=df2d,
    selections=ptsel,
    nbins=51,
    xmin=-5e-3,
    xmax=5e-3,
    fit_frac=0.9
)

# print fit results
for res in results:
    res.Print()

def plot_sigmas(selections, results, xtitle="Selection Value"):
    """
    Plot the sigma values with errors from RooFitResult list returned by fit_gaussians_to_selections.
    """
    from ROOT import TGraphErrors, TCanvas
    from array import array
    import re

    def extract_numeric(label):
        numbers = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", label)
        return float(numbers[0]) if numbers else float('nan')

    xvals = []
    yvals = []
    yerrs = []

    for (label, _), result in zip(selections, results):
        x = extract_numeric(label)
        par = result.floatParsFinal().find("sigma_" + re.sub(r"[^a-zA-Z0-9_\-]", "_", label))
        sigma = par.getVal()
        sigma_err = par.getError()
        xvals.append(x)
        yvals.append(sigma)
        yerrs.append(sigma_err)

    zeros = [0] * len(xvals)
    graph = TGraphErrors(len(xvals), array('d', xvals), array('d', yvals), array('d', zeros), array('d', yerrs))
    graph.SetTitle(f"Sigma vs {xtitle};{xtitle};#sigma")
    c = TCanvas("c_sigmas", "Sigma Plot", 800, 600)
    graph.SetMarkerStyle(20)
    graph.Draw("APL")
    c.SaveAs("sigmas_vs_selection.pdf")

plot_sigmas(ptsel, results, xtitle="#Delta#frac{1}{p_{T}}")


# plist = [0.2, 0.3, 0.5, 1]
plist = [0.2, 0.3, 0.5, 1, 5, 10, 20, 30, 100]
psel = [(f'Gen p #in ({plow, phigh})', f'IsGen && IsReco && genp > {plow} && genp < {phigh}') for (plow, phigh) in zip(plist, plist[1:])]
meanp = [(plist[i] + plist[i+1])/2 for i in range(len(plist) - 1)]
results = fit_gaussians_to_selections(
    var="DeltaPhi",
    name="#Delta#phi",
    df=df2d,
    selections=psel,
    nbins=[int(1 + 0.04 // (1e-3 / p)) for p in meanp],
    xmin=-20e-3,
    xmax=20e-3,
    fit_frac=0.9
)
