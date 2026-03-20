/*
===============================================================================
NPS_good_events.cpp
-------------------------------------------------------------------------------
Purpose:
  Standalone wrapper around NPS_selection_helper.h for Hall C NPS analysis.

Design:
  - Core event-selection logic lives in NPS_selection_helper.h
  - This wrapper adds:
      * CLI parsing
      * optional physics cuts on T
      * optional diagnostic plots
      * optional text report
      * optional reduced ROOT snapshots

Default behavior:
  - helicity mode: quartet_pm
  - current window enabled
  - I0 from BCM current mode unless --I0 is provided
  - window_frac = 0.15
  - rolling stability N = 30
  - no physics cuts
  - report enabled
  - snapshots enabled
  - plots disabled
===============================================================================
*/

#include <ROOT/RDataFrame.hxx>
#include <TString.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <regex>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TBox.h>

#include "NPS_selection_helper.h"

using namespace std;
using namespace ROOT;

static string run_seg_tag_from_input(const string& path) {
  const string base = gSystem->BaseName(path.c_str());

  // Matches e.g. nps_hms_coin_6786_7_1_-1.root
  static const std::regex re(R"(.*?_(\d+)_(\d+)_-?\d+_-?\d+\.root$)");
  std::smatch m;

  if (std::regex_match(base, m, re) && m.size() >= 3) {
    return Form("run%s_seg%s", m[1].str().c_str(), m[2].str().c_str());
  }

  // fallback if pattern is not matched
  return "runUNKNOWN_segUNKNOWN";
}

// =================================================================================
// Config
// =================================================================================
struct Config {
  // Input
  string input_root = "";

  // Helper selection settings
  string helicity_mode = "quartet_pm";   // quartet_pm | quartet_all | ignore
  bool use_current_window = true;
  double I0_fixed = -1.0;                // if > 0, use fixed I0; otherwise use peak
  double window_frac = 0.15;

  double junk_floor_uA = 2.5;
  double mean_current_calc_min = 2.5;
  double mean_current_min = 2.75;

  int stable_window_N = 30;              // N=1 effectively disables rolling stability
  double max_charge_frac_spread = 0.08;

  // Physics cuts on T (off by default)
  bool apply_physics_cuts = false;
  string physicsCuts =
    "(int)round(fEvtHdr.fEvtType) == 1 && "
    "H.dc.ntrack > 0 && "
    "H.hod.goodscinhit == 1 && "
    "H.hod.goodstarttime == 1 && "
    "H.hod.betanotrack > 0.5 && H.hod.betanotrack < 1.4 && "
    "H.cal.etottracknorm >= 0.7 && "
    "H.cal.etotnorm > 0.7 && "
    "H.cer.npeSum >= 2 && "
    "H.react.z > -10 && H.react.z < 10 && "
    "H.gtr.dp > -10 && H.gtr.dp < 10";
  
  // Wrapper outputs
  bool write_report = true;
  bool write_snapshots = true;
  bool make_plots = false;

  // Output
  string out_prefix = "nps_reduced";
  string out_dir = ".";

  // Branch names
  string branch_scaler_current = "H.BCM4A_Hel.scalerCurrent";
  string branch_scaler_charge  = "H.BCM4A_Hel.scalerCharge";
  string branch_helicity       = "actualHelicity";
  string branch_evcount        = "evcount";
  string branch_evNumber       = "evNumber";
  string branch_evnum_T        = "g.evnum";
};

static HelicityMode parse_helicity_mode(const string& s) {
  if (s == "quartet_pm")  return HelicityMode::QuartetPM;
  if (s == "quartet_all") return HelicityMode::QuartetAll;
  if (s == "ignore")      return HelicityMode::IgnoreHelicity;

  cerr << "Unknown helicity mode: " << s
       << " (use quartet_pm, quartet_all, or ignore)\n";
  exit(2);
}

static Config parse_args(int argc, char** argv) {
  Config cfg;

  for (int i = 1; i < argc; ++i) {
    string a = argv[i];

    auto require_value = [&](const string& key) -> string {
      if (i + 1 >= argc) {
        cerr << "Missing value after " << key << endl;
        exit(2);
      }
      return string(argv[++i]);
    };

    if (a == "--input") cfg.input_root = require_value(a);

    else if (a == "--helicity") cfg.helicity_mode = require_value(a);
    else if (a == "--no-current-window") cfg.use_current_window = false;
    else if (a == "--I0") cfg.I0_fixed = stod(require_value(a));
    else if (a == "--window-frac") cfg.window_frac = stod(require_value(a));

    else if (a == "--stable-window-N") cfg.stable_window_N = stoi(require_value(a));

    else if (a == "--physics-cuts") cfg.apply_physics_cuts = true;
    else if (a == "--plots") cfg.make_plots = true;

    else if (a == "--outdir") cfg.out_dir = require_value(a);
    else if (a == "--outprefix") cfg.out_prefix = require_value(a);

    else if (a == "--no-report") cfg.write_report = false;
    else if (a == "--no-snapshots") cfg.write_snapshots = false;

    else if (a == "--help") {
      cout <<
        "Usage:\n"
        "  root -l -b -q 'NPS_good_events.cpp+(\"--input\",\"/path/to/file.root\")'\n\n"
        "Options:\n"
        "  --input <path>                     Input ROOT file (required)\n"
        "  --helicity <quartet_pm|quartet_all|ignore>\n"
        "  --no-current-window                Disable current window cut\n"
        "  --I0 <uA>                          Use fixed I0; otherwise use peak automatically\n"
        "  --window-frac <f>                  Fractional current window around I0 (default 0.15)\n"
        "  --stable-window-N <int>            Rolling stability window N (default 30; N=1 ~ off)\n"
        "  --physics-cuts                     Apply combined physics cuts to T\n"
        "  --plots                            Make all diagnostic plots\n"
        "  --outdir <dir>\n"
        "  --outprefix <str>\n"
        "  --no-report\n"
        "  --no-snapshots\n";
      exit(0);
    }
    else {
      cerr << "Unknown arg: " << a << "\n";
      exit(2);
    }
  }

  return cfg;
}

static string build_input_path(const Config& cfg) {
  if (cfg.input_root.empty()) {
    cerr << "Error: need --input <path>\n";
    exit(2);
  }
  return cfg.input_root;
}

static string join_path(const string& dir, const string& file) {
  if (dir.empty() || dir == ".") return file;
  if (!dir.empty() && dir.back() == '/') return dir + file;
  return dir + "/" + file;
}

static SelectionSettings make_selection_settings(const Config& cfg) {
  SelectionSettings sel;

  sel.helicity_mode = parse_helicity_mode(cfg.helicity_mode);

  sel.junk_floor_uA = cfg.junk_floor_uA;
  sel.mean_current_calc_min = cfg.mean_current_calc_min;
  sel.mean_current_min = cfg.mean_current_min;

  sel.use_current_window = cfg.use_current_window;
  sel.i0_mode = (cfg.I0_fixed > 0.0) ? I0Mode::Fixed : I0Mode::Peak;
  sel.i0_fixed_uA = cfg.I0_fixed;
  sel.window_frac = cfg.window_frac;

  sel.stable_window_N = cfg.stable_window_N;
  sel.max_charge_frac_spread = cfg.max_charge_frac_spread;

  sel.branch_scaler_current = cfg.branch_scaler_current;
  sel.branch_scaler_charge  = cfg.branch_scaler_charge;
  sel.branch_helicity       = cfg.branch_helicity;
  sel.branch_evcount        = cfg.branch_evcount;
  sel.branch_evNumber       = cfg.branch_evNumber;
  sel.branch_evnum_T        = cfg.branch_evnum_T;

  return sel;
}


// =================================================================================
// Rolling-spread diagnostics helpers (wrapper only)
// =================================================================================
struct PreselectedScalerData {
  vector<double> evcount;
  vector<double> charge;
};

struct RollingSpreadData {
  vector<double> evcount_eval;
  vector<double> frac_spread;
  vector<double> pass_mask;   // 1 = pass, 0 = fail
};

static PreselectedScalerData build_preselected_scaler_data(const string& rootfile,
                                                           const SelectionSettings& sel)
{
  RDataFrame dH("TSHelH", rootfile);

  auto d_pre = dH.Filter(Form("%s > %f",
                              sel.branch_scaler_current.c_str(),
                              sel.junk_floor_uA));

  const string hel_expr = build_helicity_filter_expr(sel);
  if (!hel_expr.empty()) d_pre = d_pre.Filter(hel_expr);

  if (sel.use_current_window) {
    double I0 = -1.0;
    if (sel.i0_mode == I0Mode::Fixed) {
      I0 = sel.i0_fixed_uA;
    } else {
      auto hI = d_pre.Histo1D(
        {"hI_tmp_diag","BCM current;I [uA];counts", 120, 0.0, 60.0},
        sel.branch_scaler_current
      );
      I0 = hI->GetXaxis()->GetBinCenter(hI->GetMaximumBin());
    }

    const double Imin = (1.0 - sel.window_frac) * I0;
    const double Imax = (1.0 + sel.window_frac) * I0;

    d_pre = d_pre.Filter(
      Form("(%s >= %f) && (%s <= %f)",
           sel.branch_scaler_current.c_str(), Imin,
           sel.branch_scaler_current.c_str(), Imax)
    );
  }

  auto ev_ptr = d_pre.Take<double>(sel.branch_evcount);
  auto q_ptr  = d_pre.Take<double>(sel.branch_scaler_charge);

  PreselectedScalerData out;
  out.evcount.assign(ev_ptr->begin(), ev_ptr->end());
  out.charge.assign(q_ptr->begin(), q_ptr->end());

  vector<size_t> idx(out.evcount.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;

  std::sort(idx.begin(), idx.end(),
            [&](size_t a, size_t b) { return out.evcount[a] < out.evcount[b]; });

  vector<double> ev_sorted, q_sorted;
  ev_sorted.reserve(idx.size());
  q_sorted.reserve(idx.size());

  for (size_t k : idx) {
    ev_sorted.push_back(out.evcount[k]);
    q_sorted.push_back(out.charge[k]);
  }

  out.evcount.swap(ev_sorted);
  out.charge.swap(q_sorted);
  return out;
}

static RollingSpreadData compute_rolling_spread_data(const PreselectedScalerData& pre,
                                                     int N)
{
  RollingSpreadData out;

  if (N < 1) N = 1;
  const int n = (int)pre.charge.size();
  if (n <= 0) return out;

  for (int i = N - 1; i < n; ++i) {
    double qmin = pre.charge[i];
    double qmax = pre.charge[i];
    double qsum = 0.0;

    for (int j = i - N + 1; j <= i; ++j) {
      const double q = pre.charge[j];
      qsum += q;
      if (q < qmin) qmin = q;
      if (q > qmax) qmax = q;
    }

    const double qmean = qsum / double(N);
    const double frac  = (qmean > 0.0) ? (qmax - qmin) / qmean : 1e9;

    out.evcount_eval.push_back(pre.evcount[i]);
    out.frac_spread.push_back(frac);
  }

  return out;
}


// =================================================================================
// Optional diagnostics (wrapper only)
// =================================================================================
static void makeCurrentPanels(const string& rootfile,
                              const Config& cfg,
                              const SelectionSettings& sel,
                              const SelectionResult& pick)
{
  RDataFrame d_raw("TSHelH", rootfile);

  auto d_sel = d_raw.Filter(pick.evcount_cut);

  auto evMin = d_raw.Min<double>(sel.branch_evNumber).GetValue();
  auto evMax = d_raw.Max<double>(sel.branch_evNumber).GetValue();

  double yMin = 0.0;
  double yMax = std::max(30.0, d_raw.Max<double>(sel.branch_scaler_current).GetValue() * 1.05);
  int nBinsX = std::max(100, (int)((evMax - evMin) / 10.0));

  auto hRaw = d_raw.Histo2D(
    {"hRaw",
     "Raw TSHelH;evNumber;H.BCM4A_Hel.scalerCurrent [uA]",
     nBinsX, evMin, evMax, 200, yMin, yMax},
    sel.branch_evNumber,
    sel.branch_scaler_current
  );

  auto hSel = d_sel.Histo2D(
    {"hSel",
     "Selected TSHelH;evNumber;H.BCM4A_Hel.scalerCurrent [uA]",
     nBinsX, evMin, evMax, 200, yMin, yMax},
    sel.branch_evNumber,
    sel.branch_scaler_current
  );

  TCanvas c("c_current_panels","BCM current panels",1400,900);
  c.Divide(1,2);

  c.cd(1);
  gPad->SetMargin(0.10, 0.05, 0.12, 0.08);
  hRaw->Draw("colz");

  c.cd(2);
  gPad->SetMargin(0.10, 0.05, 0.12, 0.08);
  hSel->Draw("colz");

  const string tag = run_seg_tag_from_input(rootfile);
  
  const string png_path = join_path(
    cfg.out_dir,
    Form("%s_%s_diag_I_panels.png", cfg.out_prefix.c_str(), tag.c_str())
  );
  c.SaveAs(png_path.c_str());
  cout << "[diag] Saved " << png_path << "\n";
}

static void makeChargePeakPlot(const string& rootfile,
                               const Config& cfg,
                               const SelectionSettings& sel)
{
  RDataFrame dH("TSHelH", rootfile);

  auto d_pre = dH.Filter(Form("%s > %f",
                              sel.branch_scaler_current.c_str(),
                              sel.junk_floor_uA));

  const string hel_expr = build_helicity_filter_expr(sel);
  if (!hel_expr.empty()) d_pre = d_pre.Filter(hel_expr);

  if (sel.use_current_window) {
    double I0 = -1.0;
    if (sel.i0_mode == I0Mode::Fixed) {
      I0 = sel.i0_fixed_uA;
    } else {
      auto hI = d_pre.Histo1D(
        {"hI_tmp","BCM current;I [uA];counts", 120, 0.0, 60.0},
        sel.branch_scaler_current
      );
      I0 = hI->GetXaxis()->GetBinCenter(hI->GetMaximumBin());
    }

    const double Imin = (1.0 - sel.window_frac) * I0;
    const double Imax = (1.0 + sel.window_frac) * I0;

    d_pre = d_pre.Filter(
      Form("(%s >= %f) && (%s <= %f)",
           sel.branch_scaler_current.c_str(), Imin,
           sel.branch_scaler_current.c_str(), Imax)
    );
  }

  auto hQ = d_pre.Histo1D(
    {"hQ",
     "Helicity scaler charge/read after current-window cut; scaler charge/read; counts",
     300, 0.0, 1.0},
    sel.branch_scaler_charge
  );

  int ibinQ = hQ->GetMaximumBin();
  double Q_mode = hQ->GetXaxis()->GetBinCenter(ibinQ);

  double q_fit_halfwidth = 0.08;
  double q_fit_min = Q_mode - q_fit_halfwidth;
  double q_fit_max = Q_mode + q_fit_halfwidth;

  hQ->Fit("gaus", "RQ0", "", q_fit_min, q_fit_max);
  TF1 *fQ = hQ->GetFunction("gaus");

  if (!fQ) {
    cerr << "[charge-fit] ERROR: Gaussian fit failed.\n";
    return;
  }
  double Q_mu    = fQ->GetParameter(1);
  double Q_sigma = fabs(fQ->GetParameter(2));
  
  double Qmin = Q_mu * (1.0 - sel.max_charge_frac_spread);
  double Qmax = Q_mu * (1.0 + sel.max_charge_frac_spread);
  
  cout << "[charge-fit] mu            = " << Q_mu << " uC\n";
  cout << "[charge-fit] sigma         = " << Q_sigma << " uC\n";
  cout << "[charge-fit] allowed frac  = " << sel.max_charge_frac_spread << "\n";
  cout << "[charge-fit] allowed range = [" << Qmin << ", " << Qmax << "] uC\n";

  TCanvas cQ("cQ","Scaler charge peak scan",900,700);
  hQ->Draw();

  TLine lQmin(Qmin, 0.0, Qmin, hQ->GetMaximum()*1.05);
  TLine lQmax(Qmax, 0.0, Qmax, hQ->GetMaximum()*1.05);
  lQmin.SetLineColor(kRed+1);
  lQmax.SetLineColor(kRed+1);
  lQmin.SetLineStyle(2);
  lQmax.SetLineStyle(2);
  lQmin.SetLineWidth(2);
  lQmax.SetLineWidth(2);
  lQmin.Draw("same");
  lQmax.Draw("same");

  const string tag = run_seg_tag_from_input(rootfile);
  
  const string png_path = join_path(
    cfg.out_dir,
    Form("%s_%s_diag_Qscan.png", cfg.out_prefix.c_str(), tag.c_str())
  );

  cQ.SaveAs(png_path.c_str());
  cout << "[diag] Saved " << png_path << "\n";
}

static void makeRollingFracSpreadPlot(const string& rootfile,
                                      const Config& cfg,
                                      const SelectionSettings& sel)
{
  const PreselectedScalerData pre = build_preselected_scaler_data(rootfile, sel);
  const RollingSpreadData rs = compute_rolling_spread_data(pre, sel.stable_window_N);

  if (rs.evcount_eval.empty()) {
    cerr << "[diag] No rolling-spread points available for plotting.\n";
    return;
  }

  double xMin = rs.evcount_eval.front();
  double xMax = rs.evcount_eval.back();

  double yMaxData = 0.0;
  for (double y : rs.frac_spread) yMaxData = std::max(yMaxData, y);
  double yMax = std::max(1.2 * sel.max_charge_frac_spread,
                         1.05 * std::max(yMaxData, 1e-6));

  TGraph g((int)rs.evcount_eval.size(), rs.evcount_eval.data(), rs.frac_spread.data());
  g.SetTitle(Form("Rolling fractional charge spread (N=%d);evcount;fractional spread",
                  sel.stable_window_N));
  g.SetMarkerStyle(20);
  g.SetMarkerSize(0.6);

  TCanvas c("c_frac_spread","Rolling fractional spread",1200,700);
  g.Draw("AP");

  g.GetXaxis()->SetLimits(xMin, xMax);
  g.GetHistogram()->SetMinimum(0.0);
  g.GetHistogram()->SetMaximum(yMax);

  TLine lThr(xMin, sel.max_charge_frac_spread, xMax, sel.max_charge_frac_spread);
  lThr.SetLineColor(kRed+1);
  lThr.SetLineStyle(2);
  lThr.SetLineWidth(2);
  lThr.Draw("same");

  TLine lZero(xMin, 0.0, xMax, 0.0);
  lZero.SetLineColor(kBlue+1);
  lZero.SetLineStyle(3);
  lZero.SetLineWidth(2);
  lZero.Draw("same");

  const string tag = run_seg_tag_from_input(rootfile);
  const string png_path = join_path(
    cfg.out_dir,
    Form("%s_%s_diag_fracSpread.png", cfg.out_prefix.c_str(), tag.c_str())
  );

  c.SaveAs(png_path.c_str());
  cout << "[diag] Saved " << png_path << "\n";
}

static void makeRollingFracSpreadMaskPlot(const string& rootfile,
                                          const Config& cfg,
                                          const SelectionSettings& sel)
{
  const PreselectedScalerData pre = build_preselected_scaler_data(rootfile, sel);
  RollingSpreadData rs = compute_rolling_spread_data(pre, sel.stable_window_N);

  if (rs.evcount_eval.empty()) {
    cerr << "[diag] No rolling-spread mask points available for plotting.\n";
    return;
  }

  rs.pass_mask.reserve(rs.frac_spread.size());
  for (double y : rs.frac_spread) {
    rs.pass_mask.push_back((y <= sel.max_charge_frac_spread) ? 1.0 : 0.0);
  }

  double xMin = rs.evcount_eval.front();
  double xMax = rs.evcount_eval.back();

  TGraph g((int)rs.evcount_eval.size(), rs.evcount_eval.data(), rs.pass_mask.data());
  g.SetTitle(Form("Rolling-spread pass mask after junk/current/helicity preselection (N=%d);evcount;pass mask",
                  sel.stable_window_N));
  g.SetMarkerStyle(20);
  g.SetMarkerSize(0.7);

  TCanvas c("c_frac_spread_mask","Rolling spread pass mask",1200,500);
  g.Draw("AP");

  g.GetXaxis()->SetLimits(xMin, xMax);
  g.GetHistogram()->SetMinimum(-0.1);
  g.GetHistogram()->SetMaximum(1.1);

  TLine lPass(xMin, 1.0, xMax, 1.0);
  lPass.SetLineColor(kGreen+2);
  lPass.SetLineStyle(2);
  lPass.SetLineWidth(2);
  lPass.Draw("same");

  TLine lFail(xMin, 0.0, xMax, 0.0);
  lFail.SetLineColor(kRed+1);
  lFail.SetLineStyle(2);
  lFail.SetLineWidth(2);
  lFail.Draw("same");

  const string tag = run_seg_tag_from_input(rootfile);
  const string png_path = join_path(
    cfg.out_dir,
    Form("%s_%s_diag_fracSpreadMask.png", cfg.out_prefix.c_str(), tag.c_str())
  );

  c.SaveAs(png_path.c_str());
  cout << "[diag] Saved " << png_path << "\n";
}

// =================================================================================
// Report (wrapper only)
// =================================================================================
static void writeSelectionReport(const string& rootfile,
                                 const Config& cfg,
                                 const SelectionSettings& sel,
                                 const SelectionResult& pick)
{
  RDataFrame dH("TSHelH", rootfile);
  auto dH_sel = dH.Filter(pick.evcount_cut);

  auto sum_charge = [&](const string& cut) {
    return dH.Filter(cut).Sum<double>(sel.branch_scaler_charge).GetValue();
  };

  double q_raw_h0 = sum_charge(Form("%s == 0", sel.branch_helicity.c_str()));
  double q_raw_hm = sum_charge(Form("%s < 0",  sel.branch_helicity.c_str()));
  double q_raw_hp = sum_charge(Form("%s > 0",  sel.branch_helicity.c_str()));
  double q_raw_tot = q_raw_h0 + q_raw_hm + q_raw_hp;

  double q_cut_h0 = dH_sel.Filter(Form("%s == 0", sel.branch_helicity.c_str()))
                      .Sum<double>(sel.branch_scaler_charge).GetValue();
  double q_cut_hm = dH_sel.Filter(Form("%s < 0", sel.branch_helicity.c_str()))
                      .Sum<double>(sel.branch_scaler_charge).GetValue();
  double q_cut_hp = dH_sel.Filter(Form("%s > 0", sel.branch_helicity.c_str()))
                      .Sum<double>(sel.branch_scaler_charge).GetValue();
  double q_cut_tot = q_cut_h0 + q_cut_hm + q_cut_hp;

  auto frac = [](double num, double den) {
    return (den > 0.0) ? num / den : 0.0;
  };

  const string tag = run_seg_tag_from_input(rootfile);
  
  const string report_path = join_path(
   cfg.out_dir,
   Form("%s_%s_selection_report.txt", cfg.out_prefix.c_str(), tag.c_str())
  );
  
  ofstream report(report_path);
  report << "Input file: " << rootfile << "\n";
  report << "Selection summary\n";
  report << "-----------------\n";
  report << "Helicity mode: " << helicity_mode_to_string(pick.helicity_mode) << "\n";
  report << "Quartet snap applied: " << (pick.quartet_snap_applied ? "yes" : "no") << "\n";
  report << "Current window enabled: " << (sel.use_current_window ? "yes" : "no") << "\n";
  report << "I0 used [uA]: " << pick.i0_used_uA << "\n";
  report << "Imin used [uA]: " << pick.current_min_uA << "\n";
  report << "Imax used [uA]: " << pick.current_max_uA << "\n";
  report << "Mean current [uA]: " << pick.mean_current_uA << "\n";
  report << "Rolling stability window N: " << sel.stable_window_N << "\n";
  report << "Rolling fractional charge-range threshold: " << sel.max_charge_frac_spread << "\n";
  report << "Scaler reads before stability: " << pick.n_scaler_reads_pre << "\n";
  report << "Scaler reads after stability: " << pick.n_scaler_reads_post << "\n\n";

  report << "Accumulated helicity scaler charge BEFORE cuts [uC]\n";
  report << "  hel =  0 : " << q_raw_h0 << "\n";
  report << "  hel = -1 : " << q_raw_hm << "\n";
  report << "  hel = +1 : " << q_raw_hp << "\n";
  report << "  total    : " << q_raw_tot << "\n\n";

  report << "Accumulated helicity scaler charge AFTER cuts [uC]\n";
  report << "  hel =  0 : " << q_cut_h0 << "\n";
  report << "  hel = -1 : " << q_cut_hm << "\n";
  report << "  hel = +1 : " << q_cut_hp << "\n";
  report << "  total    : " << q_cut_tot << "\n\n";

  report << "Fraction of charge retained AFTER cuts\n";
  report << "  hel =  0 : " << frac(q_cut_h0, q_raw_h0) << "\n";
  report << "  hel = -1 : " << frac(q_cut_hm, q_raw_hm) << "\n";
  report << "  hel = +1 : " << frac(q_cut_hp, q_raw_hp) << "\n";
  report << "  total    : " << frac(q_cut_tot, q_raw_tot) << "\n\n";

  report << "Accepted scaler-read ranges (evcount)\n";
  for (const auto& r : pick.evcount_ranges) {
    report << "  " << r.lo << " -> " << r.hi << "\n";
  }
  report << "\n";

  report << "Accepted cumulative event-number ranges (evNumber)\n";
  for (const auto& r : pick.evNumber_ranges) {
    report << "  " << r.lo << " -> " << r.hi << "\n";
  }
  report << "\n";

  report << "Accepted T event ranges (g.evnum)\n";
  for (const auto& r : pick.gevnum_ranges) {
    report << "  " << r.lo << " -> " << r.hi << "\n";
  }
  report << "\n";

  report.close();
  cout << "[report] Wrote " << report_path << "\n";
}

// =================================================================================
// Main processing
// =================================================================================
int nps(const Config& cfg) {
  const string rootfile = build_input_path(cfg);
  cout << "Input ROOT: " << rootfile << "\n";

  gSystem->mkdir(cfg.out_dir.c_str(), /*recursive=*/true);

  SelectionSettings sel = make_selection_settings(cfg);
  SelectionResult pick = build_event_selection(rootfile, sel);

  if (!pick.ok) {
    cerr << "[selection] " << pick.message << "\n";
    return 2;
  }

  cout << "[selection] " << pick.message << "\n";
  cout << "[selection] helicity mode = " << helicity_mode_to_string(pick.helicity_mode) << "\n";
  cout << "[selection] mean current  = " << pick.mean_current_uA << " uA\n";
  cout << "[selection] I0 used       = " << pick.i0_used_uA << " uA\n";
  cout << "[selection] current range = [" << pick.current_min_uA
       << ", " << pick.current_max_uA << "] uA\n";
  cout << "[selection] scaler reads  = " << pick.n_scaler_reads_pre
       << " -> " << pick.n_scaler_reads_post << "\n";

  RDataFrame dT("T", rootfile);
  RDataFrame dH("TSHelH", rootfile);

  auto d_good_H = dH.Filter(pick.evcount_cut);

  auto d_good_T = dT.Filter(pick.gevnum_cut)
                     .Filter("(int)round(fEvtHdr.fEvtType) == 1");

  if (cfg.apply_physics_cuts) {
    d_good_T = d_good_T.Filter(cfg.physicsCuts);
  }

  if (cfg.make_plots) {
    makeCurrentPanels(rootfile, cfg, sel, pick);
    makeChargePeakPlot(rootfile, cfg, sel);
    makeRollingFracSpreadPlot(rootfile, cfg, sel);
    makeRollingFracSpreadMaskPlot(rootfile, cfg, sel);
  }
 
  if (cfg.write_report) {
    writeSelectionReport(rootfile, cfg, sel, pick);
  }

  if (cfg.write_snapshots) {
    const string tag = run_seg_tag_from_input(rootfile);
    
    TString outT = Form("%s_%s_T.root", cfg.out_prefix.c_str(), tag.c_str());
    TString outH = Form("%s_%s_TSHelH.root", cfg.out_prefix.c_str(), tag.c_str());

    string outT_path = join_path(cfg.out_dir, outT.Data());
    string outH_path = join_path(cfg.out_dir, outH.Data());

    vector<string> colsT = {
      "g.evnum",
      "fEvtHdr.fEvtType",
      "T.helicity.mps",
      "H.gtr.dp",
      "H.gtr.th",
      "H.gtr.ph",
      "H.react.z"
    };

    vector<string> colsH = {
      "evcount",
      "evNumber",
      "actualHelicity",
      "H.BCM4A_Hel.scalerCurrent",
      "H.BCM4A_Hel.scalerCharge"
    };

    d_good_T.Snapshot("T", outT_path.c_str(), colsT);
    d_good_H.Snapshot("TSHelH", outH_path.c_str(), colsH);

    cout << "Wrote:\n  " << outT_path << "\n  " << outH_path << "\n";
  }

  return 0;
}

// =================================================================================
// ROOT entry point
// =================================================================================
int NPS_good_events(const char* a0="--input", const char* a1="",
                    const char* a2=nullptr, const char* a3=nullptr,
                    const char* a4=nullptr, const char* a5=nullptr,
                    const char* a6=nullptr, const char* a7=nullptr,
                    const char* a8=nullptr, const char* a9=nullptr,
                    const char* a10=nullptr, const char* a11=nullptr,
                    const char* a12=nullptr, const char* a13=nullptr,
                    const char* a14=nullptr, const char* a15=nullptr,
                    const char* a16=nullptr, const char* a17=nullptr,
                    const char* a18=nullptr, const char* a19=nullptr) {

  vector<string> sargs;
  sargs.reserve(1 + 20);
  sargs.push_back("NPS_good_events");

  auto push_if = [&](const char* p) {
    if (p && string(p).size() > 0) sargs.emplace_back(p);
  };

  push_if(a0);  push_if(a1);
  push_if(a2);  push_if(a3);
  push_if(a4);  push_if(a5);
  push_if(a6);  push_if(a7);
  push_if(a8);  push_if(a9);
  push_if(a10); push_if(a11);
  push_if(a12); push_if(a13);
  push_if(a14); push_if(a15);
  push_if(a16); push_if(a17);
  push_if(a18); push_if(a19);

  vector<vector<char>> buffers;
  buffers.reserve(sargs.size());

  vector<char*> argv;
  argv.reserve(sargs.size());

  for (const auto& s : sargs) {
    buffers.emplace_back(s.begin(), s.end());
    buffers.back().push_back('\0');
    argv.push_back(buffers.back().data());
  }

  int argc = (int)argv.size();
  Config cfg = parse_args(argc, argv.data());

  if (cfg.input_root.empty()) {
    cerr << "Need --input\n";
    return 2;
  }

  return nps(cfg);
}
