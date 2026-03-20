#ifndef NPS_SELECTION_HELPER_H
#define NPS_SELECTION_HELPER_H

/*
===============================================================================
NPS_selection_helper.h
-------------------------------------------------------------------------------
Purpose:
  Self-contained reusable helper for beam-stability / helicity-aware event
  selection using Hall C NPS replay ROOT files containing trees "TSHelH" and "T".

What it returns:
  - accepted evcount ranges
  - accepted evNumber ranges
  - accepted g.evnum ranges
  - ready-to-use cut strings for TSHelH and T

Usage:
  #include "NPS_selection_helper.h"

  SelectionSettings sel;
  sel.helicity_mode = HelicityMode::QuartetPM;   // or QuartetAll / IgnoreHelicity
  sel.stable_window_N = 30;                      // N=1 effectively disables rolling stability
  sel.max_charge_frac_spread = 0.08;
  sel.window_frac = 0.15;

  auto pick = build_event_selection(run, seg, rootfile, sel);
  if (!pick.ok) {
    std::cerr << pick.message << std::endl;
    return;
  }
  
  ROOT::RDataFrame dT("T", rootfile);
  auto d_good_T = dT.Filter(pick.gevnum_cut);
  
  ROOT::RDataFrame dH("TSHelH", rootfile);
  auto d_good_H = dH.Filter(pick.evcount_cut);
    
Notes:
  - QuartetPM    : keep only helicity != 0 and quartet-snap accepted ranges
  - QuartetAll   : keep all helicity states and quartet-snap accepted ranges
  - IgnoreHelicity: ignore helicity entirely and do NOT quartet-snap
  - Rolling charge stability is always part of the helper.
    Set stable_window_N = 1 to make it effectively off.
  - This helper preserves the original evcount -> event-number indexing convention.
===============================================================================
*/

#include <ROOT/RDataFrame.hxx>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using std::string;
using std::vector;

// =================================================================================
// Reusable event-selection helper types
// =================================================================================
enum class HelicityMode {
  QuartetPM,      // keep only helicity ±1, then quartet-snap accepted regions
  QuartetAll,     // keep helicity -1,0,+1, then quartet-snap accepted regions
  IgnoreHelicity  // do not cut on helicity, do not quartet-snap
};

enum class I0Mode {
  Peak,   // use mode of scalerCurrent histogram above junk floor
  Fixed   // use user-provided fixed I0
};

struct SelectionSettings {
  // --- tunables ---
  HelicityMode helicity_mode = HelicityMode::QuartetPM;

  double junk_floor_uA = 2.5;         // reject junk current reads under this threshold
  double mean_current_calc_min = 2.5; // floor for computing mean current
  double mean_current_min = 2.75;     // production runs were as low as 3 uA

  bool use_current_window = true;
  I0Mode i0_mode = I0Mode::Peak;
  double i0_fixed_uA = -1.0;          // used only when i0_mode == Fixed
  double window_frac = 0.15;          // keep |I-I0|/I0 <= window_frac

  int stable_window_N = 30;              // N=1 effectively disables the rolling requirement
  double max_charge_frac_spread = 0.08;  // (qmax-qmin)/mean within rolling window

  // branch names
  string branch_scaler_current = "H.BCM4A_Hel.scalerCurrent";
  string branch_scaler_charge  = "H.BCM4A_Hel.scalerCharge";
  string branch_helicity       = "actualHelicity";
  string branch_evcount        = "evcount";
  string branch_evNumber       = "evNumber";
  string branch_evnum_T        = "g.evnum";
};

struct RangeI {
  int lo = 0;
  int hi = 0;
};

struct SelectionResult {
  bool ok = false;
  string rootfile;
  string message;

  // bookkeeping / diagnostics
  HelicityMode helicity_mode = HelicityMode::QuartetPM;
  bool quartet_snap_applied = true;

  double mean_current_uA = 0.0;
  double i0_used_uA = -1.0;
  double current_min_uA = -1.0;
  double current_max_uA = -1.0;

  int n_scaler_reads_pre = 0;
  int n_scaler_reads_post = 0;

  // accepted regions
  vector<RangeI> evcount_ranges;   // accepted evcount ranges
  vector<RangeI> evNumber_ranges;  // mapped from TSHelH evNumber
  vector<RangeI> gevnum_ranges;    // mapped from T g.evnum

  // ready-to-use filter strings
  string evcount_cut;
  string evNumber_cut;
  string gevnum_cut;
};

// =================================================================================
// Small helper utilities
// =================================================================================
static inline string helicity_mode_to_string(HelicityMode mode) {
  switch (mode) {
    case HelicityMode::QuartetPM:      return "QuartetPM";
    case HelicityMode::QuartetAll:     return "QuartetAll";
    case HelicityMode::IgnoreHelicity: return "IgnoreHelicity";
  }
  return "Unknown";
}

static inline bool helicity_requires_quartet_snap(HelicityMode mode) {
  return (mode == HelicityMode::QuartetPM || mode == HelicityMode::QuartetAll);
}

static inline string build_helicity_filter_expr(const SelectionSettings& sel) {
  const string& h = sel.branch_helicity;

  switch (sel.helicity_mode) {
    case HelicityMode::QuartetPM:
      return Form("(int)round(%s) != 0", h.c_str());

    case HelicityMode::QuartetAll:
      return ""; // keep all helicity states

    case HelicityMode::IgnoreHelicity:
      return ""; // completely ignore helicity
  }
  return "";
}

static inline vector<RangeI> toRangeI(const vector<vector<int>>& vv) {
  vector<RangeI> out;
  out.reserve(vv.size());
  for (const auto& v : vv) {
    if (v.size() >= 2) out.push_back({v[0], v[1]});
  }
  return out;
}

static inline string build_evcount_cut_string(const vector<RangeI>& ranges,
                                              const string& branch_evcount) {
  string out;
  for (const auto& r : ranges) {
    string term = Form("(%d<(int)round(%s)&&(int)round(%s)<=%d)",
                       r.lo, branch_evcount.c_str(), branch_evcount.c_str(), r.hi);
    out = out.empty() ? term : out + "||" + term;
  }
  return out;
}

static inline string build_event_cut_string(const vector<RangeI>& ranges,
                                            const string& branch_event) {
  string out;
  for (const auto& r : ranges) {
    string term = Form("(%d<=(int)round(%s)&&(int)round(%s)<%d)",
                       r.lo, branch_event.c_str(), branch_event.c_str(), r.hi);
    out = out.empty() ? term : out + "||" + term;
  }
  return out;
}

static inline vector<int> take_int_sorted(ROOT::RDataFrame& df, const string& branch) {
  auto valsD = df.Take<double>(branch).GetValue();
  vector<int> vals(valsD.begin(), valsD.end());
  std::sort(vals.begin(), vals.end());
  return vals;
}

// =================================================================================
// Quartet alignment helpers
// =================================================================================
template<typename RDF>
static inline int alignStartToQuartetEnd(int evnum, const RDF& /*dfHelper*/) {
  return evnum - (evnum % 4);
}

template<typename RDF>
static inline int alignEndToQuartetEnd(int evnum, const RDF& /*dfHelper*/) {
  return evnum + (3 - (evnum % 4));
}

template<typename RDF>
static inline vector<vector<int>> buildQuartetAlignedRanges(const vector<int>& evC_arr,
                                                            const RDF& dfHelper) {
  vector<vector<int>> evtCRanges;
  if (evC_arr.empty()) return evtCRanges;

  size_t run_start_idx = 0;

  for (size_t i = 1; i <= evC_arr.size(); ++i) {
    bool end_of_run = (i == evC_arr.size()) || (evC_arr[i] != evC_arr[i - 1] + 1);
    if (!end_of_run) continue;

    int run_start_ev = evC_arr[run_start_idx];
    int run_end_ev   = evC_arr[i - 1];
    int run_len      = run_end_ev - run_start_ev + 1;

    if (run_len > 0) {
      vector<int> range;
      range.push_back(alignStartToQuartetEnd(run_start_ev, dfHelper));
      range.push_back(alignEndToQuartetEnd(run_end_ev, dfHelper));
      evtCRanges.push_back(range);
    }

    run_start_idx = i;
  }

  return evtCRanges;
}

static inline vector<vector<int>> mergeOverlappingRanges(const vector<vector<int>>& ranges) {
  vector<vector<int>> merged;
  if (ranges.empty()) return merged;

  merged.push_back(ranges[0]);

  for (size_t i = 1; i < ranges.size(); ++i) {
    int start = ranges[i][0];
    int end   = ranges[i][1];

    if (start <= merged.back()[1] + 1) {
      merged.back()[1] = std::max(merged.back()[1], end);
    } else {
      merged.push_back(ranges[i]);
    }
  }

  return merged;
}

// =================================================================================
// Rolling charge-stability helper
// =================================================================================
static inline vector<int> stableEvcountsFromChargeWindow(const vector<int>& evC_arr,
                                                         const vector<double>& Q_arr,
                                                         int N,
                                                         double fracRangeMax)
{
  vector<int> out;

  if ((int)evC_arr.size() != (int)Q_arr.size()) return out;
  if (N < 1) return out;                  // invalid
  if ((int)evC_arr.size() < N) return out;

  int runStart = 0;
  for (int i = 1; i <= (int)evC_arr.size(); i++) {
    bool endRun = (i == (int)evC_arr.size()) || (evC_arr[i] != evC_arr[i - 1] + 1);
    if (!endRun) continue;

    int runEnd = i - 1;
    int runLen = runEnd - runStart + 1;
    if (runLen >= N) {
      double sum = 0.0;
      std::deque<int> dqMin, dqMax;

      for (int j = runStart; j <= runEnd; j++) {
        sum += Q_arr[j];

        while (!dqMin.empty() && Q_arr[dqMin.back()] >= Q_arr[j]) dqMin.pop_back();
        dqMin.push_back(j);

        while (!dqMax.empty() && Q_arr[dqMax.back()] <= Q_arr[j]) dqMax.pop_back();
        dqMax.push_back(j);

        int old = j - N;
        if (old >= runStart) {
          sum -= Q_arr[old];
          if (!dqMin.empty() && dqMin.front() == old) dqMin.pop_front();
          if (!dqMax.empty() && dqMax.front() == old) dqMax.pop_front();
        }

        if (j >= runStart + (N - 1)) {
          double mean = sum / (double)N;
          if (mean > 0) {
            double qmin = Q_arr[dqMin.front()];
            double qmax = Q_arr[dqMax.front()];
            double fracRange = (qmax - qmin) / mean;

            if (fracRange <= fracRangeMax) {
              out.push_back(evC_arr[j]);
            }
          }
        }
      }
    }

    runStart = i;
  }

  return out;
}

// =================================================================================
// Reusable selection engine
// =================================================================================
static inline SelectionResult build_event_selection(const string& rootfile, const SelectionSettings& sel)
{
  SelectionResult res;
  res.rootfile = rootfile;
  res.helicity_mode = sel.helicity_mode;
  res.quartet_snap_applied = helicity_requires_quartet_snap(sel.helicity_mode);

  ROOT::RDataFrame d_shelp("TSHelH", rootfile);
  ROOT::RDataFrame d_T("T", rootfile);

  // ---------------------------------------------------------
  // Stable window sanity check
  // ---------------------------------------------------------
  
  if (sel.stable_window_N < 1) {
    res.ok = false;
    res.message = "Invalid stable_window_N: must be >= 1.";
    return res;
  }

  // ---------------------------------------------------------
  // Mean current sanity check
  // ---------------------------------------------------------
  res.mean_current_uA = d_shelp
    .Filter(Form("%s > %f", sel.branch_scaler_current.c_str(), sel.mean_current_calc_min))
    .Mean(sel.branch_scaler_current)
    .GetValue();

  if (res.mean_current_uA < sel.mean_current_min) {
    res.ok = false;
    res.message = Form("Mean current too low (%g uA < %g uA)",
                       res.mean_current_uA, sel.mean_current_min);
    return res;
  }

  // ---------------------------------------------------------
  // Determine I0 if current window is requested
  // ---------------------------------------------------------
  if (sel.use_current_window) {
    if (sel.i0_mode == I0Mode::Fixed) {
      res.i0_used_uA = sel.i0_fixed_uA;
    } else {
      auto d_Ivalid = d_shelp.Filter(
        Form("%s > %f", sel.branch_scaler_current.c_str(), sel.junk_floor_uA)
      );

      auto hI = d_Ivalid.Histo1D(
        {"hI_helper","BCM current;I [uA];counts", 120, 0.0, 60.0},
        sel.branch_scaler_current
      );
      res.i0_used_uA = hI->GetXaxis()->GetBinCenter(hI->GetMaximumBin());
    }

    if (!(res.i0_used_uA > 0.0)) {
      res.ok = false;
      res.message = Form("Invalid I0 value (%g uA)", res.i0_used_uA);
      return res;
    }

    res.current_min_uA = (1.0 - sel.window_frac) * res.i0_used_uA;
    res.current_max_uA = (1.0 + sel.window_frac) * res.i0_used_uA;
  }

  // ---------------------------------------------------------
  // Preselection on TSHelH
  // ---------------------------------------------------------
  auto d_sel = d_shelp.Filter(
    Form("%s > %f", sel.branch_scaler_current.c_str(), sel.junk_floor_uA)
  );

  const string hel_expr = build_helicity_filter_expr(sel);
  if (!hel_expr.empty()) d_sel = d_sel.Filter(hel_expr);

  if (sel.use_current_window) {
    d_sel = d_sel.Filter(
      Form("(%s >= %f) && (%s <= %f)",
           sel.branch_scaler_current.c_str(), res.current_min_uA,
           sel.branch_scaler_current.c_str(), res.current_max_uA)
    );
  }

  {
    auto evC_pre_D = d_sel.Take<double>(sel.branch_evcount).GetValue();
    res.n_scaler_reads_pre = (int)evC_pre_D.size();
  }

  // ---------------------------------------------------------
  // Extract evcount + charge arrays and sort by evcount
  // ---------------------------------------------------------
  auto evC_arr_D = d_sel.Take<double>(sel.branch_evcount).GetValue();
  auto Q_arr_D   = d_sel.Take<double>(sel.branch_scaler_charge).GetValue();

  vector<int> evC_arr(evC_arr_D.begin(), evC_arr_D.end());
  vector<double> Q_arr(Q_arr_D.begin(), Q_arr_D.end());

  vector<int> idx(evC_arr.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b){ return evC_arr[a] < evC_arr[b]; });

  vector<int> evC_sorted;
  vector<double> Q_sorted;
  evC_sorted.reserve(evC_arr.size());
  Q_sorted.reserve(Q_arr.size());

  for (int k : idx) {
    evC_sorted.push_back(evC_arr[k]);
    Q_sorted.push_back(Q_arr[k]);
  }

  evC_arr.swap(evC_sorted);
  Q_arr.swap(Q_sorted);

  // ---------------------------------------------------------
  // Rolling charge stability (always part of helper)
  // N=1 effectively disables the rolling requirement
  // ---------------------------------------------------------
  evC_arr = stableEvcountsFromChargeWindow(
    evC_arr,
    Q_arr,
    sel.stable_window_N,
    sel.max_charge_frac_spread
  );

  res.n_scaler_reads_post = (int)evC_arr.size();

  if (evC_arr.empty()) {
    res.ok = false;
    res.message = "No scaler reads survived selection.";
    return res;
  }

  // ---------------------------------------------------------
  // Build accepted evcount ranges
  // ---------------------------------------------------------
  vector<vector<int>> evtCRanges_raw;

  if (helicity_requires_quartet_snap(sel.helicity_mode)) {
    evtCRanges_raw = buildQuartetAlignedRanges(evC_arr, d_sel);
  } else {
    size_t run_start_idx = 0;
    for (size_t i = 1; i <= evC_arr.size(); ++i) {
      bool end_of_run = (i == evC_arr.size()) || (evC_arr[i] != evC_arr[i - 1] + 1);
      if (!end_of_run) continue;

      evtCRanges_raw.push_back({evC_arr[run_start_idx], evC_arr[i - 1]});
      run_start_idx = i;
    }
  }

  auto evtCRanges_merged_vv = mergeOverlappingRanges(evtCRanges_raw);
  res.evcount_ranges = toRangeI(evtCRanges_merged_vv);

  // ---------------------------------------------------------
  // Build mapping arrays
  // ---------------------------------------------------------
  vector<int> evN_arr   = take_int_sorted(d_shelp, sel.branch_evNumber);
  vector<int> gEvnumArr = take_int_sorted(d_T, sel.branch_evnum_T);

  // ---------------------------------------------------------
  // Map evcount ranges -> evNumber and g.evnum ranges
  // NOTE: preserves original indexing convention
  // ---------------------------------------------------------
  for (const auto& r : res.evcount_ranges) {
    if (r.lo < 1 || r.hi < 1) continue;
    if (r.lo > (int)evN_arr.size() || r.hi > (int)evN_arr.size()) continue;
    if (r.lo > (int)gEvnumArr.size() || r.hi > (int)gEvnumArr.size()) continue;

    res.evNumber_ranges.push_back({evN_arr[r.lo - 1], evN_arr[r.hi - 1]});
    res.gevnum_ranges.push_back({gEvnumArr[r.lo - 1], gEvnumArr[r.hi - 1]});
  }

  // ---------------------------------------------------------
  // Build cut strings
  // ---------------------------------------------------------
  res.evcount_cut = build_evcount_cut_string(res.evcount_ranges, sel.branch_evcount);
  res.evNumber_cut = build_event_cut_string(res.evNumber_ranges, sel.branch_evNumber);
  res.gevnum_cut   = build_event_cut_string(res.gevnum_ranges, sel.branch_evnum_T);

  res.ok = true;
  res.message = "Selection built successfully.";
  return res;
}

#endif
