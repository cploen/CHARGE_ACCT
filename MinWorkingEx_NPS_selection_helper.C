// File: demo_use_selection_helper.C
#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include "NPS_selection_helper.h"

void MinWorkingEx_NPS_selection_helper(
    const char* input_root =
      "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_6786_7_1_-1.root")
{
  std::string rootfile = input_root;

  // 1. Configure the helper
  SelectionSettings sel;
  sel.helicity_mode = HelicityMode::QuartetPM;  // or QuartetAll / IgnoreHelicity
  sel.window_frac = 0.15;
  sel.stable_window_N = 30;
  sel.max_charge_frac_spread = 0.08;

  // Optional:
  // sel.use_current_window = false;
  // sel.I0_fixed_uA = 3.0;   // only if also setting sel.i0_mode = I0Mode::Fixed
  // sel.i0_mode = I0Mode::Fixed;

  // 2. Build the selection
  SelectionResult pick = build_event_selection(rootfile, sel);

  if (!pick.ok) {
    std::cerr << "[selection] " << pick.message << std::endl;
    return;
  }

  // 3. Apply returned cuts to the trees you care about
  ROOT::RDataFrame dT("T", rootfile);
  ROOT::RDataFrame dH("TSHelH", rootfile);

  auto d_good_T = dT.Filter(pick.gevnum_cut);
  auto d_good_H = dH.Filter(pick.evcount_cut);

  // 4. Example: count surviving entries
  auto nT = d_good_T.Count();
  auto nH = d_good_H.Count();

  std::cout << "Selection built successfully\n";
  std::cout << "helicity mode     = " << helicity_mode_to_string(pick.helicity_mode) << "\n";
  std::cout << "I0 used [uA]      = " << pick.i0_used_uA << "\n";
  std::cout << "good T entries    = " << nT.GetValue() << "\n";
  std::cout << "good TSHelH reads = " << nH.GetValue() << "\n";
}
