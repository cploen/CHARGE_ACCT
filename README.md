# Charge Accounting / Event Selection Helper

Small ROOT/C++ utilities for selecting "good" beam periods directly from the scaler behavior observed in a run segment.

The goal is to build a usable event selector from what the run segment actually did, rather than depending only on an external nominal-current value. The selection is configurable so users can make it stricter or looser as needed.

## Files

- `NPS_selection_helper.h`  
  Core reusable selection logic.

- `NPS_good_events.C`  
  Standalone wrapper around the helper. Adds argument parsing, optional plots, a text report, and reduced ROOT snapshots.

- `MinWorkingEx_NPS_selection_helper.C`  
  Minimal example showing how to call the helper directly.

## What the selector can use

At the scaler level (`TSHelH`), the helper can apply:

- a junk-current floor
- helicity handling (`quartet_pm`, `quartet_all`, or `ignore`)
- a current-window cut around either:
  - the observed current peak, or
  - a fixed user-supplied `I0`
- a rolling local charge/read stability cut

It then maps the accepted scaler-read regions onto:
- `evcount`
- `evNumber`
- `g.evnum`

so the accepted region can be applied to both `TSHelH` and `T`.

## Default wrapper behavior

`NPS_good_events.C` defaults to:

- helicity mode: `quartet_pm`
- current window: enabled
- `I0`: taken from the scaler-current peak unless `--I0` is given
- current window fraction: `0.15`
- rolling stability window: `N = 30`
- max rolling fractional charge spread: `0.08`
- physics cuts on `T`: off by default
- report writing: on
- reduced snapshots: on
- plots: off

## Quick use

Run the wrapper on a replay file:

```bash
root -l -b -q 'NPS_good_events.C("--input","/path/to/file.root")'

Make the diagnostic plots too:

root -l -b -q 'NPS_good_events.C("--input","/path/to/file.root","--plots")'

## Additional Options and Notes

Common options
	•	--input <path>
Input replay ROOT file.
	•	--helicity <quartet_pm|quartet_all|ignore>
Choose helicity handling mode.
	•	--no-current-window
Disable current-window selection.
	•	--I0 <uA>
Use a fixed nominal current instead of the observed peak.
	•	--window-frac <f>
Set current window size.
	•	--stable-window-N <int>
Set rolling stability window length.
	•	--physics-cuts
Apply optional HMS physics cuts to T.
	•	--plots
Write diagnostic PNGs.
	•	--no-report
Suppress the text report.
	•	--no-snapshots
Suppress reduced ROOT outputs.
	•	--outdir <dir>
Output directory.
	•	--outprefix <str>
Output filename prefix.

Minimal helper use

If you want to use the selector inside your own code, see:
	•	MinWorkingEx_NPS_selection_helper.C

The basic pattern is:
	1.	fill SelectionSettings
	2.	call build_event_selection(...)
	3.	apply the returned cut strings to TSHelH and/or T

Notes on the plots

The rolling-spread diagnostic plots are built from the scaler stream after:
	•	junk-floor cut
	•	helicity filter
	•	current-window cut

A point on the rolling fractional spread plot at a given evcount corresponds to the full trailing N-read window ending there. It is not a single-read quantity.

Intended use

This tool is useful for:
	•	precision helicity accounting
	•	charge-normalized yields
	•	checking whether a run segment contains a stable usable beam plateau
	•	producing a reproducible event selector based on the observed scaler behavior
