"""Transform and plotting functions for pulsed / ferroelectric measurements.

Companion to ``PMU.py``, in the same relationship as ``SMU_functions.py`` is to
``SMU.py``.  The module has two halves:

**Transforms** - ``extract_pund``, ``extract_dhm``, ``extract_fatigue``.  These
take the :class:`PMUdata` / :class:`FatigueRun` objects produced by the
importers, compute the processed loop and the derived scalars, and write both
back into the object (``.pund`` / ``.dhm`` / ``.extracted_params``).  Extraction
is deliberately separate from plotting so the same numbers can be plotted,
tabulated and exported without being recomputed each time.

**Plots** - ``plot_PUND`` (waveforms vs time), ``plot_PUND_polarisation``
(switching loop), ``plot_DHM`` (hysteresis loop), ``plot_fatigue``
(polarisation vs cycle count) and ``plot_custom_waveform``.  Each will run the
matching transform on demand, so a plot call on a freshly imported object just
works.
"""

# import all the libraries needed
from import_dep import *
from scipy.signal import medfilt
from scipy.ndimage import uniform_filter1d
from scipy.integrate import cumulative_trapezoid

from PMU import PMUdata, FatigueRun, PMUContainer


# =============================================================================
#  Input handling
# =============================================================================

def _as_run_list(pmu_data) -> list:
    """Normalise the many things a caller might pass into a list of runs.

    Accepts a container object, a list of containers, a list of
    :class:`PMUdata` / :class:`FatigueRun`, or a single one of those.  This is
    what lets ``plot_DHM(fatigue_run.sub_runs)`` work without any special case.
    """
    if pmu_data is None:
        return []
    if isinstance(pmu_data, (PMUdata, FatigueRun)):
        return [pmu_data]
    if isinstance(pmu_data, PMUContainer):
        return list(pmu_data)

    out = []
    for item in pmu_data:
        out.extend(_as_run_list(item))
    return out


def _get_pmu_run_nums(run_nums, index):
    """Per-object run selection, mirroring ``_get_smu_run_nums``.

    ``run_nums`` may be None (all runs), a flat list applied to every object, or
    a list of lists - one selection per object in ``pmu_data``.
    """
    if run_nums is None:
        return None
    if len(run_nums) > 0 and isinstance(run_nums[0], (list, tuple)):
        return list(run_nums[index]) if index < len(run_nums) else None
    return list(run_nums)


def _select_runs(pmu_data, run_nums, meas_type: Optional[str] = None) -> list:
    """Resolve ``pmu_data`` + ``run_nums`` into the runs to plot."""
    items = pmu_data if isinstance(pmu_data, (list, tuple)) else [pmu_data]

    runs = []
    for idx, item in enumerate(items):
        if isinstance(item, PMUContainer):
            selected, _ = item.get(_get_pmu_run_nums(run_nums, idx))
            runs.extend(selected)
        else:
            runs.extend(_as_run_list(item))

    if meas_type is not None:
        wrong = [r for r in runs
                 if getattr(r, 'meas_type', meas_type) != meas_type]
        if wrong:
            print(f"Note: skipping {len(wrong)} run(s) that are not {meas_type}.")
        runs = [r for r in runs if getattr(r, 'meas_type', meas_type) == meas_type]
    return runs


def _first_thickness(runs) -> Optional[float]:
    """Thickness for the secondary field axis; warn if the runs disagree."""
    values = {r.film_thickness_nm for r in runs if r.film_thickness_nm is not None}
    if not values:
        return None
    if len(values) > 1:
        print(f"Warning: runs have different film thicknesses {sorted(values)}; "
              f"the electric-field axis uses {min(values)} nm.")
    return min(values)


def _smooth(y, med_filt, sigma_filt):
    """Optional median then boxcar filtering of a current trace."""
    y = np.asarray(y, dtype=float)
    if med_filt is not None and med_filt > 0:
        kernel = med_filt if med_filt % 2 != 0 else med_filt + 1
        if kernel < len(y):
            y = medfilt(y, kernel_size=kernel)
    if sigma_filt is not None and sigma_filt > 0:
        y = uniform_filter1d(y, size=sigma_filt)
    return y


#: Colours for a twin-axis figure showing a *single* dataset: black on the left,
#: red on the right. With one dataset there is no colour information to carry, so
#: colour is better spent saying which trace belongs to which axis. With several
#: datasets colour has to distinguish them instead, and the project cycle is used.
SINGLE_LEFT = 'k'
SINGLE_RIGHT = '#C00000'

#: Axis symbols for the custom-waveform channels: ch2 drives the gate, ch1 the drain
_CH_SYM = {'V_Ch1': 'V_D', 'V_Ch2': 'V_G', 'I_Ch1': 'I_D', 'I_Ch2': 'I_G'}


def _twin_colours(n_datasets, prop_colors, index):
    """(left, right) colours for dataset ``index`` on a twin-axis figure."""
    if n_datasets == 1:
        return SINGLE_LEFT, SINGLE_RIGHT
    colour = prop_colors[index % len(prop_colors)]
    return colour, colour


def _colour_right_axis(ax_right, colour):
    """Tint a right-hand axis label, ticks and spine to match its trace."""
    ax_right.yaxis.label.set_color(colour)
    ax_right.tick_params(axis='y', colors=colour)
    ax_right.spines['right'].set_color(colour)


def _set_title(ax, plot_title):
    """Apply an explicit title, or clear it.

    Titles are off by default across the project: a figure caption carries the
    description in a paper, and a title duplicated above the axes is the first
    thing a journal asks to be removed.
    """
    ax.set_title(plot_title if plot_title else "")


def _plot_kwargs(display_points, line_width, markersize) -> dict:
    kwargs = {'marker': 'o' if display_points else ''}
    if line_width is not None:
        kwargs['linewidth'] = line_width
    if markersize is not None:
        kwargs['markersize'] = markersize
    return kwargs


def _finish(fig, export_data, output_PMU, fig_name, suffix, fig_format,
            plot_transparency, show, extra_csv=None):
    """Save the figure (and optionally a summary CSV), then show or close it."""
    if export_data and output_PMU is not None:
        output_dir = Path(output_PMU)
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = f"{fig_name}_{suffix}" if fig_name else suffix
        fig.savefig(str(output_dir / f"{stem}.{fig_format}"), dpi=600,
                    bbox_inches='tight', transparent=plot_transparency,
                    format=fig_format)
        if extra_csv is not None and not extra_csv.empty:
            extra_csv.to_csv(str(output_dir / f"{stem}_summary.csv"), index=False)

    if show:
        plt.show()
    else:
        plt.close(fig)


def _align_right_zero(ax_left, ax_right):
    """Line up zero on the right-hand axis with zero on the left.

    Without this the current trace floats at an arbitrary height and its peaks
    cannot be read against the polarisation steps.
    """
    l0, l1 = ax_left.get_ylim()
    if not (l0 < 0.0 < l1):
        return
    frac = (0.0 - l0) / (l1 - l0)
    if not (0.02 < frac < 0.98):
        return
    r0, r1 = ax_right.get_ylim()
    needed = [0.0]
    if r0 < 0:
        needed.append(-r0 / frac)
    if r1 > 0:
        needed.append(r1 / (1.0 - frac))
    span = max(needed) * 1.05
    if span > 0:
        ax_right.set_ylim(-frac * span, (1.0 - frac) * span)


def _field_axis(ax, thickness_nm):
    """Add a secondary electric-field axis (kV/cm) along the top."""
    if thickness_nm is None:
        return None
    d_cm = float(thickness_nm) * 1e-7

    def to_field(V):
        return np.asarray(V, dtype=float) / d_cm / 1e3

    def from_field(E):
        return np.asarray(E, dtype=float) * d_cm * 1e3

    sec = ax.secondary_xaxis('top', functions=(to_field, from_field))
    sec.set_xlabel(r"Electric field (kV/cm)")
    return sec


def _to_field(V, thickness_nm):
    """Volts -> kV/cm for a given film thickness in nm."""
    if thickness_nm is None:
        return np.nan
    return float(V) / (float(thickness_nm) * 1e-7) / 1e3


# =============================================================================
#  PUND pulse handling
# =============================================================================

def _split_pund_pulses(V, I, t, rel_threshold: float = 0.05,
                       min_frac: float = 0.02):
    """Split a PUND time-series into individual pulses.

    Pulses are detected as contiguous regions where |V| exceeds
    *rel_threshold* x max|V| (a *relative* threshold, so the same value works
    for a 2 V and a 20 V waveform).  Each detected region is then widened
    outwards until |V| falls to the measurement noise floor, so that every
    pulse starts and ends at V ~ 0.  Without this the branches would be
    clipped part-way up the ramp and the reconstructed loop could not close.

    Regions shorter than *min_frac* of the record are discarded as glitches.

    Returns a list of dicts, each containing 'V', 'I', 't' arrays for one
    pulse plus 'sign' (+1 or -1) for polarity and 'peak_idx', the index of
    the voltage extremum (which separates the rising ramp from the
    hold/falling part of a trapezoidal pulse).
    """
    Vmax = np.max(np.abs(V))
    if Vmax == 0:
        return []

    active = np.abs(V) > rel_threshold * Vmax

    # Noise floor of the quiescent (between-pulse) baseline
    quiet = np.abs(V[~active])
    v_noise = 5.0 * np.median(quiet) if quiet.size else rel_threshold * Vmax

    edges = np.diff(active.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0]
    if active[0]:
        starts = np.r_[0, starts]
    if active[-1]:
        ends = np.r_[ends, len(V) - 1]

    min_len = int(min_frac * len(V))
    keep = [(s, e) for s, e in zip(starts, ends) if (e - s) >= min_len]

    pulses = []
    for k, (s, e) in enumerate(keep):
        # Widen down to the noise floor, but never past the midpoint of the
        # gap to the neighbouring pulse.
        lo = 0 if k == 0 else (keep[k - 1][1] + s) // 2
        hi = len(V) - 1 if k == len(keep) - 1 else (e + keep[k + 1][0]) // 2
        while s > lo and abs(V[s - 1]) > v_noise:
            s -= 1
        while e < hi and abs(V[e + 1]) > v_noise:
            e += 1

        v_seg = V[s:e + 1]
        pk = int(np.argmax(np.abs(v_seg)))
        pulses.append({
            'V': v_seg,
            'I': I[s:e + 1],
            't': t[s:e + 1],
            'sign': +1 if v_seg[pk] > 0 else -1,
            'peak_idx': pk,
        })
    return pulses


def _pulses_from_frame(df, med_filt=None, sigma_filt=None) -> list:
    """Build the pulse list for one measurement.

    When the importer already knows where the pulses are - the aixACCT export
    writes each pulse as its own column group, recorded in a ``pulse`` column -
    that grouping is used directly, which is exact.  Otherwise the pulses are
    recovered from the voltage trace by :func:`_split_pund_pulses`.
    """
    V = df['V'].to_numpy(dtype=float)
    I = _smooth(df['I'].to_numpy(dtype=float), med_filt, sigma_filt)
    t = df['Time'].to_numpy(dtype=float)

    if 'pulse' not in df.columns:
        return _split_pund_pulses(V, I, t)

    pulses = []
    for tag in pd.unique(df['pulse']):
        m = (df['pulse'] == tag).to_numpy()
        v_seg = V[m]
        if v_seg.size == 0:
            continue
        pk = int(np.argmax(np.abs(v_seg)))
        pulses.append({'V': v_seg, 'I': I[m], 't': t[m],
                       'sign': +1 if v_seg[pk] > 0 else -1, 'peak_idx': pk})
    return pulses


def _assign_pund_pulses(pulses, label_str: str = '', p_pulse: str = 'last'):
    """Map a list of detected pulses onto the P, U, N, D roles.

    Handles the waveforms seen in practice:

    * 4 pulses, ``+ + - -``      -> P U N D  (Keithley 4200 ``nvm`` pundTest)
    * 5 pulses, ``- + + - -``    -> leading preset discarded, then P U N D
    * 5 pulses, ``+ + - - +``    -> P U N D P  (aixACCT TFAnalyzer)

    In the aixACCT case there are two positive switching pulses.  The *last* one
    is preferred by default: it follows D, so the film is definitively
    down-poled before it and the charge it switches is unambiguous.  The first P
    has no such guarantee - what it measures depends on how the device was left
    by whatever ran before - so ``p_pulse='first'`` should be used knowingly.

    Returns ``(P, U, N, D)``, or None if the sequence cannot be interpreted.
    """
    signs = [p['sign'] for p in pulses]

    # aixACCT: P U N D P
    if len(pulses) == 5 and signs == [1, 1, -1, -1, 1]:
        P = pulses[0] if p_pulse == 'first' else pulses[4]
        return P, pulses[1], pulses[2], pulses[3]

    if len(pulses) == 5 and signs[0] < 0:
        print(f"  '{label_str}': 5 pulses found, leading negative pulse "
              f"treated as a preset and discarded.")
        pulses = pulses[1:]
    elif len(pulses) > 4:
        print(f"  '{label_str}': {len(pulses)} pulses found, using the last 4.")
        pulses = pulses[-4:]

    if len(pulses) != 4:
        print(f"Warning: '{label_str}' - expected 4 PUND pulses, found "
              f"{len(pulses)} - skipping.")
        return None

    signs = [p['sign'] for p in pulses]
    if signs != [1, 1, -1, -1]:
        print(f"Warning: '{label_str}' - pulse polarity order {signs} is not "
              f"the expected [+, +, -, -] (P, U, N, D) - skipping.")
        return None

    return pulses[0], pulses[1], pulses[2], pulses[3]


def _pund_branch(p_switch, p_nonsw, electrode_area: float):
    """Compute one half of the PUND loop from a (switching, non-switching) pair.

    The switching and non-switching pulses are, by construction, the *same*
    waveform applied twice.  They are therefore subtracted **pointwise in
    time-within-the-pulse**, which cancels the leakage and linear-dielectric
    response exactly:

        I_sw(t) = I_P(t) - I_U(t)      (positive half)
        I_sw(t) = I_N(t) - I_D(t)      (negative half)

    Interpolation is done on the time axis, which is strictly increasing, so
    ``np.interp`` is valid here (interpolating on *voltage* would not be:
    a PUND pulse ramps up and back down, so V is not monotonic).

    Charge is then accumulated by the trapezoidal rule and normalised by the
    electrode area.  ``electrode_area`` is in m**2; the 1e2 factor converts
    C/m**2 to uC/cm**2.
    """
    t_rel_sw = p_switch['t'] - p_switch['t'][0]
    t_rel_ns = p_nonsw['t'] - p_nonsw['t'][0]

    I_ns = np.interp(t_rel_sw, t_rel_ns, p_nonsw['I'])
    V_ns = np.interp(t_rel_sw, t_rel_ns, p_nonsw['V'])

    I_sw = p_switch['I'] - I_ns
    Q_sw = cumulative_trapezoid(I_sw, t_rel_sw, initial=0.0)   # Coulomb
    P_raw = (Q_sw / electrode_area) * 1e2                      # uC/cm**2

    # Diagnostic: how well the two pulses' voltage waveforms overlay
    v_resid = float(np.max(np.abs(p_switch['V'] - V_ns)))

    return {
        'V': p_switch['V'],
        'I_sw': I_sw,
        'Q': Q_sw,
        'P_raw': P_raw,
        't_rel': t_rel_sw,
        'peak_idx': p_switch['peak_idx'],
        'v_resid': v_resid,
    }


def _branch_levels(P, V=None, rel_threshold: float = 0.05):
    """Return the (baseline, plateau) polarisation levels of one PUND branch.

    2Pr is the step in polarisation between the start of a branch and its
    saturated end.  Both levels are read at the *resting* field - the voltage
    the branch begins and ends at - because that is where remanent polarisation
    is defined.  Selecting the windows by voltage rather than by a fixed
    fraction of the record matters: the ramps occupy a large and instrument-
    dependent share of a branch, so index-based windows straddle them and bias
    the result, and the Haoran waveform rests at +0.5 V rather than at zero.

    Each level is a median over its window rather than an extremum: on a noisy
    branch ``max - min`` picks out the two worst noise excursions and overstates
    Pr. Falls back to fixed 5% / 10% windows if the voltage windows are too
    short to be meaningful.
    """
    P = np.asarray(P, dtype=float)
    n = len(P)
    if n < 6:
        return float(P[0]), float(P[-1])

    head = max(3, int(0.05 * n))
    tail = max(3, int(0.10 * n))

    if V is not None:
        V = np.asarray(V, dtype=float)
        v_rest = V[0]
        span = np.max(np.abs(V - v_rest))
        if span > 0:
            near = np.abs(V - v_rest) <= rel_threshold * span
            # Length of the run of near-resting samples at each end
            lead = int(np.argmin(near)) if not near.all() else n
            trail = int(np.argmin(near[::-1])) if not near.all() else n
            if lead >= 3:
                head = lead
            if trail >= 3:
                tail = trail

    return float(np.median(P[:head])), float(np.median(P[-tail:]))


def _pund_zero_crossings(V, P, peak_idx):
    """Voltages at which P crosses zero on the *rising* ramp of a branch.

    Restricting the search to the rising ramp (index 0 .. peak_idx) is what
    makes this robust: the hold and falling parts of the pulse sit at the
    saturated value and never approach zero.
    """
    V_r, P_r = np.asarray(V)[:peak_idx + 1], np.asarray(P)[:peak_idx + 1]
    idx = np.where(np.diff(np.sign(P_r)))[0]
    out = []
    for i in idx:
        p1, p2 = P_r[i], P_r[i + 1]
        if p1 == p2:
            continue
        out.append(V_r[i] + (0.0 - p1) * (V_r[i + 1] - V_r[i]) / (p2 - p1))
    return out


# =============================================================================
#  Transforms
# =============================================================================

def extract_pund(meas: PMUdata,
                 p_pulse: Optional[str] = None,
                 centre_branches: bool = True,
                 med_filt: Optional[int] = None,
                 sigma_filt: Optional[int] = None,
                 electrode_area: Optional[float] = None,
                 force: bool = False,
                 verbose: bool = True) -> Optional[pd.DataFrame]:
    r"""Isolate the switching polarisation of a PUND measurement.

    The PUND protocol applies pulses

        P (+V) -> U (+V) -> N (-V) -> D (-V)

    (optionally with a preset before, or a repeat P after).  *P* and *N* contain
    switching plus non-switching contributions; *U* and *D* contain only the
    non-switching response - leakage plus linear dielectric charging.  Because
    each pair is the *same* waveform applied twice, the switching component is
    isolated by subtracting them pointwise in time-within-the-pulse, and the
    polarisation follows by trapezoidal integration normalised by area:

        P(t) = (1/A) * integral I_sw(t) dt

    Each half is then re-referenced so it is symmetric about zero
    (``centre_branches``).  The cumulative integral necessarily starts at zero,
    so an uncorrected branch runs 0 -> 2Pr; subtracting half its extremum shifts
    it to -Pr -> +Pr, which is what turns the two halves into a closed,
    square-looking loop.

    Note that PUND reconstructs only the *switching* polarisation, so the
    ramp-down of each branch is flat by construction.  This is not a true P-E
    hysteresis loop - that needs a DHM or Sawyer-Tower measurement - but it does
    give the correct 2Pr and coercive voltage.

    Measurements that arrive already processed (the Haoran export) are not
    recomputed; only the derived scalars are filled in.

    Parameters
    ----------
    meas            : PMUdata with ``meas_type == 'PUND'``
    p_pulse         : 'last' | 'first' | 'both', optional
        Which positive pulse pairs with U when the waveform has two.  Defaults
        to the value recorded by the importer, else 'last'.
    centre_branches : re-reference each half about P = 0 (default True)
    med_filt        : median-filter kernel applied to I before integration
    sigma_filt      : boxcar window applied to I before integration
    electrode_area  : override the area on the measurement, in m^2
    force           : recompute even if ``meas.pund`` is already populated

    Returns
    -------
    The ``pund`` DataFrame (``Time``, ``V``, ``I_sw``, ``Q``, ``P``, ``P_raw``,
    ``branch``, ``direction``), also stored on ``meas``.  ``meas.extracted_params``
    is filled with Pr, 2Pr, Vc, Ec, imprint and Qsw.
    """
    label = meas.plot_string

    # Sources that arrive already processed are never recomputed, even with
    # force=True: their raw trace is the full waveform including any offset
    # ramp, which the pulse splitter is not meant to interpret.
    if meas.pund is not None and (not force
                                  or meas.metadata.get('pund_precomputed')):
        _pund_params(meas, verbose=verbose)
        return meas.pund

    area = electrode_area if electrode_area is not None else meas.electrode_area
    if area is None:
        print(f"Warning: '{label}' has no electrode area - cannot compute "
              f"polarisation. Pass electrode_area (m^2).")
        return None

    df = meas.raw_data
    if not {'Time', 'V', 'I'}.issubset(df.columns):
        print(f"Warning: '{label}' lacks Time/V/I columns - skipping.")
        return None

    good = np.isfinite(df['V']) & np.isfinite(df['I']) & np.isfinite(df['Time'])
    df = df[good]
    if len(df) < 16:
        print(f"Warning: '{label}' has too few valid points - skipping.")
        return None

    if p_pulse is None:
        p_pulse = meas.metadata.get('p_pulse', 'last')

    pulses = _pulses_from_frame(df, med_filt, sigma_filt)
    roles = _assign_pund_pulses(pulses, label,
                                p_pulse='first' if p_pulse == 'first' else 'last')
    if roles is None:
        return None
    pulse_P, pulse_U, pulse_N, pulse_D = roles

    branches = [('P-U', _pund_branch(pulse_P, pulse_U, area)),
                ('N-D', _pund_branch(pulse_N, pulse_D, area))]

    # With 'both', also keep the alternate positive branch for comparison
    if p_pulse == 'both' and len(pulses) == 5 and \
            [p['sign'] for p in pulses] == [1, 1, -1, -1, 1]:
        branches.append(('P(first)-U', _pund_branch(pulses[0], pulse_U, area)))

    frames = []
    for name, br in branches:
        P_raw = br['P_raw']
        if centre_branches:
            base, plateau = _branch_levels(P_raw, br['V'])
            P = P_raw - (base + plateau) / 2.0
        else:
            P = P_raw
        direction = np.where(np.arange(len(P_raw)) <= br['peak_idx'],
                             'outbound', 'return')
        frames.append(pd.DataFrame({
            'Time': br['t_rel'], 'V': br['V'], 'I_sw': br['I_sw'],
            'Q': br['Q'], 'P': P, 'P_raw': P_raw,
            'branch': name, 'direction': direction,
        }))

        v_pk = float(np.max(np.abs(br['V'])))
        if br['v_resid'] > 0.05 * v_pk:
            print(f"Note: '{label}' {name} - the switching and non-switching "
                  f"pulses differ by up to {br['v_resid']:.3g} V "
                  f"({100 * br['v_resid'] / v_pk:.0f}% of peak) once aligned in "
                  f"time. Usually sampling jitter between the two pulses on a "
                  f"steep ramp, which adds a little error to the subtraction; "
                  f"it matters more on coarsely sampled records.")

    meas.pund = pd.concat(frames, ignore_index=True)
    meas.metadata['p_pulse'] = p_pulse
    _pund_params(meas, verbose=verbose)
    return meas.pund


def _pund_params(meas: PMUdata, verbose: bool = True) -> dict:
    """Derive Pr, 2Pr, Vc, Ec and imprint from an existing ``meas.pund`` frame.

    Pr is half the step between the branch's starting level and its saturated
    plateau (see :func:`_branch_levels`), i.e. the standard Pr = Qsw / 2A.  Qsw
    is then derived from Pr and the electrode area so the two always agree.
    """
    pund = meas.pund
    if pund is None or pund.empty:
        return {}

    params = {}
    for name, key in (('P-U', '+'), ('N-D', '-')):
        br = pund[pund['branch'] == name]
        if br.empty:
            continue

        P_col = br['P_raw'] if 'P_raw' in br else br['P']
        base, plateau = _branch_levels(P_col, br['V'])
        Pr = (plateau - base) / 2.0
        params[f'Pr{key}'] = float(Pr)
        if meas.electrode_area:
            params[f'Qsw{key}'] = float(2.0 * Pr * meas.electrode_area / 1e2)

        # Vc from where the centred polarisation crosses zero on the way out
        out = br[br['direction'] == 'outbound'] if 'direction' in br else br
        if len(out) > 1:
            crossings = _pund_zero_crossings(
                out['V'].to_numpy(), out['P'].to_numpy(), len(out) - 1)
            if crossings:
                params[f'Vc{key}'] = float(max(crossings) if key == '+'
                                           else min(crossings))

    Pr_p, Pr_m = params.get('Pr+', np.nan), params.get('Pr-', np.nan)
    Vc_p, Vc_m = params.get('Vc+', np.nan), params.get('Vc-', np.nan)
    params['Pr'] = (abs(Pr_p) + abs(Pr_m)) / 2.0
    params['2Pr'] = abs(Pr_p) + abs(Pr_m)
    params['Vc'] = (abs(Vc_p) + abs(Vc_m)) / 2.0
    params['Vimprint'] = (Vc_p + Vc_m) / 2.0

    if meas.film_thickness_nm is not None:
        d = meas.film_thickness_nm
        params['Ec+'] = _to_field(Vc_p, d)
        params['Ec-'] = _to_field(Vc_m, d)
        params['Ec'] = (abs(params['Ec+']) + abs(params['Ec-'])) / 2.0
        params['Eimprint'] = (params['Ec+'] + params['Ec-']) / 2.0

    # No known ferroelectric reaches 200 uC/cm2 (BiFeO3, the largest in common
    # use, is around 100), so anything past that is an artefact rather than a
    # measurement - almost always a current range that clipped, or the wrong
    # electrode area.
    if verbose and max(abs(Pr_p), abs(Pr_m)) > 200.0:
        print(f"Warning: '{meas.plot_string}' - Pr of "
              f"{max(abs(Pr_p), abs(Pr_m)):.4g} uC/cm2 exceeds any known "
              f"ferroelectric. Check the current range was not clipped and "
              f"that the electrode area ({meas.electrode_area:.4g} m^2) is right.")

    # The two halves should mirror one another. When they do not, the usual
    # causes are that the film was not fully poled before the P pulse, or that
    # leakage is asymmetric and is not cancelling in the P-U / N-D subtraction.
    if verbose and np.isfinite(Pr_p) and np.isfinite(Pr_m):
        lo, hi = min(abs(Pr_p), abs(Pr_m)), max(abs(Pr_p), abs(Pr_m))
        if lo == 0 and hi > 0:
            print(f"Warning: '{meas.plot_string}' - one half integrated to "
                  f"zero (Pr+ = {Pr_p:.4g}, Pr- = {Pr_m:.4g} uC/cm2). This "
                  f"usually means the measurement was clipped by the current "
                  f"range.")
        elif lo > 0 and hi / lo > 1.5:
            print(f"Warning: '{meas.plot_string}' - positive and negative "
                  f"halves differ by {hi / lo:.1f}x (Pr+ = {Pr_p:.4g}, "
                  f"Pr- = {Pr_m:.4g} uC/cm2). The loop will not close. Check "
                  f"that the film is fully poled before the P pulse and that "
                  f"leakage is not dominating.")

    meas.extracted_params.update(params)
    return params


def extract_dhm(meas: PMUdata, loop: int = 1,
                force: bool = False, verbose: bool = True) -> Optional[pd.DataFrame]:
    """Assemble one hysteresis loop from a DHM measurement.

    The aixACCT DHM export records three loops per measurement together with
    both drive polarities (``V+`` and ``V-``, the latter being the negated
    drive).  Loop 1 is the main hysteresis loop - its first sample is exactly
    the reported ``Pr-`` - while loops 2 and 3 are the relaxed-remanence loops
    (``Prrel-`` and ``Prrel+``).  The drive column that belongs with the chosen
    loop is picked by correlation rather than assumed, since loop 3 is recorded
    against the opposite polarity.

    Unlike PUND, a DHM loop cannot be separated into switching and
    non-switching contributions, so no polarisation is derived here: the
    remanent and coercive values in ``extracted_params`` are the instrument's
    own, copied across from the file header.  Use the loops for shape - whether
    the response is genuinely ferroelectric rather than leaky or capacitive -
    and PUND for quantitative Pr.
    """
    if meas.dhm is not None and not force:
        return meas.dhm

    df = meas.raw_data
    p_col, i_col = f'P{loop}', f'I{loop}'
    if p_col not in df.columns:
        if verbose:
            available = [c[1:] for c in df.columns
                         if len(c) == 2 and c.startswith('P') and c[1:].isdigit()]
            print(f"Warning: '{meas.plot_string}' has no loop {loop} "
                  f"(available loops: {available or 'none'}) - skipping.")
        return None

    P = df[p_col].to_numpy(dtype=float)
    # Pick the drive polarity that this loop was measured against
    V = df['V_pos'].to_numpy(dtype=float)
    if 'V_neg' in df.columns:
        V_neg = df['V_neg'].to_numpy(dtype=float)
        if np.corrcoef(V, P)[0, 1] < np.corrcoef(V_neg, P)[0, 1]:
            V = V_neg

    meas.dhm = pd.DataFrame({
        'Time': df['Time'].to_numpy(dtype=float),
        'V': V,
        'I': df[i_col].to_numpy(dtype=float) if i_col in df.columns else np.nan,
        'P': P,
        'loop': loop,
    })
    meas.metadata['dhm_loop'] = loop

    # Carry the instrument's own scalars through as the extracted parameters
    params = {}
    for src, dst in (('Pr+', 'Pr+'), ('Pr-', 'Pr-'), ('Vc+', 'Vc+'),
                     ('Vc-', 'Vc-'), ('Psw', 'Psw'), ('Pnsw', 'Pnsw'),
                     ('dPsw', 'dPsw'), ('Pvmax+', 'Pvmax+'),
                     ('Pvmax-', 'Pvmax-'), ('Wloss', 'Wloss')):
        if src in meas.instrument_params:
            params[dst] = meas.instrument_params[src]

    if 'Pr+' in params and 'Pr-' in params:
        params['Pr'] = (abs(params['Pr+']) + abs(params['Pr-'])) / 2.0
        params['2Pr'] = abs(params['Pr+']) + abs(params['Pr-'])
    if 'Vc+' in params and 'Vc-' in params:
        params['Vc'] = (abs(params['Vc+']) + abs(params['Vc-'])) / 2.0
        params['Vimprint'] = (params['Vc+'] + params['Vc-']) / 2.0
        if meas.film_thickness_nm is not None:
            d = meas.film_thickness_nm
            params['Ec+'] = _to_field(params['Vc+'], d)
            params['Ec-'] = _to_field(params['Vc-'], d)
            params['Ec'] = (abs(params['Ec+']) + abs(params['Ec-'])) / 2.0
            params['Eimprint'] = (params['Ec+'] + params['Ec-']) / 2.0

    params['source'] = 'instrument'
    meas.extracted_params.update(params)
    return meas.dhm


def extract_fatigue(run: FatigueRun, Pr_mode: str = 'raw',
                    verbose: bool = True, **pund_kwargs) -> pd.DataFrame:
    """Build the polarisation-vs-cycles table for one fatigue run.

    ``Pr_mode='raw'`` uses the values the instrument reported in its own summary
    table.  ``Pr_mode='extracted'`` re-derives Pr from each stored loop.  The
    latter is only meaningful when the fatigue run interrogated the film with
    PUND: a DHM loop mixes switching, leakage and dielectric charge with no way
    to separate them, so for DHM runs this warns and falls back to 'raw'.
    """
    if Pr_mode not in ('raw', 'extracted'):
        raise ValueError("Pr_mode must be 'raw' or 'extracted'")

    if Pr_mode == 'extracted' and run.loop_type != 'PUND':
        if verbose:
            print(f"Note: '{run.plot_string}' interrogates the film with "
                  f"{run.loop_type}, not PUND, so polarisation cannot be "
                  f"separated from leakage - using the instrument's reported "
                  f"values instead.")
        Pr_mode = 'raw'

    if Pr_mode == 'raw':
        cols = [c for c in ('cycles', 'Pr+', 'Pr-', 'Prrel+', 'Prrel-',
                            'Psw', 'Pnsw', 'dPsw', 'Vc+', 'Vc-', 'Wloss')
                if c in run.summary.columns]
        out = run.summary[cols].copy()
    else:
        rows = []
        for sub in run.sub_runs:
            extract_pund(sub, verbose=False, **pund_kwargs)
            row = {'cycles': sub.metadata.get('cycles', np.nan)}
            row.update({k: v for k, v in sub.extracted_params.items()
                        if isinstance(v, (int, float))})
            rows.append(row)
        out = pd.DataFrame(rows)

    if 'Pr+' in out and 'Pr-' in out:
        out['2Pr'] = out['Pr+'].abs() + out['Pr-'].abs()
    out.attrs['Pr_mode'] = Pr_mode

    run.summary_extracted = out
    run.extracted_params['Pr_mode'] = Pr_mode
    if 'cycles' in out and len(out):
        run.extracted_params['cycles_max'] = float(out['cycles'].max())
        if '2Pr' in out and len(out) > 1:
            first, last = out['2Pr'].iloc[0], out['2Pr'].iloc[-1]
            if first:
                # >100% means the film woke up rather than fatigued
                run.extracted_params['2Pr_final_pct'] = float(100.0 * last / first)
    return out


# =============================================================================
#  Plots
# =============================================================================

def plot_PUND(pmu_data,
              run_nums: Optional[list] = None,
              x_lim: Optional[Tuple[float, float]] = None,
              y_lim_left: Optional[Tuple[float, float]] = None,
              y_lim_right: Optional[Tuple[float, float]] = None,
              med_filt: Optional[int] = None,
              sigma_filt: Optional[int] = None,
              display_points: bool = False,
              line_width: Optional[float] = None,
              markersize: Optional[float] = None,
              plot_key: bool = True,
              plot_title: Optional[str] = None,
              export_data: bool = False,
              output_PMU: Optional[str] = None,
              fig_name: Optional[str] = None,
              fig_format: str = 'tiff',
              plot_transparency: bool = True,
              show: bool = True) -> Figure:
    """Plot the PUND waveforms against time: voltage left, current right.

    This is the sanity check on the raw measurement - it shows the pulse train
    as applied and the current it drew, before any subtraction.
    """
    from plot_style import apply_plot_style

    fig_size = apply_plot_style(export_data=export_data)
    fig, ax_left = plt.subplots(figsize=fig_size)
    ax_right = ax_left.twinx()

    kwargs = _plot_kwargs(display_points, line_width, markersize)
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key().get(
        'color', ['C0', 'C1', 'C2', 'C3'])

    runs = _select_runs(pmu_data, run_nums, 'PUND')
    for n, meas in enumerate(runs):
        df = meas.raw_data
        if df.empty or not {'Time', 'V', 'I'}.issubset(df.columns):
            continue
        x = df['Time'].to_numpy(dtype=float)
        y_v = _smooth(df['V'].to_numpy(dtype=float), med_filt, sigma_filt)
        y_i = _smooth(df['I'].to_numpy(dtype=float), med_filt, sigma_filt)

        c_l, c_r = _twin_colours(len(runs), prop_colors, n)
        label = meas.plot_string
        ax_left.plot(x, y_v, linestyle='-', color=c_l,
                     label=rf"{label} $V$", **kwargs)
        ax_right.plot(x, y_i, linestyle='--', color=c_r,
                      label=rf"{label} $I$", **kwargs)

    ax_left.set_xlabel(r"$t$ (s)")
    ax_left.set_ylabel(r"$V$ (V)")
    ax_right.set_ylabel(r"$I$ (A)")
    if len(runs) == 1:
        _colour_right_axis(ax_right, SINGLE_RIGHT)
    _set_title(ax_left, plot_title)
    ax_left.grid(True, alpha=0.3)
    ax_right.grid(False)

    if x_lim is not None:
        ax_left.set_xlim(x_lim)
    if y_lim_left is not None:
        ax_left.set_ylim(y_lim_left)
    if y_lim_right is not None:
        ax_right.set_ylim(y_lim_right)

    if plot_key:
        h1, l1 = ax_left.get_legend_handles_labels()
        h2, l2 = ax_right.get_legend_handles_labels()
        ax_left.legend(h1 + h2, l1 + l2, loc='best')

    _finish(fig, export_data, output_PMU, fig_name, 'PUND', fig_format,
            plot_transparency, show)
    return fig


def plot_PUND_polarisation(pmu_data,
                           run_nums: Optional[list] = None,
                           electrode_area: Optional[float] = None,
                           film_thickness_nm: Optional[float] = None,
                           p_pulse: Optional[str] = None,
                           centre_branches: bool = True,
                           align_zero: bool = True,
                           x_lim: Optional[Tuple[float, float]] = None,
                           y_lim_left: Optional[Tuple[float, float]] = None,
                           y_lim_right: Optional[Tuple[float, float]] = None,
                           med_filt: Optional[int] = None,
                           sigma_filt: Optional[int] = None,
                           display_points: bool = False,
                           line_width: Optional[float] = None,
                           markersize: Optional[float] = None,
                           plot_key: bool = True,
                           plot_title: Optional[str] = None,
                           print_summary: bool = True,
                           export_data: bool = False,
                           output_PMU: Optional[str] = None,
                           fig_name: Optional[str] = None,
                           fig_format: str = 'tiff',
                           plot_transparency: bool = True,
                           show: bool = True):
    r"""Plot PUND switching polarisation and switching current against voltage.

    Runs :func:`extract_pund` on any measurement that has not been processed
    yet, so this works straight after import.  See that function for the
    physics; the short version is that the two branches are ``I_P - I_U`` and
    ``I_N - I_D`` subtracted pointwise in time, integrated to charge, and
    re-referenced about zero so that together they close into a loop.

    Parameters
    ----------
    pmu_data          : PMU container(s), or a list of PMUdata objects
    run_nums          : run selection, per object if given as a list of lists
    electrode_area    : m^2, overriding the value carried by each measurement
    film_thickness_nm : nm, overriding the value carried by each measurement;
                        enables the secondary electric-field axis and Ec
    p_pulse           : 'last' | 'first' | 'both' - which P pairs with U
    centre_branches   : re-reference each half about P = 0 (default True)
    align_zero        : put I = 0 level with P = 0 (default True)

    Returns
    -------
    (Figure, pandas.DataFrame)
        The figure and a summary table of Pr, 2Pr, Vc, Ec and imprint per run.
        With ``export_data`` the table is written alongside the figure.
    """
    from plot_style import apply_plot_style

    runs = _select_runs(pmu_data, run_nums, 'PUND')

    # The importer has already extracted every run. Only re-extract when this
    # call overrides something the transform depends on.
    reprocess = any(a is not None for a in
                    (p_pulse, med_filt, sigma_filt, electrode_area)) \
        or not centre_branches
    if film_thickness_nm is not None:
        for meas in runs:
            meas.film_thickness_nm = film_thickness_nm

    fig_size = apply_plot_style(export_data=export_data)
    fig, ax_left = plt.subplots(figsize=fig_size)
    ax_right = ax_left.twinx()

    prop_colors = plt.rcParams['axes.prop_cycle'].by_key().get(
        'color', ['C0', 'C1', 'C2', 'C3'])
    kwargs = _plot_kwargs(display_points, line_width, markersize)

    summary_rows, n_plotted = [], 0
    for meas in runs:
        pund = extract_pund(meas, p_pulse=p_pulse,
                            centre_branches=centre_branches,
                            med_filt=med_filt, sigma_filt=sigma_filt,
                            electrode_area=electrode_area,
                            force=reprocess, verbose=False)
        if pund is None or pund.empty:
            continue
        if film_thickness_nm is not None and not reprocess:
            _pund_params(meas, verbose=False)   # refresh Ec for the new thickness

        label = meas.plot_string
        c_l, c_r = _twin_colours(len(runs), prop_colors, n_plotted)
        n_plotted += 1

        # Both halves share a colour so each dataset reads as one loop; only
        # the first branch carries a legend entry.
        for k, name in enumerate(pund['branch'].unique()):
            br = pund[pund['branch'] == name]
            style = '-' if not name.startswith('P(first)') else ':'
            ax_left.plot(br['V'], br['P'], linestyle=style, color=c_l,
                         label=rf"{label} $P$" if k == 0 else None,
                         **kwargs)
            ax_right.plot(br['V'], br['I_sw'] * 1e6, linestyle='--', color=c_r,
                          alpha=0.7,
                          label=rf"{label} $I$" if k == 0 else None,
                          **kwargs)

        row = {'Dataset': label}
        row.update({k: v for k, v in meas.extracted_params.items()
                    if isinstance(v, (int, float))})
        summary_rows.append(row)

    ax_left.set_xlabel(r"$V$ (V)")
    ax_left.set_ylabel(r"$P$ ($\mu$C/cm$^2$)")
    ax_right.set_ylabel(r"$I$ ($\mu$A)")
    if len(runs) == 1:
        _colour_right_axis(ax_right, SINGLE_RIGHT)
    _set_title(ax_left, plot_title)
    ax_left.grid(True, alpha=0.3)
    ax_right.grid(False)
    ax_left.axhline(0.0, color='grey', linestyle='--', linewidth=0.5, zorder=0)
    ax_left.axvline(0.0, color='grey', linestyle='--', linewidth=0.5, zorder=0)

    if x_lim is not None:
        ax_left.set_xlim(x_lim)
    if y_lim_left is not None:
        ax_left.set_ylim(y_lim_left)
    if y_lim_right is not None:
        ax_right.set_ylim(y_lim_right)
    elif align_zero:
        _align_right_zero(ax_left, ax_right)

    _field_axis(ax_left, _first_thickness(runs))

    if plot_key:
        h1, l1 = ax_left.get_legend_handles_labels()
        h2, l2 = ax_right.get_legend_handles_labels()
        ax_left.legend(h1 + h2, l1 + l2, loc='best')

    summary_df = pd.DataFrame(summary_rows)
    if print_summary and not summary_df.empty:
        print("\nPUND extracted parameters:")
        print(summary_df.to_string(index=False,
                                   float_format=lambda v: f"{v:.4g}"))

    _finish(fig, export_data, output_PMU, fig_name, 'PUND_polarisation',
            fig_format, plot_transparency, show, extra_csv=summary_df)
    return fig, summary_df


def plot_DHM(pmu_data,
             run_nums: Optional[list] = None,
             loop: int = 1,
             film_thickness_nm: Optional[float] = None,
             plot_current: bool = True,
             align_zero: bool = True,
             annotate: bool = True,
             x_lim: Optional[Tuple[float, float]] = None,
             y_lim_left: Optional[Tuple[float, float]] = None,
             y_lim_right: Optional[Tuple[float, float]] = None,
             med_filt: Optional[int] = None,
             sigma_filt: Optional[int] = None,
             display_points: bool = False,
             line_width: Optional[float] = None,
             markersize: Optional[float] = None,
             plot_key: bool = True,
             plot_title: Optional[str] = None,
             print_summary: bool = True,
             export_data: bool = False,
             output_PMU: Optional[str] = None,
             fig_name: Optional[str] = None,
             fig_format: str = 'tiff',
             plot_transparency: bool = True,
             show: bool = True):
    r"""Plot a dynamic hysteresis (P-V) loop with its switching current.

    The convention used in the ferroelectrics literature: polarisation against
    applied voltage on the left axis, and the current on a twinned right axis so
    the switching peaks line up with the steep part of the loop.  A genuine
    ferroelectric shows two well-separated current peaks near +/-Vc; a leaky or
    purely capacitive film shows a single broad hump or none at all, which is
    what makes the paired axes worth having.  Where a film thickness is known a
    secondary electric-field axis is added along the top.

    The remanent and coercive values printed and annotated are the
    **instrument's own**, not re-derived here - a DHM loop cannot be separated
    into switching and non-switching charge, so a recomputed Pr would silently
    include leakage. Use PUND when a quantitative Pr is needed.

    ``loop`` selects which of the recorded loops to draw: 1 is the main
    hysteresis loop, 2 and 3 are the relaxed-remanence loops.  ``loop='all'``
    overlays every loop present.
    """
    from plot_style import apply_plot_style

    runs = _select_runs(pmu_data, run_nums, 'DHM')
    if film_thickness_nm is not None:
        for meas in runs:
            meas.film_thickness_nm = film_thickness_nm

    fig_size = apply_plot_style(export_data=export_data)
    fig, ax_left = plt.subplots(figsize=fig_size)
    ax_right = ax_left.twinx() if plot_current else None

    prop_colors = plt.rcParams['axes.prop_cycle'].by_key().get(
        'color', ['C0', 'C1', 'C2', 'C3'])
    kwargs = _plot_kwargs(display_points, line_width, markersize)

    # Number of curves the figure will hold, so a single one can be drawn in
    # the black/red twin-axis convention
    n_curves = len(runs) * (3 if loop == 'all' else 1)

    summary_rows, n_plotted = [], 0
    for meas in runs:
        if loop == 'all':
            loops = sorted(int(c[1:]) for c in meas.raw_data.columns
                           if len(c) == 2 and c.startswith('P') and c[1:].isdigit())
        else:
            loops = [loop]

        for lp in loops:
            # Already extracted at import; only redo it if a different loop
            # is asked for than the one that was stored.
            df = extract_dhm(meas, loop=lp, verbose=False,
                             force=meas.metadata.get('dhm_loop') != lp)
            if df is None or df.empty:
                continue

            label = meas.plot_string + (f" loop {lp}" if len(loops) > 1 else "")
            c_l, c_r = _twin_colours(n_curves, prop_colors, n_plotted)
            n_plotted += 1

            ax_left.plot(df['V'], df['P'], linestyle='-', color=c_l,
                         label=rf"{label} $P$", **kwargs)
            if ax_right is not None:
                ax_right.plot(df['V'],
                              _smooth(df['I'], med_filt, sigma_filt) * 1e6,
                              linestyle='--', color=c_r, alpha=0.7,
                              label=rf"{label} $I$", **kwargs)

            row = {'Dataset': label}
            row.update({k: v for k, v in meas.extracted_params.items()
                        if isinstance(v, (int, float))})
            summary_rows.append(row)

    ax_left.set_xlabel(r"$V$ (V)")
    ax_left.set_ylabel(r"$P$ ($\mu$C/cm$^2$)")
    _set_title(ax_left, plot_title)
    ax_left.grid(True, alpha=0.3)
    ax_left.axhline(0.0, color='grey', linestyle='--', linewidth=0.5, zorder=0)
    ax_left.axvline(0.0, color='grey', linestyle='--', linewidth=0.5, zorder=0)
    if ax_right is not None:
        ax_right.set_ylabel(r"$I$ ($\mu$A)")
        ax_right.grid(False)
        if n_curves == 1:
            _colour_right_axis(ax_right, SINGLE_RIGHT)

    if x_lim is not None:
        ax_left.set_xlim(x_lim)
    if y_lim_left is not None:
        ax_left.set_ylim(y_lim_left)
    if ax_right is not None:
        if y_lim_right is not None:
            ax_right.set_ylim(y_lim_right)
        elif align_zero:
            _align_right_zero(ax_left, ax_right)

    # Mark the instrument's Pr and Vc on the loop, for a single dataset only -
    # with several overlaid the annotations stop being readable.
    if annotate and len(summary_rows) == 1:
        params = runs[0].extracted_params if runs else {}
        for key, axis in (('Pr+', 'h'), ('Pr-', 'h'), ('Vc+', 'v'), ('Vc-', 'v')):
            val = params.get(key)
            if val is None or not np.isfinite(val):
                continue
            if axis == 'h':
                ax_left.axhline(val, color='grey', linewidth=0.5, alpha=0.6)
                ax_left.annotate(rf"$P_r^{{{key[-1]}}}$ = {val:.3g}",
                                 xy=(ax_left.get_xlim()[0], val),
                                 xytext=(4, 3), textcoords='offset points',
                                 fontsize='x-small', color='grey')
            else:
                ax_left.axvline(val, color='grey', linewidth=0.5, alpha=0.6)

    _field_axis(ax_left, _first_thickness(runs))

    if plot_key:
        h1, l1 = ax_left.get_legend_handles_labels()
        if ax_right is not None:
            h2, l2 = ax_right.get_legend_handles_labels()
            h1, l1 = h1 + h2, l1 + l2
        ax_left.legend(h1, l1, loc='best')

    summary_df = pd.DataFrame(summary_rows)
    if print_summary and not summary_df.empty:
        print("\nDHM parameters (as reported by the instrument):")
        print(summary_df.to_string(index=False,
                                   float_format=lambda v: f"{v:.4g}"))

    _finish(fig, export_data, output_PMU, fig_name, 'DHM', fig_format,
            plot_transparency, show, extra_csv=summary_df)
    return fig, summary_df


def plot_fatigue(pmu_data,
                 run_nums: Optional[list] = None,
                 Pr_mode: str = 'raw',
                 sub_plots: bool = False,
                 concatenate_stages: bool = False,
                 plot_2Pr: bool = False,
                 normalise: bool = False,
                 x_lim: Optional[Tuple[float, float]] = None,
                 y_lim: Optional[Tuple[float, float]] = None,
                 display_points: bool = True,
                 line_width: Optional[float] = None,
                 markersize: Optional[float] = None,
                 plot_key: bool = True,
                 plot_title: Optional[str] = None,
                 print_summary: bool = True,
                 export_data: bool = False,
                 output_PMU: Optional[str] = None,
                 fig_name: Optional[str] = None,
                 fig_format: str = 'tiff',
                 plot_transparency: bool = True,
                 show: bool = True,
                 sub_plot_kwargs: Optional[dict] = None):
    r"""Plot remanent polarisation against fatigue cycle count.

    ``Pr_mode='raw'`` (default) plots the values the instrument reported.
    ``Pr_mode='extracted'`` re-derives Pr from each stored loop, which is only
    valid when the fatigue run interrogates the film with PUND; for DHM runs it
    warns and falls back to 'raw'.

    ``sub_plots=True`` additionally writes one figure per cycle point - the
    actual loop measured there - into ``output_PMU``, without rendering them in
    the notebook.  That is the way to check whether an apparent Pr is a real
    switching loop or an artefact, which matters most for DHM-based fatigue.

    ``concatenate_stages=True`` treats successive runs as one continuous
    experiment, offsetting each stage's cycle axis by the previous stage's final
    count.  Off by default because aixACCT restarts the cycle count in each
    Result Table, and whether the stages are really sequential on the same
    device is a judgement about the experiment, not about the file.
    """
    from plot_style import apply_plot_style

    runs = [r for r in _select_runs(pmu_data, run_nums) if isinstance(r, FatigueRun)]
    if not runs:
        print("Warning: no fatigue runs found.")

    fig_size = apply_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    prop_colors = plt.rcParams['axes.prop_cycle'].by_key().get(
        'color', ['C0', 'C1', 'C2', 'C3'])
    kwargs = _plot_kwargs(display_points, line_width, markersize)

    summary_rows, offset = [], 0.0
    for n, run in enumerate(runs):
        # Built at import; only rebuilt if a different Pr_mode is asked for.
        table = run.summary_extracted
        if table is None or run.extracted_params.get('Pr_mode') != Pr_mode:
            table = extract_fatigue(run, Pr_mode=Pr_mode)
        if table is None or table.empty or 'cycles' not in table:
            continue

        cycles = table['cycles'].to_numpy(dtype=float) + offset
        if concatenate_stages:
            offset = float(np.nanmax(cycles))

        colour = prop_colors[n % len(prop_colors)]
        label = run.plot_string

        if plot_2Pr and '2Pr' in table:
            y = table['2Pr'].to_numpy(dtype=float)
            if normalise and len(y) and y[0]:
                y = 100.0 * y / y[0]
            ax.plot(cycles, y, linestyle='-', color=colour,
                    label=rf"{label} $2P_r$", **kwargs)
        else:
            for key, style, marker in (('Pr+', '-', 'o'), ('Pr-', '--', 's')):
                if key not in table:
                    continue
                y = table[key].to_numpy(dtype=float)
                if normalise and len(y) and y[0]:
                    y = 100.0 * y / abs(y[0])
                ax.plot(cycles, y, linestyle=style, color=colour,
                        label=rf"{label} $P_r^{{{key[-1]}}}$",
                        **{**kwargs, 'marker': marker if display_points else ''})

        row = {'Dataset': label, 'loop type': run.loop_type,
               'Pr mode': Pr_mode, 'points': len(table),
               'cycles max': float(np.nanmax(cycles))}
        row.update({k: v for k, v in run.extracted_params.items()
                    if isinstance(v, (int, float))})
        summary_rows.append(row)

        # ---- per-cycle loop figures ----
        if sub_plots:
            if output_PMU is None:
                print("Warning: sub_plots=True needs output_PMU - skipping.")
            else:
                _export_fatigue_sub_plots(run, output_PMU, fig_name, fig_format,
                                          plot_transparency,
                                          sub_plot_kwargs or {})

    ax.set_xscale('log')
    ax.set_xlabel(r"$N$ (cycles)")
    ax.set_ylabel(r"$P_r$ (\%)" if normalise
                  else r"$P_r$ ($\mu$C/cm$^2$)")
    _set_title(ax, plot_title)
    ax.grid(True, alpha=0.3)
    ax.axhline(0.0, color='grey', linestyle='--', linewidth=0.5, zorder=0)

    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)
    if plot_key:
        ax.legend(loc='best')

    summary_df = pd.DataFrame(summary_rows)
    if print_summary and not summary_df.empty:
        print("\nFatigue summary:")
        print(summary_df.to_string(index=False,
                                   float_format=lambda v: f"{v:.4g}"))

    _finish(fig, export_data, output_PMU, fig_name, 'fatigue', fig_format,
            plot_transparency, show, extra_csv=summary_df)
    return fig, summary_df


def _export_fatigue_sub_plots(run: FatigueRun, output_PMU, fig_name,
                              fig_format, plot_transparency, extra_kwargs):
    """Write one loop figure per cycle point, without rendering any of them."""
    out_dir = Path(output_PMU) / (f"{fig_name}_fatigue_loops" if fig_name
                                  else "fatigue_loops")
    out_dir.mkdir(parents=True, exist_ok=True)

    rt = run.metadata.get('result_table', run.run_number)
    written = 0
    for sub in run.sub_runs:
        cycles = sub.metadata.get('cycles', sub.run_number)
        stem = f"rt{rt}_cyc{cycles:g}" if isinstance(cycles, float) \
            else f"rt{rt}_pt{sub.run_number}"

        title = rf"{run.plot_string} --- {cycles:g} cycles" \
            if isinstance(cycles, float) else run.plot_string

        if run.loop_type == 'PUND':
            plot_PUND_polarisation(
                [sub], export_data=True, output_PMU=str(out_dir),
                fig_name=stem, fig_format=fig_format,
                plot_transparency=plot_transparency, print_summary=False,
                show=False, **extra_kwargs)
        else:
            plot_DHM(
                [sub], export_data=True, output_PMU=str(out_dir),
                fig_name=stem, fig_format=fig_format,
                plot_transparency=plot_transparency, print_summary=False,
                plot_title=title, show=False, **extra_kwargs)
        written += 1

    print(f"  wrote {written} loop figure(s) to {out_dir}")


def plot_custom_waveform(pmu_data,
                         run_nums: Optional[list] = None,
                         v_channel: str = 'V_Ch2',
                         i_channel: str = 'I_Ch1',
                         x_lim: Optional[Tuple[float, float]] = None,
                         y_lim_left: Optional[Tuple[float, float]] = None,
                         y_lim_right: Optional[Tuple[float, float]] = None,
                         med_filt: Optional[int] = None,
                         sigma_filt: Optional[int] = None,
                         display_points: bool = False,
                         line_width: Optional[float] = None,
                         markersize: Optional[float] = None,
                         plot_key: bool = True,
                         plot_title: Optional[str] = None,
                         export_data: bool = False,
                         output_PMU: Optional[str] = None,
                         fig_name: Optional[str] = None,
                         fig_format: str = 'tiff',
                         plot_transparency: bool = True,
                         show: bool = True) -> Figure:
    """Plot a custom-waveform FeFET measurement: gate V and drain I vs time.

    Channel 2 drives the gate with a pulse train while channel 1 monitors the
    drain, so the default axes are gate voltage on the left and drain current on
    the right.  ``v_channel`` / ``i_channel`` select other combinations.
    """
    from plot_style import apply_plot_style

    fig_size = apply_plot_style(export_data=export_data)
    fig, ax_left = plt.subplots(figsize=fig_size)
    ax_right = ax_left.twinx()

    kwargs = _plot_kwargs(display_points, line_width, markersize)
    prop_colors = plt.rcParams['axes.prop_cycle'].by_key().get(
        'color', ['C0', 'C1', 'C2', 'C3'])

    runs = _select_runs(pmu_data, run_nums, 'CUSTOM')
    for n, meas in enumerate(runs):
        df = meas.raw_data
        if df.empty or not {v_channel, i_channel, 'Time'}.issubset(df.columns):
            print(f"Warning: '{meas.plot_string}' lacks "
                  f"{v_channel}/{i_channel}/Time - skipping.")
            continue

        x = df['Time'].to_numpy(dtype=float)
        y_v = _smooth(df[v_channel].to_numpy(dtype=float), med_filt, sigma_filt)
        y_i = _smooth(df[i_channel].to_numpy(dtype=float), med_filt, sigma_filt)

        c_l, c_r = _twin_colours(len(runs), prop_colors, n)
        label = meas.plot_string
        ax_left.plot(x, y_v, linestyle='-', color=c_l,
                     label=rf"{label} $V$", **kwargs)
        ax_right.plot(x, y_i, linestyle='--', color=c_r, alpha=0.7,
                      label=rf"{label} $I$", **kwargs)

    ax_left.set_xlabel(r"$t$ (s)")
    ax_left.set_ylabel(rf"${_CH_SYM[v_channel]}$ (V)")
    ax_right.set_ylabel(rf"${_CH_SYM[i_channel]}$ (A)")
    if len(runs) == 1:
        _colour_right_axis(ax_right, SINGLE_RIGHT)
    _set_title(ax_left, plot_title)
    ax_left.grid(True, alpha=0.3)
    ax_right.grid(False)

    if x_lim is not None:
        ax_left.set_xlim(x_lim)
    if y_lim_left is not None:
        ax_left.set_ylim(y_lim_left)
    if y_lim_right is not None:
        ax_right.set_ylim(y_lim_right)

    if plot_key:
        h1, l1 = ax_left.get_legend_handles_labels()
        h2, l2 = ax_right.get_legend_handles_labels()
        ax_left.legend(h1 + h2, l1 + l2, loc='best')

    _finish(fig, export_data, output_PMU, fig_name, 'custom_waveform',
            fig_format, plot_transparency, show)
    return fig


def update_plot_string(pmu_objects: list) -> list:
    """Print each run's index and label so they can be relabelled by hand.

    Mirrors ``SMU_functions.update_plot_string``: set
    ``obj[run].plot_string = 'my label'`` to override.
    """
    for obj in pmu_objects:
        print(f"\n{obj!r}")
        for meas in obj:
            print(f"  [{meas.run_number}] {meas.plot_string}")
    return pmu_objects
