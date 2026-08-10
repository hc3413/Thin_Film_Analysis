#import all the libraries needed
from import_dep import *
from scipy.signal import medfilt, savgol_filter
from scipy.ndimage import uniform_filter1d


def _get_smu_run_nums(run_nums, index):
    if run_nums is None:
        return None
    # Check if run_nums is a list/tuple of lists/tuples/np.ndarray
    if len(run_nums) > 0 and isinstance(run_nums[0], (list, tuple, np.ndarray)):
        if index < len(run_nums):
            return run_nums[index]
        else:
            return [] # No run numbers for this SMU object
    return run_nums

def _get_smu_run_labels(run_labels, index, original_labels):
    if run_labels is None:
        return original_labels
    if len(run_labels) > 0 and isinstance(run_labels[0], (list, tuple, np.ndarray)):
        if index < len(run_labels):
            return run_labels[index]
        else:
            return original_labels
    return run_labels


def smooth_smu_data(smu_obj, columns=('D_I_A',), method='savgol',
                    window=11, polyorder=2):
    """
    Apply in-place smoothing to the current columns of all 3T DataFrames
    stored in an SMU object.

    Parameters
    ----------
    smu_obj : SMU
        The SMU object whose iv3T_data will be modified in place.
    columns : tuple of str
        Column names to smooth (default: drain current only).
    method : str
        'savgol'   — Savitzky–Golay (preserves peaks, good for derivatives)
        'median'   — Median filter (robust to spikes / outliers)
        'boxcar'   — Moving-average (uniform) filter
        'none'     — No smoothing (no-op, useful for toggling off)
    window : int
        Filter window length in points (must be odd for savgol / median).
    polyorder : int
        Polynomial order for Savitzky–Golay only (ignored by others).

    Notes
    -----
    Modifies smu_obj.iv3T_data in place.  Re-import the SMU object to
    restore the original unsmoothed data.
    """
    if method == 'none':
        return

    # Ensure odd window for savgol / median
    if method in ('savgol', 'median') and window % 2 == 0:
        window += 1

    for df in smu_obj.iv3T_data:
        for col in columns:
            if col not in df.columns:
                continue
            y = df[col].values.astype(float)
            if len(y) < window:
                continue
            if method == 'savgol':
                df[col] = savgol_filter(y, window, min(polyorder, window - 1))
            elif method == 'median':
                df[col] = medfilt(y, kernel_size=window)
            elif method == 'boxcar':
                df[col] = uniform_filter1d(y, size=window)
            else:
                raise ValueError(f"Unknown smoothing method: {method!r}. "
                                 f"Use 'savgol', 'median', 'boxcar', or 'none'.")


def _deduplicate_and_gradient(vg, id_, smooth_window=None):
    """
    Compute transconductance gm = dId/dVg robustly.

    1. Sort by Vg.
    2. Remove duplicate Vg values (average Id at duplicates).
    3. Compute np.gradient on the clean arrays.
    4. Optionally smooth with Savitzky-Golay.

    Returns (vg_clean, id_clean, gm) — all same length, sorted by Vg.
    """
    order = np.argsort(vg)
    vg_s = vg[order].astype(float)
    id_s = id_[order].astype(float)

    # Average Id at duplicate Vg points
    unique_vg, inverse = np.unique(vg_s, return_inverse=True)
    avg_id = np.zeros_like(unique_vg)
    counts = np.zeros_like(unique_vg)
    np.add.at(avg_id, inverse, id_s)
    np.add.at(counts, inverse, 1)
    avg_id /= counts

    vg_c = unique_vg
    id_c = avg_id

    if len(vg_c) < 3:
        return vg_c, id_c, np.zeros_like(vg_c)

    gm = np.gradient(id_c, vg_c)

    if smooth_window and smooth_window >= 3 and len(gm) >= smooth_window:
        gm = savgol_filter(gm, smooth_window, polyorder=2)

    # Replace any remaining NaN/inf with 0
    gm = np.where(np.isfinite(gm), gm, 0.0)

    return vg_c, id_c, gm


def _extract_Vth_from_gm(vg, id_, gm, baseline_current=None):
    """
    Extract threshold voltage via linear extrapolation at peak transconductance.

    Standard method (baseline_current=None):
        V_th = Vg_peak - Id_peak / gm_peak          (tangent → Id = 0)

    Baseline-corrected (baseline_current given):
        V_th = Vg_peak - (Id_peak - baseline) / gm_peak   (tangent → baseline)

    The baseline-corrected version should be used when the device has a
    leakage floor that it never drops below (e.g. depletion-mode FETs).

    Returns V_th (float) or np.nan.
    """
    if len(gm) == 0 or np.all(gm == 0):
        return np.nan

    idx = np.argmax(np.abs(gm))
    gm_peak = gm[idx]

    if abs(gm_peak) < 1e-30:
        return np.nan

    baseline = baseline_current if baseline_current is not None else 0.0
    return vg[idx] - (id_[idx] - baseline) / gm_peak


def _extract_Vth_onset(vg, id_, gm, window_frac=0.5):
    """
    Linear-extrapolation Vth for a depletion-mode HEMT / 2DEG.

    The device conducts at Vg = 0 and is depleted by negative Vg.
    Residual leakage means Id never reaches zero, so the extraction is:

    1.  Estimate the leakage floor (off-state baseline) from the most
        negative Vg points where the channel is fully depleted.
    2.  Subtract baseline → ΔId = Id − baseline  (gate-modulated current).
    3.  Find peak |gm| and select a fitting window where
        |gm| ≥ window_frac × |gm_peak|.
    4.  Fit a line  ΔId = a·Vg + b  through that window.
    5.  Vth = −b/a  (gate voltage where the channel current extrapolates
        to zero, i.e. the pinch-off / threshold voltage).

    Parameters
    ----------
    vg, id_ : 1-D arrays, sorted by Vg (ascending)
    gm      : transconductance array (same length)
    window_frac : fraction of peak |gm| used to define the fitting window
                  (default 0.5 → use all points where |gm| ≥ 50 % of peak)

    Returns
    -------
    Vth      : float  — threshold voltage (V), or np.nan
    baseline : float  — estimated baseline current (A)
    """
    if len(vg) < 5:
        return np.nan, 0.0

    baseline = _estimate_baseline_current(vg, id_, gm)
    delta_id = id_ - baseline

    abs_gm = np.abs(gm)
    gm_peak = np.max(abs_gm)

    if gm_peak < 1e-30:
        return np.nan, baseline

    # Select a contiguous window around peak gm
    above = abs_gm >= window_frac * gm_peak
    idx_peak = np.argmax(abs_gm)

    # Expand left from peak while above threshold
    lo = idx_peak
    while lo > 0 and above[lo - 1]:
        lo -= 1
    # Expand right from peak while above threshold
    hi = idx_peak
    while hi < len(above) - 1 and above[hi + 1]:
        hi += 1

    if hi - lo + 1 < 2:
        # Fall back to single-point tangent
        if abs(gm[idx_peak]) < 1e-30:
            return np.nan, baseline
        return vg[idx_peak] - delta_id[idx_peak] / gm[idx_peak], baseline

    mask = slice(lo, hi + 1)
    vg_win = vg[mask]
    did_win = delta_id[mask]

    # Fit line: ΔId = a * Vg + b
    a, b = np.polyfit(vg_win, did_win, 1)

    if abs(a) < 1e-30:
        return np.nan, baseline

    Vth = -b / a
    return Vth, baseline


def _estimate_baseline_current(vg, id_, gm, n_edge=5):
    """
    Estimate the leakage-floor (off-state) current for a depletion-mode
    HEMT half-sweep.

    The device is ON at Vg = 0 and depleted by negative Vg.  Data are
    sorted Vg-ascending, so the first few points (most negative Vg) are
    the depleted end.  The mean of those points gives the residual
    leakage that cannot be gated away.

    Returns baseline (float) or 0.0 if it cannot be determined.
    """
    if len(id_) < 2:
        return 0.0

    n = min(n_edge, len(id_))
    return np.mean(id_[:n])


def _extract_SS(vg, id_):
    """
    Extract subthreshold swing SS = min(dVg / d(log10|Id|)) in mV/dec.

    Only considers regions where |Id| > 0 and the slope is positive
    (increasing current with increasing voltage).
    """
    abs_id = np.abs(id_)
    ok = abs_id > 0
    if np.sum(ok) < 5:
        return np.nan

    log_id = np.log10(abs_id[ok])
    vg_ok = vg[ok]

    # Remove duplicate log_id values to avoid gradient divide-by-zero
    unique_log, inv = np.unique(log_id, return_inverse=True)
    if len(unique_log) < 3:
        return np.nan
    avg_vg = np.zeros_like(unique_log)
    counts = np.zeros_like(unique_log)
    np.add.at(avg_vg, inv, vg_ok)
    np.add.at(counts, inv, 1)
    avg_vg /= counts

    dVg_dlog = np.gradient(avg_vg, unique_log)
    dVg_dlog = dVg_dlog[np.isfinite(dVg_dlog)]
    pos = dVg_dlog[dVg_dlog > 0]

    return np.min(pos) * 1e3 if len(pos) > 0 else np.nan

def auto_scale_current(current_data):
    """
    Automatically scale current data to appropriate units (A, mA, μA, nA).
    
    Parameters:
        current_data : array-like
            Current data in Amperes
    
    Returns:
        tuple: (scaled_data, unit_string, scale_factor)
    """
    max_abs_current = np.max(np.abs(current_data))
    
    if max_abs_current >= 1:        # >= 1 A
        return current_data, "A", 1
    elif max_abs_current >= 1e-3:   # >= 1 mA
        return current_data * 1e3, "mA", 1e3
    elif max_abs_current >= 1e-6:   # >= 1 μA
        return current_data * 1e6, r"$\mu$A", 1e6
    elif max_abs_current >= 1e-9:   # >= 1 nA
        return current_data * 1e9, "nA", 1e9
    else:  # < 1 nA
        return current_data * 1e12, "pA", 1e12

def auto_scale_resistance(resistance_data):
    """
    Automatically scale resistance data to appropriate units (Ω, kΩ, MΩ, GΩ).
    
    Parameters:
        resistance_data : array-like
            Resistance data in Ohms
    
    Returns:
        tuple: (scaled_data, unit_string, scale_factor)
    """
    max_abs_resistance = np.max(np.abs(resistance_data[np.isfinite(resistance_data)]))
    
    if max_abs_resistance >= 1e9:  # >= 1 GΩ
        return resistance_data / 1e9, r"G$\Omega$", 1e-9
    elif max_abs_resistance >= 1e6:  # >= 1 MΩ
        return resistance_data / 1e6, r"M$\Omega$", 1e-6
    elif max_abs_resistance >= 1e3:  # >= 1 kΩ
        return resistance_data / 1e3, r"k$\Omega$", 1e-3
    else:  # < 1 kΩ
        return resistance_data, r"$\Omega$", 1

def plot_IV(smu_data: list, 
           run_nums: Optional[list] = None,
           x_lim: Optional[Tuple[float, float]] = None,
           y_lim: Optional[Tuple[float, float]] = None,
           plot_mode: str = 'linear',
           plot_arrows: bool = False,
           display_points: bool = False,
           line_width: Optional[float] = None,
           markersize: Optional[float] = None,
           plot_key: bool = True,
           export_data: bool = False,
           output_SMU: Optional[str] = None,
           fig_name: Optional[str] = None,
           fig_format: str = 'tiff',
           plot_title: Optional[str] = None,
           plot_transparency: bool = True,
           show: bool = True) -> Figure:
    """
    Plot I-V characteristics with different conduction mechanism analysis modes.
    
    Supported plot modes for identifying conduction mechanisms in thin-film /
    memristive devices:
    
    'linear' (default):
        Standard I vs V plot.
        I(V) = V / R  (Ohmic conduction)
        Useful for basic characterisation and identifying linear/ohmic regions.

    'log':
        Plots log10|I| vs V.
        Useful for viewing overall leakage currents across a wide dynamic range.
    
    'SCLC' (Space-Charge Limited Current):
        Plots log10|I| vs log10|V|.
        I(V) = (9/8) * epsilon_0 * epsilon_r * mu * A / d^3  *  V^alpha
        The slope alpha on the log-log plot identifies the regime:
          alpha ~ 1 : Ohmic regime (I ~ V), thermally generated carriers dominate.
          alpha ~ 2 : Child's law / trap-free SCLC (I ~ V^2), injected carriers dominate.
          alpha > 2 : Trap-filled limit (TFL), steep rise as traps saturate.
          alpha ~ 2 (again) : Trap-saturated SCLC, all traps filled.
        Standard analysis for oxide thin films and memristive devices.
    
    'schottky' (Schottky / thermionic emission):
        Plots ln|I| vs sqrt|V|.
        I(V) = A* T^2  exp[ -q (phi_B - sqrt(qV / (4 pi epsilon_0 epsilon_r d))) / kT ]
        Linearised: ln(I) = ln(I_0) + (q/kT) sqrt(q / (4 pi eps_0 eps_r d)) * sqrt(V)
        A straight line indicates thermionic emission over the metal-insulator
        interface barrier. The slope yields the dynamic permittivity and the
        intercept gives the barrier height phi_B.
    
    'PF' (Poole-Frenkel emission):
        Plots ln(|I|/|V|) vs sqrt|V|.
        I(V) = (q mu N_c A / d) V exp[ -q(phi_T - sqrt(qV / (pi eps_0 eps_r d))) / kT ]
        Linearised: ln(I/V) = C + (q/kT) sqrt(q / (pi eps_0 eps_r d)) * sqrt(V)
        A straight line indicates bulk-limited Poole-Frenkel emission:
        field-enhanced thermal detrapping from Coulombic trap centres.
        Note: beta_PF = 2 * beta_Schottky, which distinguishes PF from
        Schottky emission via the slope. 
    
    'FN' (Fowler-Nordheim tunnelling):
        Plots ln(|I|/V^2) vs 1/|V|.
        I(V) = (q^3 A / (8 pi h phi_B d^2)) V^2 exp[-8 pi sqrt(2 m*) phi_B^(3/2) d / (3 q h V)]
        Linearised: ln(I/V^2) = C - [8 pi sqrt(2 m*) phi_B^(3/2) d / (3 q h)] * (1/V)
        A straight line with negative slope indicates Fowler-Nordheim quantum
        mechanical tunnelling through a triangular barrier at high fields.
        The slope gives the barrier height phi_B if the film thickness d is known.
    
    'TAT' (Trap-Assisted Tunnelling):
        Plots ln(|I|/V^2) vs 1/|V|  (same axes as FN).
        I(V) ~ V^2 exp[-8 pi sqrt(2 m*) phi_T^(3/2) d_trap / (3 q h V)]
        Same linearised form as FN but with trap depth phi_T replacing the
        full barrier height phi_B, and inter-trap distance d_trap replacing
        film thickness d. Distinguished from FN by: (i) lower extracted
        barrier (trap depth vs full barrier), (ii) weaker temperature
        dependence, (iii) dominance at moderate fields in defect-rich films.
        Particularly important for VCM memristors where oxygen-vacancy
        filaments provide trap-assisted tunnelling pathways in the
        high-resistance state (HRS).
    
    Parameters:
        smu_data          : List of SMU objects containing IV data
        run_nums          : Optional list of run numbers to plot. If None, plots all.
        x_lim             : Optional tuple specifying the x-axis limits (min, max)
                            applied to raw voltage before transformation.
        y_lim             : Optional tuple specifying the y-axis limits (min, max)
                            in transformed coordinates.
        plot_mode         : Conduction mechanism analysis mode. One of:
                            'linear', 'log', 'SCLC', 'schottky', 'PF', 'FN', 'TAT'.
                            Default is 'linear'.
        plot_arrows       : Whether to add directional arrows showing sweep direction.
        line_width        : Optional line width for the plot. If None, uses default.
        markersize        : Optional marker size for the plot. If None, uses default.
        plot_key          : Whether to show the legend/key. Default is True.
        export_data       : Whether to use export styling and generate output path.
        output_SMU        : Directory path for saving figures when export_data=True.
        fig_name          : Optional prefix for the filename. Output file will be
                            '{fig_name}_IV_{plot_mode}.{fig_format}' if provided,
                            or 'IV_{plot_mode}.{fig_format}' if None.
        fig_format        : Format for saved figure ('png', 'svg', 'tiff', etc.).
        plot_transparency : Whether to save with transparent background.
    
    Returns:
        The matplotlib Figure object.
    
    Notes:
        For non-linear plot modes, absolute values of voltage and current are
        used. Points where V=0 or I=0 are excluded as required by the
        transformation. Use x_lim to select a single voltage polarity for
        cleaner analysis.
        
        Additional conduction mechanisms relevant to VCM memristive devices
        that require temperature-dependent measurements to identify:
          - Ionic conduction: ln(I/V) vs 1/T  (Arrhenius, oxygen vacancy migration)
          - Variable-range hopping: ln(sigma) vs T^(-1/4) (Mott) or T^(-1/2) (Efros-Shklovskii)
          - Nearest-neighbour hopping: ln(sigma) vs 1/T
    """
    
    valid_modes = {'linear', 'log', 'SCLC', 'schottky', 'PF', 'FN', 'TAT'}
    if plot_mode not in valid_modes:
        raise ValueError(f"plot_mode must be one of {valid_modes}, got '{plot_mode}'")
    
    # Import here to avoid circular imports
    from plot_style import apply_plot_style
    
    # Set plot style based on export_data
    fig_size = apply_plot_style(export_data=export_data)
    
    fig, ax = plt.subplots(figsize=fig_size)
    
    # For 'linear' mode, auto-scale current units; analysis modes use raw Amperes
    current_unit = "A"
    current_scale = 1
    if plot_mode == 'linear':
        all_current_data = []
        for idx, smu_obj in enumerate(smu_data):
            current_run_nums = _get_smu_run_nums(run_nums, idx)
            iv_data, plot_strings = smu_obj.get_iv_data(current_run_nums)
            for df, label_str in zip(iv_data, plot_strings):
                if not df.empty:
                    current_data = df['I_A']
                    if x_lim is not None:
                        x = df['V_V']
                        mask = (x >= x_lim[0]) & (x <= x_lim[1])
                        current_data = current_data[mask]
                    all_current_data.extend(current_data.values)
        if all_current_data:
            _, current_unit, current_scale = auto_scale_current(np.array(all_current_data))
    
    # Now plot with the determined scaling
    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        iv_data, plot_strings = smu_obj.get_iv_data(current_run_nums)
        
        for df, label_str in zip(iv_data, plot_strings):
            if not df.empty:
                x_raw = df['V_V'].copy()
                y_raw = df['I_A'].copy()
                
                # Apply x_lim masking on raw voltage before transformation
                if x_lim is not None:
                    mask = (x_raw >= x_lim[0]) & (x_raw <= x_lim[1])
                    x_raw = x_raw[mask]
                    y_raw = y_raw[mask]
                    if len(x_raw) == 0:
                        continue
                
                # Transform data based on plot_mode
                if plot_mode == 'linear':
                    x = x_raw
                    y = y_raw * current_scale
                elif plot_mode == 'log':
                    valid = np.abs(y_raw) > 0
                    x = x_raw[valid]
                    y = np.abs(y_raw[valid])
                elif plot_mode == 'SCLC':
                    valid = (np.abs(x_raw) > 0) & (np.abs(y_raw) > 0)
                    x = np.log10(np.abs(x_raw[valid]))
                    y = np.log10(np.abs(y_raw[valid]))
                elif plot_mode == 'schottky':
                    valid = np.abs(y_raw) > 0
                    x = np.sqrt(np.abs(x_raw[valid]))
                    y = np.log(np.abs(y_raw[valid]))
                elif plot_mode == 'PF':
                    valid = (np.abs(x_raw) > 0) & (np.abs(y_raw) > 0)
                    x = np.sqrt(np.abs(x_raw[valid]))
                    y = np.log(np.abs(y_raw[valid]) / np.abs(x_raw[valid]))
                elif plot_mode in ('FN', 'TAT'):
                    valid = (np.abs(x_raw) > 0) & (np.abs(y_raw) > 0)
                    x = 1.0 / np.abs(x_raw[valid])
                    y = np.log(np.abs(y_raw[valid]) / x_raw[valid]**2)
                
                if len(x) == 0:
                    continue
                
                # Plot the main line with optional line width and markersize
                plot_kwargs = {'label': label_str, 'marker': 'o' if display_points else ''}
                if line_width is not None:
                    plot_kwargs['linewidth'] = line_width
                if markersize is not None:
                    plot_kwargs['markersize'] = markersize
                
                line = ax.plot(x, y, **plot_kwargs)
                
                # Add directional arrows if requested
                if plot_arrows and len(x) > 10:  # Only add arrows if we have enough points
                    color = line[0].get_color()  # Get the color of the line
                    
                    # Apply 3-point median filter to smooth data for arrow calculation only
                    # Convert to numpy arrays for easier processing
                    x_arr = np.asarray(x)
                    y_arr = np.asarray(y)
                    x_smooth = medfilt(x_arr, kernel_size=3)
                    y_smooth = medfilt(y_arr, kernel_size=3)
                    
                    # Add arrows at several points along the curve to show direction
                    
                    # If positive + negative double loop then 4 arrows are added
                    if np.max(x_smooth) > 0 and np.min(x_smooth) < 0:
                        arrow_positions = [0.05, 0.4, 0.65, 0.9]
                    else:
                        arrow_positions = [0.15, 0.65]

                    for i, pos in enumerate(arrow_positions):
                        idx = int(pos * (len(x_smooth) - 1))
                        if idx < len(x_smooth) - 5:  # Make sure we have enough points for direction
                            # Use a larger step for direction calculation on smoothed data
                            step = min(5, len(x_smooth) - idx - 1)
                            
                            # Calculate arrow direction using smoothed data
                            dx = x_smooth[idx + step] - x_smooth[idx]
                            dy = y_smooth[idx + step] - y_smooth[idx]
                            
                            # Skip if the step is too small
                            if abs(dx) < 1e-10 and abs(dy) < 1e-10:
                                continue
                            
                            # Normalize the direction vector
                            length = np.sqrt(dx**2 + dy**2)
                            if length > 0:
                                dx_norm = dx / length
                                dy_norm = dy / length
                                
                                # Scale arrow size based on data range
                                x_range = np.ptp(x_arr)
                                y_range = np.ptp(y_arr)
                                arrow_length = min(x_range, y_range) * 0.05  # 5% of data range
                                
                                # Position the arrow using original (unsmoothed) data
                                arrow_x = x_arr[idx]
                                arrow_y = y_arr[idx]
                                
                                # Use different arrow style for first arrow (start direction indicator)
                                if i == 0:  # First arrow - 1.7x larger to show starting direction
                                    ax.annotate('', 
                                              xy=(arrow_x + dx_norm * arrow_length, 
                                                  arrow_y + dy_norm * arrow_length), 
                                              xytext=(arrow_x, arrow_y),
                                              arrowprops=dict(arrowstyle='->', 
                                                            color=color, 
                                                            lw=2,  # Same line thickness as regular arrows
                                                            alpha=0.8,  # Same alpha as regular arrows
                                                            mutation_scale=int(15 * 1.7)))  # 1.7x larger arrowhead
                                else:  # Regular arrows for direction
                                    ax.annotate('', 
                                              xy=(arrow_x + dx_norm * arrow_length, 
                                                  arrow_y + dy_norm * arrow_length), 
                                              xytext=(arrow_x, arrow_y),
                                              arrowprops=dict(arrowstyle='->', 
                                                            color=color, 
                                                            lw=2, 
                                                            alpha=0.8,
                                                            mutation_scale=15))
    
    # Set axis labels and title based on plot_mode
    if plot_mode == 'linear':
        ax.set_xlabel(r"Voltage (V)")
        ax.set_ylabel(rf"Current ({current_unit})")
        ax.set_title(r"I-V Characteristics")
    elif plot_mode == 'log':
        ax.set_xlabel(r"Voltage (V)")
        ax.set_ylabel(r"Current (A)")
        ax.set_title(r"I-V Characteristics")
        ax.set_yscale('log')
    elif plot_mode == 'SCLC':
        ax.set_xlabel(r"$\log_{10}|V|$")
        ax.set_ylabel(r"$\log_{10}|I|$")
        ax.set_title(r"SCLC Analysis")
    elif plot_mode == 'schottky':
        ax.set_xlabel(r"$\sqrt{|V|}$ ($\mathrm{V^{1/2}}$)")
        ax.set_ylabel(r"$\ln|I|$ (A)")
        ax.set_title(r"Schottky Emission Analysis")
    elif plot_mode == 'PF':
        ax.set_xlabel(r"$\sqrt{|V|}$ ($\mathrm{V^{1/2}}$)")
        ax.set_ylabel(r"$\ln(|I|/|V|)$ (A/V)")
        ax.set_title(r"Poole-Frenkel Analysis")
    elif plot_mode in ('FN', 'TAT'):
        ax.set_xlabel(r"$1/|V|$ ($\mathrm{V^{-1}}$)")
        ax.set_ylabel(r"$\ln(|I|/V^2)$ ($\mathrm{A/V^2}$)")
        title = r"Fowler-Nordheim Analysis" if plot_mode == 'FN' else r"Trap-Assisted Tunnelling Analysis"
        ax.set_title(title)
    
    # Only apply y_lim if explicitly specified
    if y_lim is not None:
        ax.set_ylim(y_lim)
    
    # Show legend/key only if requested
    if plot_key:
        ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Save figure if export_data is True and output directory is provided
    if export_data and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        fname = f"{fig_name}_IV_{plot_mode}.{fig_format}" if fig_name else f"IV_{plot_mode}.{fig_format}"
        output_file = str(output_dir / fname)
        fig.savefig(output_file, dpi=600, bbox_inches='tight', 
                   transparent=plot_transparency, format=fig_format)
    
    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_CV(smu_data: list,
           run_nums: Optional[list] = None, 
           x_lim: Optional[Tuple[float, float]] = None,
           y_lim: Optional[Tuple[float, float]] = None,
           log_plot: bool = False,
           plot_Cp: bool = True,
           plot_Gp: bool = False,
           display_points: bool = False,
           line_width: Optional[float] = None,
           markersize: Optional[float] = None,
           plot_key: bool = True,
           export_data: bool = False,
           output_SMU: Optional[str] = None,
           fig_name: Optional[str] = None,
           fig_format: str = 'tiff',
           plot_title: Optional[str] = None,
           plot_transparency: bool = True,
           show: bool = True) -> Figure:
    """
    Plot C-V (capacitance-voltage) characteristics from SMU data.
    
    Parameters:
        smu_data          : List of SMU objects containing CV data
        run_nums          : Optional list of run numbers to plot. If None, plots all.
        x_lim             : Optional tuple specifying the x-axis limits (min, max).
        y_lim             : Optional tuple specifying the y-axis limits (min, max).
        log_plot          : Whether to use logarithmic scale for y-axis.
        plot_Cp           : Whether to plot capacitance (Cp).
        plot_Gp           : Whether to plot conductance (Gp).
        line_width        : Optional line width for the plot. If None, uses default.
        markersize        : Optional marker size for the plot. If None, uses rcParams or default.
        plot_key          : Whether to show the legend/key. Default is True.
        export_data       : Whether to use export styling and generate output path.
        output_SMU        : Directory path for saving figures when export_data=True.
        fig_name          : Optional prefix for the filename. Output file will be
                            '{fig_name}_CV.{fig_format}' if provided,
                            or 'CV.{fig_format}' if None.
        fig_format        : Format for saved figure ('png', 'svg', 'tiff', etc.).
        plot_transparency : Whether to save with transparent background.
    
    Returns:
        The matplotlib Figure object.
    """
    
    # Import here to avoid circular imports
    from plot_style import apply_plot_style
    
    # Set plot style based on export_data
    fig_size = apply_plot_style(export_data=export_data)
    
    fig, ax = plt.subplots(figsize=fig_size)
    
    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        cv_data, plot_strings = smu_obj.get_cv_data(current_run_nums)
        
        for df, label_str in zip(cv_data, plot_strings):
            if not df.empty:
                x = df['V_DC']
                
                # Apply x_lim masking if specified to rescale y-axis appropriately
                if x_lim is not None:
                    mask = (x >= x_lim[0]) & (x <= x_lim[1])
                    x = x[mask]
                    
                    # Skip if no data points in the specified range
                    if len(x) == 0:
                        continue
                else:
                    # Create a mask that includes all data if no x_lim specified
                    mask = x.notna()  # Just use all non-NaN values
                
                if plot_Cp:
                    y_cp = df['Cp'][mask] if x_lim is not None else df['Cp']
                    plot_kwargs_cp = {'label': rf"{label_str} ($C_p$)", 'marker': 'o' if display_points else ''}
                    if line_width is not None:
                        plot_kwargs_cp['linewidth'] = line_width
                    if markersize is not None:
                        plot_kwargs_cp['markersize'] = markersize
                    ax.plot(x, y_cp, **plot_kwargs_cp)
                
                if plot_Gp:
                    y_gp = df['Gp'][mask] if x_lim is not None else df['Gp']
                    plot_kwargs_gp = {'label': rf"{label_str} ($G_p$)", 'marker': 's' if display_points else '', 'linestyle': '--'}
                    if line_width is not None:
                        plot_kwargs_gp['linewidth'] = line_width
                    if markersize is not None:
                        plot_kwargs_gp['markersize'] = markersize
                    ax.plot(x, y_gp, **plot_kwargs_gp)
    
    ax.set_xlabel(r"DC Voltage (V)")
    
    if plot_Cp and plot_Gp:
        ax.set_ylabel(r"Capacitance (F) / Conductance (S)")
    elif plot_Cp:
        ax.set_ylabel(r"Capacitance (F)")
    else:
        ax.set_ylabel(r"Conductance (S)")
        
    ax.set_title(r"C-V Characteristics")
    
    if log_plot:
        ax.set_yscale('log')
    
    # x_lim is already applied through data masking, so y-axis will auto-scale
    # Only apply y_lim if explicitly specified
    if y_lim is not None:
        ax.set_ylim(y_lim)
    
    # Show legend/key only if requested
    if plot_key:
        ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Save figure if export_data is True and output directory is provided
    if export_data and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        fname = f"{fig_name}_CV.{fig_format}" if fig_name else f"CV.{fig_format}"
        output_file = str(output_dir / fname)
        fig.savefig(output_file, dpi=600, bbox_inches='tight', 
                   transparent=plot_transparency, format=fig_format)
    
    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_permittivity(smu_data: list,
                      vac_cap: float,
                      run_nums: Optional[list] = None,
                      x_lim: Optional[Tuple[float, float]] = None,
                      y_lim_real: Optional[Tuple[float, float]] = None,
                      y_lim_imag: Optional[Tuple[float, float]] = None,
                      log_plot: bool = False,
                      display_points: bool = False,
                      line_width: Optional[float] = None,
                      markersize: Optional[float] = None,
                      plot_key: bool = True,
                      export_data: bool = False,
                      output_SMU: Optional[str] = None,
                      fig_name: Optional[str] = None,
                      fig_format: str = 'tiff',
                      plot_title: Optional[str] = None,
                      plot_transparency: bool = True,
                      show: bool = True) -> Tuple[Figure, Figure]:
    r"""
    Plot dielectric permittivity derived from C-V / G-V measurements.

    Given a vacuum (geometric) capacitance C_0 = \varepsilon_0 A / d the
    real and imaginary parts of the relative permittivity are computed as:

        \varepsilon'  = C_p / C_0
        \varepsilon'' = G_p / (\omega \, C_0)    where  \omega = 2\pi f

    Two separate figures are produced: one for \varepsilon' and one for
    \varepsilon'', both plotted against the DC bias voltage.

    Parameters
    ----------
    smu_data          : list
        List of SMU objects containing CV data.
    vac_cap           : float
        Vacuum capacitance C_0 = eps_0 * A / d  (in Farads).  This is the
        parallel-plate capacitance *without* the dielectric constant.
    run_nums          : list, optional
        Run numbers to plot.  If None, plots all.
    x_lim             : tuple, optional
        x-axis (voltage) limits (min, max).
    y_lim_real        : tuple, optional
        y-axis limits for the eps' figure.
    y_lim_imag        : tuple, optional
        y-axis limits for the eps'' figure.
    log_plot          : bool
        Whether to use a logarithmic y-axis.
    display_points    : bool
        Whether to show data-point markers.
    line_width        : float, optional
        Override line width.
    markersize        : float, optional
        Override marker size.
    plot_key          : bool
        Whether to show the legend.
    export_data       : bool
        Whether to use export styling and save figures.
    output_SMU        : str, optional
        Directory path for saving figures when export_data is True.
    fig_name          : str, optional
        Optional prefix for saved filenames.
    fig_format        : str
        Figure file format ('png', 'svg', 'tiff', …).
    plot_transparency : bool
        Whether to save with a transparent background.

    Returns
    -------
    (fig_real, fig_imag) : tuple of Figure
        The two matplotlib Figure objects.
    """
    from plot_style import apply_plot_style

    fig_size = apply_plot_style(export_data=export_data)

    fig_real, ax_real = plt.subplots(figsize=fig_size)
    fig_imag, ax_imag = plt.subplots(figsize=fig_size)

    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        cv_data, plot_strings = smu_obj.get_cv_data(current_run_nums)

        for df, label_str in zip(cv_data, plot_strings):
            if df.empty:
                continue

            x = df['V_DC']

            # Apply x_lim masking
            if x_lim is not None:
                mask = (x >= x_lim[0]) & (x <= x_lim[1])
                x = x[mask]
                if len(x) == 0:
                    continue
            else:
                mask = x.notna()

            cp = df['Cp'][mask]
            gp = df['Gp'][mask]
            freq = df['Frequency_Hz'][mask]

            # Compute permittivity components
            eps_real = cp / vac_cap
            omega = 2.0 * np.pi * freq
            eps_imag = gp / (omega * vac_cap)

            # --- ε' plot ---
            kwargs_real = {
                'label': rf"{label_str}",
                'marker': 'o' if display_points else '',
            }
            if line_width is not None:
                kwargs_real['linewidth'] = line_width
            if markersize is not None:
                kwargs_real['markersize'] = markersize
            ax_real.plot(x, eps_real, **kwargs_real)

            # --- ε'' plot ---
            kwargs_imag = {
                'label': rf"{label_str}",
                'marker': 'o' if display_points else '',
            }
            if line_width is not None:
                kwargs_imag['linewidth'] = line_width
            if markersize is not None:
                kwargs_imag['markersize'] = markersize
            ax_imag.plot(x, eps_imag, **kwargs_imag)

    # ---- Finalise ε' figure ----
    ax_real.set_xlabel(r"DC Voltage (V)")
    ax_real.set_ylabel(r"$\varepsilon'$")
    ax_real.set_title(r"Real Permittivity ($\varepsilon'$)")
    if log_plot:
        ax_real.set_yscale('log')
    if y_lim_real is not None:
        ax_real.set_ylim(y_lim_real)
    if plot_key:
        ax_real.legend()
    ax_real.grid(True, alpha=0.3)

    # ---- Finalise ε'' figure ----
    ax_imag.set_xlabel(r"DC Voltage (V)")
    ax_imag.set_ylabel(r"$\varepsilon''$")
    ax_imag.set_title(r"Imaginary Permittivity ($\varepsilon''$)")
    if log_plot:
        ax_imag.set_yscale('log')
    if y_lim_imag is not None:
        ax_imag.set_ylim(y_lim_imag)
    if plot_key:
        ax_imag.legend()
    ax_imag.grid(True, alpha=0.3)

    # ---- Export ----
    if export_data and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)

        fname_real = f"{fig_name}_eps_real.{fig_format}" if fig_name else f"eps_real.{fig_format}"
        fig_real.savefig(str(output_dir / fname_real), dpi=600,
                         bbox_inches='tight', transparent=plot_transparency,
                         format=fig_format)

        fname_imag = f"{fig_name}_eps_imag.{fig_format}" if fig_name else f"eps_imag.{fig_format}"
        fig_imag.savefig(str(output_dir / fname_imag), dpi=600,
                         bbox_inches='tight', transparent=plot_transparency,
                         format=fig_format)

    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax_real.set_title(plot_title if plot_title else '')
    ax_imag.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig_real)
        plt.close(fig_imag)
    return fig_real, fig_imag


def plot_resistance(smu_data: list,
                   run_nums: Optional[list] = None,
                   x_lim: Optional[Tuple[float, float]] = None,
                   y_lim: Optional[Tuple[float, float]] = None,
                   log_plot: bool = True,
                   plot_arrows: bool = False,
                   display_points: bool = False,
                   line_width: Optional[float] = None,
                   markersize: Optional[float] = None,
                   plot_key: bool = True,
                   export_data: bool = False,
                   output_SMU: Optional[str] = None,
                   fig_name: Optional[str] = None,
                   fig_format: str = 'tiff',
                   plot_title: Optional[str] = None,
                   plot_transparency: bool = True,
                   show: bool = True) -> Figure:
    """
    Plot resistance vs voltage from IV data.
    
    Parameters:
        smu_data          : List of SMU objects containing IV data
        run_nums          : Optional list of run numbers to plot. If None, plots all.
        x_lim             : Optional tuple specifying the x-axis limits (min, max).
        y_lim             : Optional tuple specifying the y-axis limits (min, max).
        log_plot          : Whether to use logarithmic scale for y-axis.
        plot_arrows       : Whether to add directional arrows showing voltage sweep direction.
        line_width        : Optional line width for the plot. If None, uses default.
        markersize        : Optional marker size for the plot. If None, uses rcParams or default.
        plot_key          : Whether to show the legend/key. Default is True.
        export_data       : Whether to use export styling and generate output path.
        output_SMU        : Directory path for saving figures when export_data=True.
        fig_name          : Optional prefix for the filename. Output file will be
                            '{fig_name}_resistance.{fig_format}' if provided,
                            or 'resistance.{fig_format}' if None.
        fig_format        : Format for saved figure ('png', 'svg', 'tiff', etc.).
        plot_transparency : Whether to save with transparent background.
    
    Returns:
        The matplotlib Figure object.
    """
    
    # Import here to avoid circular imports
    from plot_style import apply_plot_style
    
    # Set plot style based on export_data
    fig_size = apply_plot_style(export_data=export_data)
    
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Collect all resistance data to determine appropriate scaling
    all_resistance_data = []
    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        iv_data, plot_strings = smu_obj.get_iv_data(current_run_nums)
        for df, label_str in zip(iv_data, plot_strings):
            if not df.empty:
                resistance_data = df['R_Ohm']
                # Filter out infinite and NaN values
                mask = np.isfinite(resistance_data)
                resistance_clean = resistance_data[mask]
                # Apply x_lim masking if specified
                if x_lim is not None:
                    x = df['V_V'][mask]
                    x_mask = (x >= x_lim[0]) & (x <= x_lim[1])
                    resistance_clean = resistance_clean[x_mask]
                all_resistance_data.extend(resistance_clean.values)
    
    # Determine the best resistance scaling
    if all_resistance_data:
        _, resistance_unit, resistance_scale = auto_scale_resistance(np.array(all_resistance_data))
    else:
        resistance_unit = "Ω"
        resistance_scale = 1
    
    # Now plot with the determined scaling
    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        iv_data, plot_strings = smu_obj.get_iv_data(current_run_nums)
        
        for df, label_str in zip(iv_data, plot_strings):
            if not df.empty:
                x = df['V_V']
                y = df['R_Ohm'] * resistance_scale  # Apply scaling for better units
                # Filter out infinite and NaN values
                mask = np.isfinite(y)
                x_clean = x[mask]
                y_clean = y[mask]
                
                # Apply x_lim masking if specified to rescale y-axis appropriately
                if x_lim is not None:
                    x_mask = (x_clean >= x_lim[0]) & (x_clean <= x_lim[1])
                    x_clean = x_clean[x_mask]
                    y_clean = y_clean[x_mask]
                    
                    # Skip if no data points in the specified range
                    if len(x_clean) == 0:
                        continue
                
                # Plot the main line with optional line width and markersize
                plot_kwargs = {'label': label_str, 'marker': 'o' if display_points else ''}
                if line_width is not None:
                    plot_kwargs['linewidth'] = line_width
                if markersize is not None:
                    plot_kwargs['markersize'] = markersize
                
                line = ax.plot(x_clean, y_clean, **plot_kwargs)
                
                # Add directional arrows if requested
                if plot_arrows and len(x_clean) > 10:  # Only add arrows if we have enough points
                    color = line[0].get_color()  # Get the color of the line
                    
                    # Apply 5-point median filter to smooth data for arrow calculation only
                    # Convert to numpy arrays for easier processing
                    x_smooth = medfilt(x_clean.values, kernel_size=3)
                    y_smooth = medfilt(y_clean.values, kernel_size=3)
                    
                    # If positive + negative double loop then 4 arrows are added
                    if np.max(x_smooth) > 0 and np.min(x_smooth) < 0:
                        arrow_positions = [0.05, 0.4, 0.65, 0.9]
                    else:
                        arrow_positions = [0.15, 0.65]
                    
                    
                    for i, pos in enumerate(arrow_positions):
                        idx = int(pos * (len(x_smooth) - 1))
                        if idx < len(x_smooth) - 5:  # Make sure we have enough points for direction
                            # Use a larger step for direction calculation on smoothed data
                            step = min(5, len(x_smooth) - idx - 1)
                            
                            # Calculate arrow direction using smoothed data
                            dx = x_smooth[idx + step] - x_smooth[idx]
                            dy = y_smooth[idx + step] - y_smooth[idx]
                            
                            # Skip if the step is too small
                            if abs(dx) < 1e-10 and abs(dy) < 1e-10:
                                continue
                            
                            # Normalize the direction vector
                            length = np.sqrt(dx**2 + dy**2)
                            if length > 0:
                                dx_norm = dx / length
                                dy_norm = dy / length
                                
                                # Scale arrow size based on data range
                                x_range = np.ptp(x_clean)
                                y_range = np.ptp(y_clean)
                                arrow_length = min(x_range, y_range) * 0.05  # 5% of data range
                                
                                # Position the arrow using original (unsmoothed) data
                                arrow_x = x_clean.iloc[idx]
                                arrow_y = y_clean.iloc[idx]
                                
                                # Use different arrow style for first arrow (start direction indicator)
                                if i == 0:  # First arrow - 1.7x larger to show starting direction
                                    ax.annotate('', 
                                              xy=(arrow_x + dx_norm * arrow_length, 
                                                  arrow_y + dy_norm * arrow_length), 
                                              xytext=(arrow_x, arrow_y),
                                              arrowprops=dict(arrowstyle='->', 
                                                            color=color, 
                                                            lw=2,  # Same line thickness as regular arrows
                                                            alpha=0.8,  # Same alpha as regular arrows
                                                            mutation_scale=int(15 * 1.7)))  # 1.7x larger arrowhead
                                else:  # Regular arrows for direction
                                    ax.annotate('', 
                                              xy=(arrow_x + dx_norm * arrow_length, 
                                                  arrow_y + dy_norm * arrow_length), 
                                              xytext=(arrow_x, arrow_y),
                                              arrowprops=dict(arrowstyle='->', 
                                                            color=color, 
                                                            lw=2, 
                                                            alpha=0.8,
                                                            mutation_scale=15))
    
    ax.set_xlabel(r"Voltage (V)")
    ax.set_ylabel(rf"Resistance ({resistance_unit})")
    ax.set_title(r"Resistance vs Voltage")
    
    if log_plot:
        ax.set_yscale('log')
    
    # x_lim is already applied through data masking, so y-axis will auto-scale
    # Only apply y_lim if explicitly specified
    if y_lim is not None:
        ax.set_ylim(y_lim)
    
    # Show legend/key only if requested
    if plot_key:
        ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Save figure if export_data is True and output directory is provided
    if export_data and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        fname = f"{fig_name}_resistance.{fig_format}" if fig_name else f"resistance.{fig_format}"
        output_file = str(output_dir / fname)
        fig.savefig(output_file, dpi=600, bbox_inches='tight', 
                   transparent=plot_transparency, format=fig_format)
    
    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


# =============================================================================
# Terminal-to-column mapping for 3-terminal data
# =============================================================================
_TERMINAL_MAP = {
    'drain':  {'I': 'D_I_A', 'V': 'D_V_V', 'R': 'D_R_Ohm',
               'label': r'$I_D$', 'linestyle': '-'},
    'source': {'I': 'S_I_A', 'V': 'S_V_V', 'R': None,
               'label': r'$I_S$', 'linestyle': '--'},
    'gate':   {'I': 'G_I_A', 'V': 'G_V_V', 'R': 'G_R_Ohm',
               'label': r'$I_G$', 'linestyle': (0, (1, 3))},
}

_TERMINAL_X_LABELS = {
    'gate':   r"Gate Voltage $V_G$ (V)",
    'drain':  r"Drain Voltage $V_D$ (V)",
    'source': r"Source Voltage $V_S$ (V)",
}


def plot_IV_3T(smu_data: list,
               run_nums: Optional[list] = None,
               run_labels: Optional[list] = None,
               x_terminal: str = 'gate',
               terminal_display: Optional[list] = None, #'gate', 'drain', or 'source'
               plot_mode: str = 'linear',
               norm_current: bool = False,
               x_lim: Optional[Tuple[float, float]] = None,
               y_lim: Optional[Tuple[float, float]] = None,
               plot_arrows: bool = False,
               display_points: bool = False,
               line_width: Optional[float] = None,
               markersize: Optional[float] = None,
               smooth_window: Optional[int] = None,
               plot_key: bool = True,
               export_data: bool = False,
               output_SMU: Optional[str] = None,
               fig_name: Optional[str] = None,
               fig_format: str = 'tiff',
               plot_title: Optional[str] = None,
               plot_transparency: bool = True,
               show: bool = True,
               drain_correct: bool = False) -> Figure:
    """
    Plot 3-terminal FET I-V characteristics with selectable terminals and
    conduction-mechanism analysis modes.

    Line styles per terminal (within each run's colour):
      - **Drain** : solid  (primary measurement)
      - **Source** : dashed
      - **Gate**   : fine-dotted  (leakage monitor)

    Supported plot modes
    --------------------
    'linear'          : I vs V  (standard transfer / output curve)
    'log'             : |I| vs V  (log y-axis – on/off ratio, subthreshold)
    'SCLC'            : log₁₀|I|  vs  log₁₀|V|
    'schottky'        : ln|I|  vs  √|V|
    'PF'              : ln(|I|/|V|)  vs  √|V|
    'FN' / 'TAT'      : ln(|I|/V²)  vs  1/|V|
    'transconductance' : dI/dV  vs  V   (numerical derivative)
    'sqrt_Id'          : √|I|  vs  V   (saturation mobility extraction)

    Parameters
    ----------
    smu_data          : list of SMU objects containing iv3T data.
    run_nums          : run numbers to plot (None → all).
    x_terminal        : 'gate', 'drain', or 'source' (case-insensitive).
                        Selects which terminal voltage appears on the x-axis.
    terminal_display  : list of terminals whose *current* is plotted on y.
                        Default ['drain'].  E.g. ['drain', 'gate'] to monitor
                        gate leakage alongside drain current.
    plot_mode         : analysis mode (see above).
    x_lim / y_lim     : axis limits (x_lim masks raw voltage before transform).
    plot_arrows       : directional sweep arrows.
    line_width        : override line width.
    markersize        : override marker size.
    smooth_window     : Savitzky-Golay window length for 'transconductance'
                        mode.  Must be odd.  None → no smoothing.
    plot_key          : show legend.
    export_data       : use export styling & save figure.
    output_SMU        : output directory.
    fig_name          : filename prefix.
    fig_format        : figure format ('png', 'svg', 'tiff', …).
    plot_transparency : transparent background.
    drain_correct     : if True, applies correction to drain current: I_drain_corrected = I_drain + I_gate / 2. Default False.

    Returns
    -------
    matplotlib Figure.
    """

    valid_modes = {'linear', 'log', 'SCLC', 'schottky', 'PF', 'FN', 'TAT',
                   'transconductance', 'sqrt_Id'}
    if plot_mode not in valid_modes:
        raise ValueError(f"plot_mode must be one of {valid_modes}, got '{plot_mode}'")

    x_terminal = x_terminal.lower()
    if x_terminal not in _TERMINAL_MAP:
        raise ValueError(f"x_terminal must be 'gate', 'drain', or 'source', got '{x_terminal}'")

    if terminal_display is None:
        terminal_display = ['drain']
    terminal_display = [t.lower() for t in terminal_display]
    for t in terminal_display:
        if t not in _TERMINAL_MAP:
            raise ValueError(f"terminal_display entries must be 'gate', 'drain', or 'source', got '{t}'")

    from plot_style import apply_plot_style
    fig_size = apply_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    x_v_col = _TERMINAL_MAP[x_terminal]['V']

    # ---- auto-scale current for linear / log / transconductance / sqrt_Id ---
    current_unit = "A"
    current_scale = 1
    if plot_mode in ('linear', 'sqrt_Id') and not norm_current:
        all_current = []
        for idx, smu_obj in enumerate(smu_data):
            current_run_nums = _get_smu_run_nums(run_nums, idx)
            data_list, _ = smu_obj.get_iv3T_data(current_run_nums)
            for df_orig in data_list:
                if df_orig.empty:
                    continue
                df = df_orig.copy()
                if drain_correct:
                    df[_TERMINAL_MAP['drain']['I']] = df[_TERMINAL_MAP['drain']['I']] + df[_TERMINAL_MAP['gate']['I']] / 2
                for term in terminal_display:
                    if norm_current and term in ('drain', 'source'):
                        continue
                    i_col = _TERMINAL_MAP[term]['I']
                    c = df[i_col]
                    if x_lim is not None:
                        xv = df[x_v_col]
                        m = (xv >= x_lim[0]) & (xv <= x_lim[1])
                        c = c[m]
                    all_current.extend(c.values)
        if all_current:
            _, current_unit, current_scale = auto_scale_current(np.array(all_current))

    # ---- main plotting loop -------------------------------------------------
    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        
        # Override run labels if provided
        current_run_labels = _get_smu_run_labels(run_labels, idx, None)

        data_list, default_plot_strings = smu_obj.get_iv3T_data(current_run_nums)
        plot_strings = current_run_labels if current_run_labels is not None else default_plot_strings

        for df_orig, run_label in zip(data_list, plot_strings):
            if df_orig.empty:
                continue

            df = df_orig.copy()
            if drain_correct:
                df[_TERMINAL_MAP['drain']['I']] = df[_TERMINAL_MAP['drain']['I']] + df[_TERMINAL_MAP['gate']['I']] / 2

            x_voltage_raw = df[x_v_col].copy()

            for term in terminal_display:
                info = _TERMINAL_MAP[term]
                i_col = info['I']
                y_current_raw = df[i_col].copy()

                if norm_current and term in ('drain', 'source'):
                    # Find current closest to Vg=0. Usually the first point in fresh data, 
                    # but closest to 0 is safer.
                    first_idx = np.argmin(np.abs(x_voltage_raw.values))
                    i_norm = y_current_raw.iloc[first_idx] if hasattr(y_current_raw, 'iloc') else y_current_raw.values[first_idx]
                    if abs(i_norm) > 1e-15:
                        y_current_raw = y_current_raw / i_norm

                # ---- x_lim masking on raw voltage ----
                xr = x_voltage_raw.copy()
                yr = y_current_raw.copy()
                if x_lim is not None:
                    mask = (xr >= x_lim[0]) & (xr <= x_lim[1])
                    xr = xr[mask]
                    yr = yr[mask]
                    if len(xr) == 0:
                        continue

                # ---- transform data per plot_mode ----
                if plot_mode == 'linear':
                    x = xr
                    # if normalized, we don't scale it. if gate, it uses current_scale
                    y_scale = 1 if (norm_current and term in ('drain', 'source')) else current_scale
                    y = yr * y_scale
                elif plot_mode == 'log':
                    valid = np.abs(yr) > 0
                    x = xr[valid]
                    y = np.abs(yr[valid])
                elif plot_mode == 'SCLC':
                    valid = (np.abs(xr) > 0) & (np.abs(yr) > 0)
                    x = np.log10(np.abs(xr[valid]))
                    y = np.log10(np.abs(yr[valid]))
                elif plot_mode == 'schottky':
                    valid = np.abs(yr) > 0
                    x = np.sqrt(np.abs(xr[valid]))
                    y = np.log(np.abs(yr[valid]))
                elif plot_mode == 'PF':
                    valid = (np.abs(xr) > 0) & (np.abs(yr) > 0)
                    x = np.sqrt(np.abs(xr[valid]))
                    y = np.log(np.abs(yr[valid]) / np.abs(xr[valid]))
                elif plot_mode in ('FN', 'TAT'):
                    valid = (np.abs(xr) > 0) & (np.abs(yr) > 0)
                    x = 1.0 / np.abs(xr[valid])
                    y = np.log(np.abs(yr[valid]) / xr[valid]**2)
                elif plot_mode == 'transconductance':
                    x_arr = np.asarray(xr)
                    y_arr = np.asarray(yr)
                    gm = np.gradient(y_arr, x_arr)
                    if smooth_window is not None and smooth_window >= 3 and len(gm) >= smooth_window:
                        gm = savgol_filter(gm, smooth_window, polyorder=2)
                    x = xr
                    y = gm
                elif plot_mode == 'sqrt_Id':
                    valid = np.abs(yr) > 0
                    x = xr[valid]
                    y = np.sqrt(np.abs(yr[valid]) * current_scale)

                if len(x) == 0:
                    continue

                # ---- label & style ----
                label_str = rf"{run_label} {info['label']}"
                plot_kwargs = {
                    'label': label_str,
                    'marker': 'o' if display_points else '',
                    'linestyle': info['linestyle'],
                }
                if line_width is not None:
                    plot_kwargs['linewidth'] = line_width
                if markersize is not None:
                    plot_kwargs['markersize'] = markersize

                line = ax.plot(x, y, **plot_kwargs)

                # ---- directional arrows ----
                if plot_arrows and len(x) > 10:
                    color = line[0].get_color()
                    x_arr = np.asarray(x)
                    y_arr = np.asarray(y)
                    x_smooth = medfilt(x_arr, kernel_size=3)
                    y_smooth = medfilt(y_arr, kernel_size=3)

                    if np.max(x_smooth) > 0 and np.min(x_smooth) < 0:
                        arrow_positions = [0.05, 0.4, 0.65, 0.9]
                    else:
                        arrow_positions = [0.15, 0.65]

                    for i, pos in enumerate(arrow_positions):
                        idx = int(pos * (len(x_smooth) - 1))
                        if idx < len(x_smooth) - 5:
                            step = min(5, len(x_smooth) - idx - 1)
                            dx = x_smooth[idx + step] - x_smooth[idx]
                            dy = y_smooth[idx + step] - y_smooth[idx]
                            if abs(dx) < 1e-10 and abs(dy) < 1e-10:
                                continue
                            length = np.sqrt(dx**2 + dy**2)
                            if length > 0:
                                dx_n = dx / length
                                dy_n = dy / length
                                x_range = np.ptp(x_arr)
                                y_range = np.ptp(y_arr)
                                a_len = min(x_range, y_range) * 0.05
                                ax_x = x_arr[idx]
                                ay_y = y_arr[idx]
                                mscale = int(15 * 1.7) if i == 0 else 15
                                ax.annotate('',
                                    xy=(ax_x + dx_n * a_len,
                                        ay_y + dy_n * a_len),
                                    xytext=(ax_x, ay_y),
                                    arrowprops=dict(arrowstyle='->',
                                                    color=color, lw=2,
                                                    alpha=0.8,
                                                    mutation_scale=mscale))

    # ---- axis labels & title ------------------------------------------------
    x_label_base = _TERMINAL_X_LABELS[x_terminal]

    y_label_current = r"$I/I_{V_G=0}$" if norm_current else rf"Current ({current_unit})"
    y_label_log = r"$I/I_{V_G=0}$" if norm_current else r"Current $I$ (A)"

    if plot_mode == 'linear':
        ax.set_xlabel(x_label_base)
        ax.set_ylabel(y_label_current)
        ax.set_title(r"FET Transfer Characteristics")
    elif plot_mode == 'log':
        ax.set_xlabel(x_label_base)
        ax.set_ylabel(y_label_log)
        ax.set_yscale('log')
        ax.set_title(r"FET Transfer Characteristics")
    elif plot_mode == 'SCLC':
        ax.set_xlabel(r"$\log_{10}|V|$")
        ax.set_ylabel(r"$\log_{10}|I|$")
        ax.set_title(r"SCLC Analysis (3T)")
    elif plot_mode == 'schottky':
        ax.set_xlabel(r"$\sqrt{|V|}$ ($\mathrm{V^{1/2}}$)")
        ax.set_ylabel(r"$\ln|I|$ (A)")
        ax.set_title(r"Schottky Emission Analysis (3T)")
    elif plot_mode == 'PF':
        ax.set_xlabel(r"$\sqrt{|V|}$ ($\mathrm{V^{1/2}}$)")
        ax.set_ylabel(r"$\ln(|I|/|V|)$ (A/V)")
        ax.set_title(r"Poole-Frenkel Analysis (3T)")
    elif plot_mode in ('FN', 'TAT'):
        ax.set_xlabel(r"$1/|V|$ ($\mathrm{V^{-1}}$)")
        ax.set_ylabel(r"$\ln(|I|/V^2)$ ($\mathrm{A/V^2}$)")
        title = r"Fowler-Nordheim Analysis (3T)" if plot_mode == 'FN' else r"Trap-Assisted Tunnelling Analysis (3T)"
        ax.set_title(title)
    elif plot_mode == 'transconductance':
        ax.set_xlabel(x_label_base)
        ax.set_ylabel(r"Transconductance $g_m$ (S)")
        ax.set_title(r"Transconductance $g_m = dI/dV$")
    elif plot_mode == 'sqrt_Id':
        ax.set_xlabel(x_label_base)
        if current_unit == r"$\mu$A":
            ax.set_ylabel(r"$\sqrt{|I|}$ ($\sqrt{\mu\,\mathrm{A}}$)")
        else:
            ax.set_ylabel(rf"$\sqrt{{|I|}}$ ($\sqrt{{\mathrm{{{current_unit}}}}}$)")
        ax.set_title(r"$\sqrt{|I_D|}$ vs $V_G$ (Saturation Mobility)")

    if y_lim is not None:
        ax.set_ylim(y_lim)
    if plot_key:
        ax.legend()
    ax.grid(True, alpha=0.3)

    if export_data and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        fname = f"{fig_name}_IV3T_{plot_mode}.{fig_format}" if fig_name else f"IV3T_{plot_mode}.{fig_format}"
        fig.savefig(str(output_dir / fname), dpi=600, bbox_inches='tight',
                    transparent=plot_transparency, format=fig_format)

    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_resistance_3T(smu_data: list,
                       run_nums: Optional[list] = None,
                       x_terminal: str = 'gate',
                       terminal_display: Optional[list] = None,
                       x_lim: Optional[Tuple[float, float]] = None,
                       y_lim: Optional[Tuple[float, float]] = None,
                       log_plot: bool = True,
                       plot_arrows: bool = False,
                       display_points: bool = False,
                       line_width: Optional[float] = None,
                       markersize: Optional[float] = None,
                       plot_key: bool = True,
                       export_data: bool = False,
                       output_SMU: Optional[str] = None,
                       fig_name: Optional[str] = None,
                       fig_format: str = 'tiff',
                       plot_title: Optional[str] = None,
                       plot_transparency: bool = True,
                       show: bool = True) -> Figure:
    """
    Plot resistance (R = V/I) for selected terminals of 3-terminal FET data.

    Parameters match plot_IV_3T; see its docstring for shared parameters.
    """

    x_terminal = x_terminal.lower()
    if terminal_display is None:
        terminal_display = ['drain']
    terminal_display = [t.lower() for t in terminal_display]

    from plot_style import apply_plot_style
    fig_size = apply_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    x_v_col = _TERMINAL_MAP[x_terminal]['V']

    # Collect all resistance data for auto-scaling
    all_res = []
    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        data_list, _ = smu_obj.get_iv3T_data(current_run_nums)
        for df in data_list:
            if df.empty:
                continue
            for term in terminal_display:
                info = _TERMINAL_MAP[term]
                i_col = info['I']
                v_col = info['V']
                r = df[v_col] / df[i_col]
                r = r.replace([np.inf, -np.inf], np.nan)
                mask = np.isfinite(r)
                r_clean = r[mask]
                if x_lim is not None:
                    xv = df[x_v_col][mask]
                    xm = (xv >= x_lim[0]) & (xv <= x_lim[1])
                    r_clean = r_clean[xm]
                all_res.extend(r_clean.values)

    if all_res:
        _, res_unit, res_scale = auto_scale_resistance(np.array(all_res))
    else:
        res_unit, res_scale = "Ω", 1

    for idx, smu_obj in enumerate(smu_data):
        current_run_nums = _get_smu_run_nums(run_nums, idx)
        data_list, plot_strings = smu_obj.get_iv3T_data(current_run_nums)
        for df, run_label in zip(data_list, plot_strings):
            if df.empty:
                continue
            x_raw = df[x_v_col].copy()
            for term in terminal_display:
                info = _TERMINAL_MAP[term]
                r = df[info['V']] / df[info['I']]
                r = r.replace([np.inf, -np.inf], np.nan)
                y = r * res_scale
                mask = np.isfinite(y)
                xc = x_raw[mask]
                yc = y[mask]

                if x_lim is not None:
                    xm = (xc >= x_lim[0]) & (xc <= x_lim[1])
                    xc = xc[xm]
                    yc = yc[xm]
                    if len(xc) == 0:
                        continue

                label_str = rf"{run_label} {info['label'].replace('I', 'R')}"
                plot_kwargs = {'label': label_str, 'marker': 'o' if display_points else '',
                               'linestyle': info['linestyle']}
                if line_width is not None:
                    plot_kwargs['linewidth'] = line_width
                if markersize is not None:
                    plot_kwargs['markersize'] = markersize

                line = ax.plot(xc, yc, **plot_kwargs)

                # Arrows
                if plot_arrows and len(xc) > 10:
                    color = line[0].get_color()
                    x_smooth = medfilt(xc.values, kernel_size=3)
                    y_smooth = medfilt(yc.values, kernel_size=3)
                    if np.max(x_smooth) > 0 and np.min(x_smooth) < 0:
                        arrow_positions = [0.05, 0.4, 0.65, 0.9]
                    else:
                        arrow_positions = [0.15, 0.65]
                    for i, pos in enumerate(arrow_positions):
                        idx = int(pos * (len(x_smooth) - 1))
                        if idx < len(x_smooth) - 5:
                            step = min(5, len(x_smooth) - idx - 1)
                            dx = x_smooth[idx + step] - x_smooth[idx]
                            dy = y_smooth[idx + step] - y_smooth[idx]
                            if abs(dx) < 1e-10 and abs(dy) < 1e-10:
                                continue
                            length = np.sqrt(dx**2 + dy**2)
                            if length > 0:
                                dx_n = dx / length
                                dy_n = dy / length
                                a_len = min(np.ptp(xc), np.ptp(yc)) * 0.05
                                mscale = int(15 * 1.7) if i == 0 else 15
                                ax.annotate('',
                                    xy=(xc.iloc[idx] + dx_n * a_len,
                                        yc.iloc[idx] + dy_n * a_len),
                                    xytext=(xc.iloc[idx], yc.iloc[idx]),
                                    arrowprops=dict(arrowstyle='->', color=color,
                                                    lw=2, alpha=0.8,
                                                    mutation_scale=mscale))

    ax.set_xlabel(_TERMINAL_X_LABELS[x_terminal])
    ax.set_ylabel(rf"Resistance ({res_unit})")
    ax.set_title(r"Resistance vs Voltage (3T)")
    if log_plot:
        ax.set_yscale('log')
    if y_lim is not None:
        ax.set_ylim(y_lim)
    if plot_key:
        ax.legend()
    ax.grid(True, alpha=0.3)

    if export_data and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        fname = f"{fig_name}_resistance_3T.{fig_format}" if fig_name else f"resistance_3T.{fig_format}"
        fig.savefig(str(output_dir / fname), dpi=600, bbox_inches='tight',
                    transparent=plot_transparency, format=fig_format)

    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def extract_fet_params(smu_data: list,
                       run_nums: Optional[list] = None,
                       Vd: Optional[float] = None,
                       Cox: Optional[float] = None,
                       W: Optional[float] = None,
                       L: Optional[float] = None) -> list:
    """
    Extract standard FET figures-of-merit from transfer curves (Id vs Vg).

    Extracted per run (always):
      - V_th          : threshold voltage (V) — linear extrapolation from
                        the Vg at peak transconductance.
      - SS            : subthreshold swing (mV/dec) — dVg / d(log10|Id|)
                        in the steepest subthreshold region.
      - on_off_ratio  : max|Id| / min|Id| (finite values only).
      - gm_max        : peak transconductance dId/dVg  (S).

    Extracted if *all* geometry parameters are provided:
      - mu_FE         : field-effect mobility  gm_max · L / (W · Cox · Vd).

    Parameters
    ----------
    smu_data : list of SMU objects.
    run_nums : run numbers to analyse (None → all).
    Vd       : drain voltage during the transfer sweep (V).
    Cox      : gate-oxide capacitance per unit area (F/m²).
    W, L     : channel width / length (m).

    Returns
    -------
    List of dicts, one per run, each containing the extracted parameters
    and a 'run_label' key.
    """
    results = []
    for smu_obj in smu_data:
        data_list, plot_strings = smu_obj.get_iv3T_data(run_nums)
        for df, label in zip(data_list, plot_strings):
            if df.empty:
                continue

            Vg = np.asarray(df['G_V_V'])
            Id = np.asarray(df['D_I_A'])

            # --- transconductance (robust: deduplicate + smooth) ---
            vg_c, id_c, gm = _deduplicate_and_gradient(Vg, Id, smooth_window=11)
            gm_max = np.max(np.abs(gm))

            # --- threshold voltage (onset method — where current leaves baseline) ---
            V_th, _ = _extract_Vth_onset(vg_c, id_c, gm)

            # --- subthreshold swing ---
            SS = _extract_SS(vg_c, id_c)

            # --- on/off ratio ---
            abs_Id = np.abs(id_c)
            finite_Id = abs_Id[abs_Id > 0]
            on_off = np.max(finite_Id) / np.min(finite_Id) if len(finite_Id) > 1 else np.nan

            res = {
                'run_label': label,
                'V_th': V_th,
                'SS_mV_dec': SS,
                'on_off_ratio': on_off,
                'gm_max_S': gm_max,
            }

            # --- field-effect mobility ---
            if Vd is not None and Cox is not None and W is not None and L is not None:
                res['mu_FE_cm2_Vs'] = (gm_max * L / (W * Cox * abs(Vd))) * 1e4  # m²→cm²

            results.append(res)

    return results


# =============================================================================
# Hysteresis analysis utilities for 3-terminal FET data
# =============================================================================

def _split_sweep(Vg: np.ndarray):
    """
    Split a gate-voltage trace at its main turning point into forward and
    backward half-sweeps.

    The sweep is assumed to ramp to one Vg extreme, then reverse to the
    other extreme, then return.  We identify the two global extrema
    (argmax and argmin) and use the region between them as the forward
    sweep and the region after the second extreme as the backward sweep.

    Returns
    -------
    idx_turn : int
        Index of the main turning point (the second Vg extreme reached).
    fwd_slice, bwd_slice : slice objects
        Indices for forward and backward half-sweeps.
    """
    if len(Vg) < 4:
        mid = len(Vg) // 2
        return mid, slice(0, mid + 1), slice(mid, len(Vg))

    i_max = int(np.argmax(Vg))
    i_min = int(np.argmin(Vg))

    # The first extreme reached is the start of the main sweep;
    # the second extreme is the turning point.
    i_first = min(i_max, i_min)
    i_second = max(i_max, i_min)

    # Guard: if both extrema are at the same index or at boundaries
    if i_first == i_second or i_second >= len(Vg) - 1:
        mid = len(Vg) // 2
        return mid, slice(0, mid + 1), slice(mid, len(Vg))

    idx_turn = i_second
    fwd_slice = slice(i_first, i_second + 1)
    bwd_slice = slice(i_second, len(Vg))

    return idx_turn, fwd_slice, bwd_slice


def extract_hysteresis(smu_data: list,
                       run_nums: Optional[list] = None,
                       smooth_gm: int = 11) -> list:
    """
    Extract intra-cycle hysteresis parameters from 3T FET transfer curves.

    For each run the gate-voltage trace is split at its turning point into
    forward and backward half-sweeps.  Key parameters are extracted from
    each half independently, allowing quantitative comparison.

    Returned per run
    ----------------
    run_label       : str       — plot label
    run_num         : int       — run number
    V_th_fwd / bwd  : float     — threshold voltage (V)
    delta_V_th      : float     — V_th_bwd − V_th_fwd  (hysteresis window)
    SS_fwd / bwd    : float     — subthreshold swing (mV/dec)
    gm_max_fwd/bwd  : float     — peak transconductance (S)
    on_off_fwd/bwd  : float     — on/off ratio
    Q_gate          : float     — cumulative absolute gate charge ∫|I_G| dV_G (C)
    I_G_max         : float     — peak |I_G| during sweep (A)
    """
    from scipy.interpolate import interp1d

    results = []
    for smu_obj in smu_data:
        data_list, plot_strings = smu_obj.get_iv3T_data(run_nums)
        run_numbers = smu_obj.iv3T_run_nums if run_nums is None else run_nums

        for (df, label, rn) in zip(data_list, plot_strings, run_numbers):
            if df.empty:
                continue

            Vg = np.asarray(df['G_V_V'])
            Id = np.asarray(df['D_I_A'])
            Ig = np.asarray(df['G_I_A'])

            idx_turn, fwd, bwd = _split_sweep(Vg)

            res = {'run_label': label, 'run_num': rn}

            for tag, sl in [('fwd', fwd), ('bwd', bwd)]:
                vg = Vg[sl]
                id_ = Id[sl]

                if len(vg) < 5:
                    res[f'V_th_{tag}'] = np.nan
                    res[f'SS_{tag}'] = np.nan
                    res[f'gm_max_{tag}'] = np.nan
                    res[f'on_off_{tag}'] = np.nan
                    continue

                # Deduplicate, sort, and compute gm robustly
                vg_c, id_c, gm = _deduplicate_and_gradient(vg, id_, smooth_window=smooth_gm)

                if len(vg_c) < 3:
                    res[f'V_th_{tag}'] = np.nan
                    res[f'SS_{tag}'] = np.nan
                    res[f'gm_max_{tag}'] = np.nan
                    res[f'on_off_{tag}'] = np.nan
                    continue

                gm_max = np.max(np.abs(gm))
                res[f'gm_max_{tag}'] = gm_max

                # V_th via onset method (where current first rises above baseline)
                Vth_onset, baseline = _extract_Vth_onset(vg_c, id_c, gm)
                res[f'V_th_{tag}'] = Vth_onset
                res[f'baseline_{tag}'] = baseline

                # Subthreshold swing
                res[f'SS_{tag}'] = _extract_SS(vg_c, id_c)

                # On/off
                abs_id = np.abs(id_c)
                fin = abs_id[abs_id > 0]
                res[f'on_off_{tag}'] = np.max(fin) / np.min(fin) if len(fin) > 1 else np.nan

            # Hysteresis window
            res['delta_V_th'] = res.get('V_th_bwd', np.nan) - res.get('V_th_fwd', np.nan)

            # Gate charge (cumulative |Ig| integrated over Vg)
            dVg = np.abs(np.diff(Vg))
            avg_Ig = 0.5 * (np.abs(Ig[:-1]) + np.abs(Ig[1:]))
            res['Q_gate'] = np.sum(avg_Ig * dVg)
            res['I_G_max'] = np.max(np.abs(Ig))

            results.append(res)

    return results


def debug_Vth_extraction(smu_data: list,
                         run_nums: Optional[list] = None,
                         smooth_gm: int = 11,
                         window_frac: float = 0.5,
                         fig_size: list = [18, 5]):
    """
    Diagnostic plot for baseline-corrected max-slope Vth extraction.

    For every requested run, produces a 3-panel figure:
      (a) delta_Id vs Vg  — baseline-subtracted current with the fitted line
          through the peak-gm window, and Vth where that line crosses zero.
      (b) gm vs Vg  — with gm window region shaded.
      (c) Raw Id vs Vg  — original data with baseline and Vth marked.
    """
    for smu_obj in smu_data:
        data_list, plot_strings = smu_obj.get_iv3T_data(run_nums)
        run_numbers = smu_obj.iv3T_run_nums if run_nums is None else run_nums

        for (df, label, rn) in zip(data_list, plot_strings, run_numbers):
            if df.empty:
                continue

            Vg = np.asarray(df['G_V_V'])
            Id = np.asarray(df['D_I_A'])

            idx_turn, fwd, bwd = _split_sweep(Vg)

            fig, axes = plt.subplots(1, 3, figsize=fig_size)
            fig.suptitle(f'{label}  ---  Vth extraction diagnostic', fontsize=11)

            for tag, sl, color in [('fwd', fwd, 'C0'), ('bwd', bwd, 'C1')]:
                vg_raw = Vg[sl]
                id_raw = Id[sl]

                if len(vg_raw) < 5:
                    continue

                vg_c, id_c, gm = _deduplicate_and_gradient(vg_raw, id_raw,
                                                            smooth_window=smooth_gm)
                if len(vg_c) < 3:
                    continue

                # Compute Vth with the new method
                Vth, baseline = _extract_Vth_onset(vg_c, id_c, gm,
                                                    window_frac=window_frac)
                delta_id = id_c - baseline

                # Identify the contiguous gm window around peak
                abs_gm = np.abs(gm)
                gm_peak_val = np.max(abs_gm)
                above = abs_gm >= window_frac * gm_peak_val
                idx_peak = np.argmax(abs_gm)

                lo = idx_peak
                while lo > 0 and above[lo - 1]:
                    lo -= 1
                hi = idx_peak
                while hi < len(above) - 1 and above[hi + 1]:
                    hi += 1
                win_slice = slice(lo, hi + 1)

                # --- (a) delta_Id vs Vg with fitted line ---
                ax = axes[0]
                ax.plot(vg_c, delta_id * 1e6, '-', color=color,
                        label=f'{tag} $\\Delta$Id', lw=1.5)

                if hi - lo + 1 >= 2:
                    vg_win = vg_c[win_slice]
                    did_win = delta_id[win_slice]
                    a, b_coeff = np.polyfit(vg_win, did_win, 1)

                    # Draw fitted line extended to Vth
                    vg_lo = min(vg_c.min(), Vth if not np.isnan(Vth) else vg_c.min()) - 0.3
                    vg_ext = np.linspace(vg_lo, vg_c.max() + 0.3, 200)
                    ax.plot(vg_ext, (a * vg_ext + b_coeff) * 1e6, '--', color=color,
                            alpha=0.6, label=f'{tag} fit (gm window)')

                    # Shade the fitting window
                    ax.axvspan(vg_win.min(), vg_win.max(), alpha=0.08, color=color)

                ax.axhline(0, color='grey', lw=0.5)
                if not np.isnan(Vth):
                    ax.axvline(Vth, ls=':', color=color, alpha=0.6)
                    ax.plot(Vth, 0, 'x', color=color, ms=10, mew=2,
                            label=f'{tag} Vth={Vth:.2f} V')

                ax.set_xlabel(r'$V_G$ (V)')
                ax.set_ylabel(r'$\Delta I_D$ ($\mu$A)')
                ax.set_title(r'(a) Baseline-subtracted + max-slope fit')
                ax.legend(fontsize=7, loc='best')
                ax.grid(True, alpha=0.3)

                # --- (b) gm vs Vg with window ---
                ax = axes[1]
                ax.plot(vg_c, gm * 1e6, '-', color=color,
                        label=f'{tag} gm', lw=1.5)
                idx_peak = np.argmax(abs_gm)
                ax.plot(vg_c[idx_peak], gm[idx_peak] * 1e6, 'o', color=color,
                        ms=8, label=f'{tag} peak gm')
                gm_thresh = window_frac * gm_peak_val
                ax.axhline(gm_thresh * 1e6, ls='--', color=color, alpha=0.3,
                           lw=1, label=f'{tag} {window_frac:.0%} threshold')
                ax.axhline(-gm_thresh * 1e6, ls='--', color=color, alpha=0.3, lw=1)
                if hi - lo + 1 >= 2:
                    ax.axvspan(vg_c[win_slice].min(), vg_c[win_slice].max(),
                               alpha=0.08, color=color)
                if not np.isnan(Vth):
                    ax.axvline(Vth, ls=':', color=color, alpha=0.5)
                ax.set_xlabel(r'$V_G$ (V)')
                ax.set_ylabel(r'$g_m$ ($\mu$A/V)')
                ax.set_title('(b) Transconductance + fit window')
                ax.legend(fontsize=7, loc='best')
                ax.grid(True, alpha=0.3)

                # --- (c) Raw Id vs Vg with baseline + Vth ---
                ax = axes[2]
                ax.plot(vg_c, id_c * 1e6, '-', color=color,
                        label=f'{tag} Id', lw=1.5)
                ax.axhline(baseline * 1e6, ls='-', color=color, alpha=0.3,
                           lw=1, label=f'{tag} baseline={baseline*1e6:.2f} $\\mu$A')
                if not np.isnan(Vth):
                    ax.axvline(Vth, ls=':', color=color, alpha=0.6)
                    ax.plot(Vth, baseline * 1e6, 'x', color=color, ms=10, mew=2,
                            label=f'{tag} Vth={Vth:.2f} V')
                ax.set_xlabel(r'$V_G$ (V)')
                ax.set_ylabel(r'$I_D$ ($\mu$A)')
                ax.set_title('(c) Raw data + baseline + Vth')
                ax.legend(fontsize=7, loc='best')
                ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.show()


def plot_hysteresis_evolution(smu_data: list,
                              run_nums: Optional[list] = None,
                              smooth_gm: int = 11,
                              display_points: bool = True,
                              export_data: bool = False,
                              output_SMU: Optional[str] = None,
                              fig_format: str = 'tiff',
                              plot_title: Optional[str] = None,
                              plot_transparency: bool = True,
                              show: bool = True) -> Tuple[Figure, list]:
    """
    Multi-panel figure showing how hysteresis evolves across successive sweeps.

    Panel layout (2 × 2):
      (a) V_th forward & backward vs run number
      (b) ΔV_th (hysteresis window) vs run number
      (c) SS forward & backward vs run number
      (d) Cumulative gate charge Q_G vs run number

    Returns the Figure and the list of hysteresis dicts.
    """
    from plot_style import apply_plot_style
    fig_size = apply_plot_style(export_data=export_data)

    hyst = extract_hysteresis(smu_data, run_nums=run_nums, smooth_gm=smooth_gm)
    if not hyst:
        raise ValueError("No hysteresis data extracted — check run_nums")

    runs = np.array([h['run_num'] for h in hyst])
    Vth_f = np.array([h['V_th_fwd'] for h in hyst])
    Vth_b = np.array([h['V_th_bwd'] for h in hyst])
    dVth = np.array([h['delta_V_th'] for h in hyst])
    SS_f = np.array([h['SS_fwd'] for h in hyst])
    SS_b = np.array([h['SS_bwd'] for h in hyst])
    Qg = np.array([h['Q_gate'] for h in hyst])

    marker = 'o' if display_points else ''
    fig, axes = plt.subplots(2, 2, figsize=(fig_size[0] * 2, fig_size[1] * 2))

    # (a) V_th
    ax = axes[0, 0]
    ax.plot(runs, Vth_f, marker=marker, label=r'$V_{\mathrm{th}}$ fwd')
    ax.plot(runs, Vth_b, marker=marker, label=r'$V_{\mathrm{th}}$ bwd')
    ax.set_xlabel('Run number')
    ax.set_ylabel(r'$V_{\mathrm{th}}$ (V)')
    ax.set_title(r'(a) Threshold Voltage')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (b) ΔV_th
    ax = axes[0, 1]
    ax.plot(runs, dVth * 1e3, marker=marker, color='C2')  # mV
    ax.set_xlabel('Run number')
    ax.set_ylabel(r'$\Delta V_{\mathrm{th}}$ (mV)')
    ax.set_title(r'(b) Hysteresis Window')
    ax.axhline(0, color='k', lw=0.5, ls='--')
    ax.grid(True, alpha=0.3)

    # (c) SS
    ax = axes[1, 0]
    ax.plot(runs, SS_f, marker=marker, label='SS fwd')
    ax.plot(runs, SS_b, marker=marker, label='SS bwd')
    ax.set_xlabel('Run number')
    ax.set_ylabel('SS (mV/dec)')
    ax.set_title(r'(c) Subthreshold Swing')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (d) Cumulative gate charge
    ax = axes[1, 1]
    ax.plot(runs, Qg, marker=marker, color='C3')
    ax.set_xlabel('Run number')
    ax.set_ylabel(r'$Q_G = \int |I_G|\,d|V_G|$ (C)')
    ax.set_title(r'(d) Cumulative Gate Charge')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax.grid(True, alpha=0.3)

    # The (a)-(d) panel labels stay - they identify the panels rather than
    # titling the figure. plot_title adds a figure-level title if one is wanted.
    if plot_title:
        fig.suptitle(plot_title)

    fig.tight_layout()

    if export_data and output_SMU is not None:
        from pathlib import Path
        out = Path(output_SMU)
        out.mkdir(exist_ok=True)
        fig.savefig(str(out / f'hysteresis_evolution.{fig_format}'), dpi=600,
                    bbox_inches='tight', transparent=plot_transparency,
                    format=fig_format)

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig, hyst


def fit_hysteresis_model(hyst_data: list,
                         display_points: bool = True,
                         beta_limit: Union[bool, tuple] = False,
                         export_data: bool = False,
                         output_SMU: Optional[str] = None,
                         fig_format: str = 'tiff',
                         plot_transparency: bool = True) -> Tuple[Figure, dict]:
    """
    Fit competing physical models to the inter-cycle V_th drift and ΔV_th
    evolution to identify the dominant hysteresis mechanism.

    Models fitted
    -------------
    1. **Stretched exponential** (charge trapping):
       ΔV(n) = ΔV_∞ · (1 − exp[−(n/τ)^β])
       Characteristic of distributed trap time-constants at the
       semiconductor–dielectric interface or in the bulk oxide.
       β < 1 → broad distribution of trap energies.

    2. **Power-law / diffusion** (mobile ion / oxygen-vacancy drift):
       ΔV(n) = A · n^α
       α ≈ 0.5 → Fickian diffusion (e.g. oxygen vacancy redistribution),
       α ≈ 1.0 → drift-dominated (electric-field-driven ion migration).

    3. **Linear** (constant-rate degradation):
       ΔV(n) = a · n + b
       Steady defect creation or irreversible electrochemical reaction.

    The function fits all three to the *forward* V_th shift and to the
    |ΔV_th| evolution. AIC (Akaike Information Criterion) selects the best.

    Parameters
    ----------
    hyst_data : list of dicts from extract_hysteresis().
    beta_limit : False or tuple (lo, hi)
        If False (default), beta is bounded (0.01, 5.0).
        If a tuple, e.g. (0, 1), constrains beta to that range.
    display_points, export_data, output_SMU, fig_format, plot_transparency :
        usual plotting parameters.

    Returns
    -------
    fig   : 1 × 2 Figure (V_th drift + ΔV_th fits with model overlays).
    stats : dict with fit parameters, AIC, and best-model label.
    """
    from scipy.optimize import curve_fit

    from plot_style import apply_plot_style
    fig_size = apply_plot_style(export_data=export_data)

    runs = np.array([h['run_num'] for h in hyst_data], dtype=float)
    Vth_f = np.array([h['V_th_fwd'] for h in hyst_data])
    dVth = np.array([h['delta_V_th'] for h in hyst_data])

    # Use sweep index (0-based) as the "time" proxy
    n = np.arange(len(runs), dtype=float)

    # --- V_th shift from first sweep ---
    Vth_shift = Vth_f - Vth_f[0]

    # ----- model definitions -----
    def stretched_exp(n, dV_inf, tau, beta):
        return dV_inf * (1.0 - np.exp(-((n + 1) / tau) ** beta))

    def power_law(n, A, alpha):
        return A * (n + 1) ** alpha

    def linear_model(n, a, b):
        return a * n + b

    # ----- fit helper -----
    def _fit(xdata, ydata, func, p0, bounds):
        try:
            popt, pcov = curve_fit(func, xdata, ydata, p0=p0, bounds=bounds,
                                   maxfev=10000)
            residuals = ydata - func(xdata, *popt)
            ss_res = np.sum(residuals ** 2)
            k = len(popt)
            nn = len(ydata)
            aic = nn * np.log(ss_res / nn + 1e-30) + 2 * k
            return popt, aic, ss_res
        except Exception:
            return None, np.inf, np.inf

    stats = {}

    marker = 'o' if display_points else ''
    fig, axes = plt.subplots(1, 2, figsize=(fig_size[0] * 2, fig_size[1]))
    n_dense = np.linspace(0, n[-1], 200)

    for ax, ydata, ylabel, tag in [
        (axes[0], Vth_shift, r'$\Delta V_{\mathrm{th, fwd}}$ from sweep 1 (V)', 'Vth_shift'),
        (axes[1], np.abs(dVth), r'$|\Delta V_{\mathrm{th}}|$ intra-cycle (V)', 'abs_dVth'),
    ]:
        ax.plot(n, ydata, marker=marker, ls='', color='k', label='Data', zorder=5)

        y_range = np.ptp(ydata) if np.ptp(ydata) > 0 else 1e-6

        # 1. Stretched exponential
        beta_lo = beta_limit[0] if beta_limit else 0.01
        beta_hi = beta_limit[1] if beta_limit else 5.0
        p0_beta = np.clip(0.5, beta_lo, beta_hi)
        p0_se = [ydata[-1] if abs(ydata[-1]) > 1e-12 else y_range, len(n) / 2, p0_beta]
        bounds_se = ([-np.inf, 0.01, beta_lo], [np.inf, len(n) * 100, beta_hi])
        popt_se, aic_se, _ = _fit(n, ydata, stretched_exp, p0_se, bounds_se)

        # 2. Power-law
        p0_pl = [y_range, 0.5]
        bounds_pl = ([-np.inf, 0.0], [np.inf, 5.0])
        popt_pl, aic_pl, _ = _fit(n, ydata, power_law, p0_pl, bounds_pl)

        # 3. Linear
        p0_lin = [0.0, ydata[0]]
        popt_lin, aic_lin, _ = _fit(n, ydata, linear_model, p0_lin, (-np.inf, np.inf))

        # Plot fits
        fits = {}
        if popt_se is not None:
            ax.plot(n_dense, stretched_exp(n_dense, *popt_se), '--',
                    label=rf'Stretched exp ($\beta$={popt_se[2]:.2f}, $\tau$={popt_se[1]:.1f})')
            fits['stretched_exp'] = {'params': popt_se, 'AIC': aic_se,
                                     'names': ['dV_inf', 'tau', 'beta']}
        if popt_pl is not None:
            ax.plot(n_dense, power_law(n_dense, *popt_pl), '-.',
                    label=rf'Power law ($\alpha$={popt_pl[1]:.2f})')
            fits['power_law'] = {'params': popt_pl, 'AIC': aic_pl,
                                 'names': ['A', 'alpha']}
        if popt_lin is not None:
            ax.plot(n_dense, linear_model(n_dense, *popt_lin), ':',
                    label=rf'Linear (slope={popt_lin[0]:.2e})')
            fits['linear'] = {'params': popt_lin, 'AIC': aic_lin,
                              'names': ['a', 'b']}

        # Select best by AIC
        best = min(fits, key=lambda k: fits[k]['AIC']) if fits else 'none'
        stats[tag] = {'fits': fits, 'best_model': best}

        ax.set_xlabel('Sweep index')
        ax.set_ylabel(ylabel)
        ax.set_title(f'Best: {best}')
        ax.legend(fontsize=6)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if export_data and output_SMU is not None:
        from pathlib import Path
        out = Path(output_SMU)
        out.mkdir(exist_ok=True)
        fig.savefig(str(out / f'hysteresis_model_fits.{fig_format}'), dpi=600,
                    bbox_inches='tight', transparent=plot_transparency,
                    format=fig_format)

    plt.show()
    return fig, stats


def plot_charge_correlation(hyst_data: list,
                            display_points: bool = True,
                            export_data: bool = False,
                            output_SMU: Optional[str] = None,
                            fig_format: str = 'tiff',
                            plot_title: Optional[str] = None,
                            plot_transparency: bool = True,
                            show: bool = True) -> Figure:
    """
    Test the electrochemical hypothesis: if ΔV_th is driven by ionic
    migration (oxygen vacancies), it should scale with the cumulative
    charge injected through the gate, Q_G = Σ ∫|I_G| d|V_G|.

    Plots |ΔV_th| vs cumulative Q_G with a linear fit.  A strong
    correlation (R² → 1) supports an electrochemical / ion-migration
    mechanism; weak correlation favours interface-trap or bulk-trap
    mechanisms.

    Returns the Figure.
    """
    from plot_style import apply_plot_style
    fig_size = apply_plot_style(export_data=export_data)

    dVth = np.array([h['delta_V_th'] for h in hyst_data])
    Qg = np.array([h['Q_gate'] for h in hyst_data])

    # Cumulative gate charge across sweeps
    Qg_cum = np.cumsum(Qg)

    marker = 'o' if display_points else ''
    fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(Qg_cum, np.abs(dVth) * 1e3, marker=marker, ls='', color='C0', label='Data')

    # Linear fit
    if len(Qg_cum) > 2:
        coeffs = np.polyfit(Qg_cum, np.abs(dVth) * 1e3, 1)
        fit_line = np.poly1d(coeffs)
        Qg_dense = np.linspace(Qg_cum[0], Qg_cum[-1], 100)
        ax.plot(Qg_dense, fit_line(Qg_dense), '--', color='C1')

        # R²
        y_pred = fit_line(Qg_cum)
        ss_res = np.sum((np.abs(dVth) * 1e3 - y_pred) ** 2)
        ss_tot = np.sum((np.abs(dVth) * 1e3 - np.mean(np.abs(dVth) * 1e3)) ** 2)
        R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        ax.plot([], [], ' ', label=rf'Linear fit $R^2$={R2:.3f}')

    ax.set_xlabel(r'Cumulative $Q_G$ (C)')
    ax.set_ylabel(r'$|\Delta V_{\mathrm{th}}|$ (mV)')
    ax.set_title(r'Gate Charge vs Hysteresis Window')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))

    fig.tight_layout()

    if export_data and output_SMU is not None:
        from pathlib import Path
        out = Path(output_SMU)
        out.mkdir(exist_ok=True)
        fig.savefig(str(out / f'charge_correlation.{fig_format}'), dpi=600,
                    bbox_inches='tight', transparent=plot_transparency,
                    format=fig_format)

    # Titles off unless one is asked for - a figure caption carries
    # the description in a paper.
    ax.set_title(plot_title if plot_title else '')

    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def update_plot_string(smu_objects: list) -> list:
    '''Update the plot_string variable for each measurement in the SMU objects.'''
    # Loop over the SMU objects in the list
    for smu_obj in smu_objects: 
        # Update CV plot strings
        for i in range(len(smu_obj.cv_plot_strings)):
            current_string = smu_obj.cv_plot_strings[i]
            new_string = input(f"Current CV label: '{current_string}'. Enter new label (or press Enter to keep): ")
            if new_string.strip():
                smu_obj.cv_plot_strings[i] = new_string.strip()
        
        # Update IV plot strings  
        for i in range(len(smu_obj.iv_plot_strings)):
            current_string = smu_obj.iv_plot_strings[i]
            new_string = input(f"Current IV label: '{current_string}'. Enter new label (or press Enter to keep): ")
            if new_string.strip():
                smu_obj.iv_plot_strings[i] = new_string.strip()

        # Update IV3T plot strings
        for i in range(len(smu_obj.iv3T_plot_strings)):
            current_string = smu_obj.iv3T_plot_strings[i]
            new_string = input(f"Current IV3T label: '{current_string}'. Enter new label (or press Enter to keep): ")
            if new_string.strip():
                smu_obj.iv3T_plot_strings[i] = new_string.strip()
                
    return smu_objects
