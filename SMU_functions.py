#import all the libraries needed
from import_dep import *
from scipy.signal import medfilt

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
    
    if max_abs_current >= 1e-3:  # >= 1 mA
        return current_data, "A", 1
    elif max_abs_current >= 1e-6:  # >= 1 μA  
        return current_data * 1e3, "mA", 1e3
    elif max_abs_current >= 1e-9:  # >= 1 nA
        return current_data * 1e6, "μA", 1e6
    else:  # < 1 nA
        return current_data * 1e9, "nA", 1e9

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
           line_width: Optional[float] = None,
           markersize: Optional[float] = None,
           plot_key: bool = True,
           export_data: bool = False,
           output_SMU: Optional[str] = None,
           fig_name: Optional[str] = None,
           fig_format: str = 'png',
           plot_transparency: bool = False) -> Figure:
    """
    Plot I-V characteristics with different conduction mechanism analysis modes.
    
    Supported plot modes for identifying conduction mechanisms in thin-film /
    memristive devices:
    
    'linear' (default):
        Standard I vs V plot.
        I(V) = V / R  (Ohmic conduction)
        Useful for basic characterisation and identifying linear/ohmic regions.
    
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
                            'linear', 'SCLC', 'schottky', 'PF', 'FN', 'TAT'.
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
    
    valid_modes = {'linear', 'SCLC', 'schottky', 'PF', 'FN', 'TAT'}
    if plot_mode not in valid_modes:
        raise ValueError(f"plot_mode must be one of {valid_modes}, got '{plot_mode}'")
    
    # Import here to avoid circular imports
    from plot_style import set_plot_style
    
    # Set plot style based on export_data
    fig_size = set_plot_style(export_data=export_data, use_tex=True)
    
    fig, ax = plt.subplots(figsize=fig_size)
    
    # For 'linear' mode, auto-scale current units; analysis modes use raw Amperes
    current_unit = "A"
    current_scale = 1
    if plot_mode == 'linear':
        all_current_data = []
        for smu_obj in smu_data:
            iv_data, plot_strings = smu_obj.get_iv_data(run_nums)
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
    for smu_obj in smu_data:
        iv_data, plot_strings = smu_obj.get_iv_data(run_nums)
        
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
                plot_kwargs = {'label': label_str, 'marker': 'o'}
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
        fig.savefig(output_file, dpi=300, bbox_inches='tight', 
                   transparent=plot_transparency, format=fig_format)
    
    plt.show()
    return fig


def plot_CV(smu_data: list,
           run_nums: Optional[list] = None, 
           x_lim: Optional[Tuple[float, float]] = None,
           y_lim: Optional[Tuple[float, float]] = None,
           log_plot: bool = False,
           plot_Cp: bool = True,
           plot_Gp: bool = False,
           line_width: Optional[float] = None,
           markersize: Optional[float] = None,
           plot_key: bool = True,
           export_data: bool = False,
           output_SMU: Optional[str] = None,
           fig_name: Optional[str] = None,
           fig_format: str = 'png',
           plot_transparency: bool = False) -> Figure:
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
    from plot_style import set_plot_style
    
    # Set plot style based on export_data
    fig_size = set_plot_style(export_data=export_data, use_tex=True)
    
    fig, ax = plt.subplots(figsize=fig_size)
    
    for smu_obj in smu_data:
        cv_data, plot_strings = smu_obj.get_cv_data(run_nums)
        
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
                    plot_kwargs_cp = {'label': rf"{label_str} ($C_p$)", 'marker': 'o'}
                    if line_width is not None:
                        plot_kwargs_cp['linewidth'] = line_width
                    if markersize is not None:
                        plot_kwargs_cp['markersize'] = markersize
                    ax.plot(x, y_cp, **plot_kwargs_cp)
                
                if plot_Gp:
                    y_gp = df['Gp'][mask] if x_lim is not None else df['Gp']
                    plot_kwargs_gp = {'label': rf"{label_str} ($G_p$)", 'marker': 's', 'linestyle': '--'}
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
        fig.savefig(output_file, dpi=300, bbox_inches='tight', 
                   transparent=plot_transparency, format=fig_format)
    
    plt.show()
    return fig


def plot_resistance(smu_data: list,
                   run_nums: Optional[list] = None,
                   x_lim: Optional[Tuple[float, float]] = None,
                   y_lim: Optional[Tuple[float, float]] = None,
                   log_plot: bool = True,
                   plot_arrows: bool = False,
                   line_width: Optional[float] = None,
                   markersize: Optional[float] = None,
                   plot_key: bool = True,
                   export_data: bool = False,
                   output_SMU: Optional[str] = None,
                   fig_name: Optional[str] = None,
                   fig_format: str = 'tiff',
                   plot_transparency: bool = False) -> Figure:
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
    from plot_style import set_plot_style
    
    # Set plot style based on export_data
    fig_size = set_plot_style(export_data=export_data, use_tex=True)
    
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Collect all resistance data to determine appropriate scaling
    all_resistance_data = []
    for smu_obj in smu_data:
        iv_data, plot_strings = smu_obj.get_iv_data(run_nums)
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
    for smu_obj in smu_data:
        iv_data, plot_strings = smu_obj.get_iv_data(run_nums)
        
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
                plot_kwargs = {'label': label_str, 'marker': 'o'}
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
        fig.savefig(output_file, dpi=300, bbox_inches='tight', 
                   transparent=plot_transparency, format=fig_format)
    
    plt.show()
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
                
    return smu_objects
