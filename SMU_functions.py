#import all the libraries needed
from import_dep import *
from scipy.signal import medfilt

def plot_IV(smu_data: list, 
           run_nums: Optional[list] = None,
           x_lim: Optional[Tuple[float, float]] = None,
           y_lim: Optional[Tuple[float, float]] = None,
           log_plot: bool = False,
           plot_arrows: bool = False,
           line_width: Optional[float] = None,
           plot_key: bool = True,
           export_data: bool = False,
           output_SMU: Optional[str] = None,
           output_path: Optional[str] = None,
           fig_format: str = 'png',
           plot_transparency: bool = False) -> Figure:
    """
    Plot I-V (current-voltage) characteristics from SMU data.
    
    Parameters:
        smu_data          : List of SMU objects containing IV data
        run_nums          : Optional list of run numbers to plot. If None, plots all.
        x_lim             : Optional tuple specifying the x-axis limits (min, max).
        y_lim             : Optional tuple specifying the y-axis limits (min, max).
        log_plot          : Whether to use logarithmic scale for y-axis.
        plot_arrows       : Whether to add directional arrows showing voltage sweep direction.
        line_width        : Optional line width for the plot. If None, uses default.
        plot_key          : Whether to show the legend/key. Default is True.
        export_data       : Whether to use export styling and generate output path.
        output_SMU        : Directory path for saving figures when export_data=True.
        output_path       : Specific path to save the figure (overrides auto-generation).
        fig_format        : Format for saved figure ('png', 'svg', etc.).
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
        iv_data, plot_strings = smu_obj.get_iv_data(run_nums)
        
        for df, label_str in zip(iv_data, plot_strings):
            if not df.empty:
                x = df['V_V']
                y = df['I_A']
                
                # Apply x_lim masking if specified to rescale y-axis appropriately
                if x_lim is not None:
                    mask = (x >= x_lim[0]) & (x <= x_lim[1])
                    x = x[mask]
                    y = y[mask]
                    
                    # Skip if no data points in the specified range
                    if len(x) == 0:
                        continue
                
                # For log plots, use absolute value to avoid log(negative) issues
                if log_plot:
                    y = np.abs(y)
                
                # Plot the main line with optional line width
                plot_kwargs = {'label': label_str, 'marker': 'o'}
                if line_width is not None:
                    plot_kwargs['linewidth'] = line_width
                
                line = ax.plot(x, y, **plot_kwargs)
                
                # Add directional arrows if requested
                if plot_arrows and len(x) > 10:  # Only add arrows if we have enough points
                    color = line[0].get_color()  # Get the color of the line
                    
                    # Apply 5-point median filter to smooth data for arrow calculation only
                    # Convert to numpy arrays for easier processing
                    x_smooth = medfilt(x.values, kernel_size=3)
                    y_smooth = medfilt(y.values, kernel_size=3)
                    
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
                                x_range = np.ptp(x)
                                y_range = np.ptp(y)
                                arrow_length = min(x_range, y_range) * 0.05  # 5% of data range
                                
                                # Position the arrow using original (unsmoothed) data
                                arrow_x = x.iloc[idx]
                                arrow_y = y.iloc[idx]
                                
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
    ax.set_ylabel(r"Current (A)")
    ax.set_title(r"I-V Characteristics")
    
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
    
    # Auto-generate output path if export_data is True but no specific path provided
    if export_data and output_path is None and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        output_path = str(output_dir / f"SMU_IV_plot.{fig_format}")
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', 
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
           plot_key: bool = True,
           export_data: bool = False,
           output_SMU: Optional[str] = None,
           output_path: Optional[str] = None,
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
        plot_key          : Whether to show the legend/key. Default is True.
        export_data       : Whether to use export styling and generate output path.
        output_SMU        : Directory path for saving figures when export_data=True.
        output_path       : Specific path to save the figure (overrides auto-generation).
        fig_format        : Format for saved figure ('png', 'svg', etc.).
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
                    ax.plot(x, y_cp, **plot_kwargs_cp)
                
                if plot_Gp:
                    y_gp = df['Gp'][mask] if x_lim is not None else df['Gp']
                    plot_kwargs_gp = {'label': rf"{label_str} ($G_p$)", 'marker': 's', 'linestyle': '--'}
                    if line_width is not None:
                        plot_kwargs_gp['linewidth'] = line_width
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
    
    # Auto-generate output path if export_data is True but no specific path provided
    if export_data and output_path is None and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        output_path = str(output_dir / f"SMU_CV_plot.{fig_format}")
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', 
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
                   plot_key: bool = True,
                   export_data: bool = False,
                   output_SMU: Optional[str] = None,
                   output_path: Optional[str] = None,
                   fig_format: str = 'png',
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
        plot_key          : Whether to show the legend/key. Default is True.
        export_data       : Whether to use export styling and generate output path.
        output_SMU        : Directory path for saving figures when export_data=True.
        output_path       : Specific path to save the figure (overrides auto-generation).
        fig_format        : Format for saved figure ('png', 'svg', etc.).
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
        iv_data, plot_strings = smu_obj.get_iv_data(run_nums)
        
        for df, label_str in zip(iv_data, plot_strings):
            if not df.empty:
                x = df['V_V']
                y = df['R_Ohm']
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
                
                # Plot the main line with optional line width
                plot_kwargs = {'label': label_str, 'marker': 'o'}
                if line_width is not None:
                    plot_kwargs['linewidth'] = line_width
                
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
    ax.set_ylabel(r"Resistance ($\Omega$)")
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
    
    # Auto-generate output path if export_data is True but no specific path provided
    if export_data and output_path is None and output_SMU is not None:
        from pathlib import Path
        output_dir = Path(output_SMU)
        output_dir.mkdir(exist_ok=True)
        output_path = str(output_dir / f"SMU_resistance_plot.{fig_format}")
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', 
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
