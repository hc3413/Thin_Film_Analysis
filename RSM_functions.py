#import all the libraries needed
from import_dep import *
from plot_style import set_plot_style
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

def _plot_on_ax(
    ax: plt.Axes,
    d: pd.DataFrame,
    class_obj: Any,
    count_d: int,
    threshold_filt: float = 0,
    grid_filt: float = 100,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal',
    colormap: str = 'plasma',
    label_str: Optional[str] = None,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    contour_color: bool = False,
    show_colorbar: bool = True,
    show_key: bool = True,
    first_index: bool = False,  # New parameter to prevent double plotting of labels
    resolution_scaling: Optional[float] = None, # New parameter to scale the resolution of the meshgrid
    intensity_lim: Optional[list[float]] = None, # New parameter to set intensity limits
    show_title: bool = False, # New parameter to show/hide plot title
    fitted_peaks: Optional[pd.DataFrame] = None, # Overlay fitted peak markers
) -> None:
    """
    A helper function to plot a single dataset on the provided axis.

    Parameters:
        ax           : The matplotlib Axes instance to draw on.
        d            : Dictionary containing the data (e.g., 'qx', 'qz', 'Intensity').
        class_obj    : An object that may hold additional data (like lattice parameters).
        count_d      : Index to reference specific dataset properties from class_obj.
        threshold_filt: Intensity threshold for filtering data.
        grid_filt:    Filter the grid data based on threshold.
        plot_type    : Type of plot to generate: 'scatter', 'contour', or 'mesh'.
        ax_unit      : Type of axis unit to use, e.g., 'reciprocal' or 'lattice'.
        colormap     : The colormap to use for plotting.
        label_str    : Label for the plot legend.
        x_lim        : Optional x-axis limits.
        y_lim        : Optional y-axis limits.
        intensity_lim: Optional intensity limits [min, max].
    """
    # Extract data based on axis unit
    if ax_unit == 'reciprocal':
        # df format
        x = class_obj.qxqz_df[count_d]['qx']
        y = class_obj.qxqz_df[count_d]['qz']
        z = class_obj.qxqz_df[count_d]['Intensity']
        # Numpy format where it is reshaped into meshgrid
        grid_x = class_obj.qxqz_np[count_d][:,:,0]
        grid_y = class_obj.qxqz_np[count_d][:,:,1]
        grid_intensity = class_obj.qxqz_np[count_d][:,:,2]
        
        ax.set_xlabel(r"$Q_{\parallel} \,(2\pi/\mathrm{\AA})$")
        ax.set_ylabel(r"$Q_{\perp} \,(2\pi/\mathrm{\AA})$")
        if show_title:
            ax.set_title(r"Reciprocal Space Map ($Q_x$ vs $Q_z$)")
    elif ax_unit == 'lattice':
        # df format
        x = class_obj.lat_param_df[count_d]['a']
        y = class_obj.lat_param_df[count_d]['c']
        z = class_obj.lat_param_df[count_d]['Intensity']
        # Numpy format where it is reshaped into meshgrid
        grid_x = class_obj.lat_param_np[count_d][:,:,0]
        grid_y = class_obj.lat_param_np[count_d][:,:,1]
        grid_intensity = class_obj.lat_param_np[count_d][:,:,2]
        
        ax.set_xlabel(r"$a \,(\mathrm{\AA})$")
        ax.set_ylabel(r"$c \,(\mathrm{\AA})$")
        if show_title:
            ax.set_title(r"Lattice Parameters ($a$ vs $c$)")
    else:
        raise ValueError("Unknown ax_unit: choose 'reciprocal' or 'lattice'")
    
    
     # Interpolate data if resolution_scaling is specified
    if resolution_scaling is not None:
        # Scale and interpolate the meshgrid data
        nx = int(grid_x.shape[1] * resolution_scaling)
        ny = int(grid_y.shape[0] * resolution_scaling)
        xi = np.linspace(grid_x.min(), grid_x.max(), nx)
        yi = np.linspace(grid_y.min(), grid_y.max(), ny)
        xi, yi = np.meshgrid(xi, yi)
        grid_intensity = griddata((grid_x.flatten(), grid_y.flatten()), grid_intensity.flatten(), (xi, yi), method='nearest')
        grid_x, grid_y = xi, yi
        
        # Scale and interpolate the scatter data
        x = xi.flatten()
        y = yi.flatten()
        z = grid_intensity.flatten()
        
    

    # Filter data based on threshold
    filt_indices = z > threshold_filt
    x_filt = x[filt_indices]
    y_filt = y[filt_indices]
    plot_intensity = z[filt_indices]

        
    # Color contour lines based on colormap or default to black
    c_map = cm.get_cmap(colormap)
    middle_color = c_map(0.5)  # This returns an RGBA tuple 
    if contour_color == True:
        c_color =[middle_color] 
    else:
        c_color = 'black'

    # Filtering for the contours only
    grid_intensity[grid_intensity < grid_filt] = 0
    
    # Apply logarithmic scaling to plot_intensity for alpha mapping
    log_plot_intensity = np.log10(plot_intensity + 1e-99)

    # Normalize log_plot_intensity to range [0, 1] for alpha mapping
    norm_log_intensity = (log_plot_intensity - log_plot_intensity.min()) / (log_plot_intensity.max() - log_plot_intensity.min())    
    
    norm_log_intensity[norm_log_intensity > 0.1] = 1
    norm_log_intensity[norm_log_intensity < 0.1] += 0.06
    
    # Invert the normalized intensity to get higher transparency for smaller values
    alpha_values = norm_log_intensity
    
    vmin = None
    vmax = None
    if intensity_lim is not None:
        vmin = 10**intensity_lim[0]
        vmax = 10**intensity_lim[1]

    # Generate plot based on plot_type
    if plot_type == 'scatter':
        p_out = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=colormap, s=2.5, edgecolors='none', label=label_str, norm=LogNorm(vmin=vmin, vmax=vmax),alpha=alpha_values)
    
    elif plot_type == 'contour':
        
        p_out = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=colormap, s=2.5, edgecolors='none',  norm=LogNorm(vmin=vmin, vmax=vmax),alpha=alpha_values)
        
        if intensity_lim is not None:
             contour_levels = np.logspace(intensity_lim[0], intensity_lim[1], 10)
        else:
             contour_levels = np.logspace(np.log10(1 + 1e-99), np.log10(plot_intensity.max()), 10)
        
        # generate contours from the meshgrid data to prevent errors/streaking
        ax.contour(grid_x, grid_y, grid_intensity, levels=contour_levels, colors=c_color, linewidths=0.3)
        
        scatter_for_legend = ax.scatter([], [], c=c_color, label=label_str)
        
    elif plot_type == 'mesh':

        
        p_out = ax.pcolormesh(grid_x, grid_y, grid_intensity, cmap=colormap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax))
        scatter_for_legend = ax.scatter([], [], c=c_color, label=label_str)
        
        # Apply Gaussian filter to the interpolated data for contour plot
        #smoothed_intensity = gaussian_filter(grid_intensity, sigma=0.00001)
        
        if intensity_lim is not None:
             contour_levels = np.logspace(intensity_lim[0], intensity_lim[1], 10)
        else:
             contour_levels = np.logspace(np.log10(1 + 1e-99), np.log10(plot_intensity.max()), 10) 
        ax.contour(grid_x, grid_y, grid_intensity, levels=contour_levels, colors=c_color, linewidths=0.2)
        
    elif plot_type == 'surf':
        # Set up the 3D plot
        
        p_out = ax.plot_surface(grid_x, grid_y, grid_intensity, cmap=colormap, linewidth=0, antialiased=False, norm=LogNorm(vmin=vmin, vmax=vmax))
        
        
    elif plot_type == 'triangle':
        
        ax.plot_trisurf(x_filt, y_filt, plot_intensity, cmap=colormap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax))
        # Add a color bar
        p_out = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=colormap, norm=LogNorm(vmin=vmin, vmax=vmax))
        
    
    else:
        raise ValueError("Unknown plot_type: choose 'scatter', 'contour', or 'mesh'")        

    # Generate a colorbar for the plot
    if show_colorbar: 
            cbar = plt.colorbar(p_out, ax=ax, pad=0.02, format=LogFormatterSciNotation())  
            cbar.set_label("Intensity")
            
    # Add guidelines for lattice plots
    if ax_unit == 'lattice' and first_index == True:
        ax.scatter(4.12, 4.12, color='purple', marker='+', s=100, label='BaSnO3 bulk')
        ax.scatter(4.035, 4.033, color='orange', marker='+', s=100, label='SrSnO3 bulk')
        ax.scatter(4.05, 4.05, color='yellow', marker='+', s=100, label='LaScO3 bulk')
        ax.axvline(x=3.905, color='green', linestyle='--', linewidth=1)
        ax.axhline(y=3.905, color='green', linestyle='--', linewidth=1)
        
    # Add guidelines for lattice plots
    if ax_unit == 'reciprocal' and first_index == True:
        #ax.scatter(1.526, 4.571, color='red', marker='+', s=100, label='BaSnO3 on STO')
        #ax.scatter(1.581, 4.555, color='green', marker='+', s=100, label='SrSnO3 on STO')
        #ax.scatter(2*np.pi/4.05, 6*np.pi/4.05, color='yellow', marker='+', s=100, label='LaScO3 bulk')
        ax.scatter(2*np.pi/4.12, 6*np.pi/4.12, color='purple', marker='+', s=100, label='BaSnO3 bulk')
        ax.scatter(2*np.pi/4.036, 6*np.pi/4.036, color='orange', marker='+', s=100, label='SrSnO3 bulk')
        #ax.axvline(x=1.61, color='green', linestyle='--', linewidth=1)
        #ax.axhline(y=1.61*3, color='green', linestyle='--', linewidth=1)

    # Set axis limits if provided
    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    # Overlay fitted peak markers if provided
    if fitted_peaks is not None and not fitted_peaks.empty:
        for _, pk in fitted_peaks.iterrows():
            if ax_unit == 'reciprocal':
                px, pz = pk['qx'], pk['qz']
            elif ax_unit == 'lattice':
                px, pz = pk['a'], pk['c']
            else:
                continue
            marker_color = 'lime' if pk.get('is_substrate', False) else 'white'
            ax.scatter(px, pz, marker='+', s=120, linewidths=1.5,
                       color=marker_color, edgecolors='black', zorder=10)
            label_text = pk.get('material', None)
            if label_text:
                ax.annotate(label_text, (px, pz), textcoords='offset points',
                            xytext=(6, 6), fontsize=7, color=marker_color,
                            path_effects=[pe.withStroke(linewidth=1.5, foreground='black')],
                            zorder=11)

    # Add a legend if a label is provided
    #if label_str is not None:
    if show_key == True:
        ax.legend(loc='upper left', framealpha=0.4)


def create_custom_colormap(base_cmap='jet'):
    'Need to add white to the base of the map'
    # Get the colormap
    cmap = plt.get_cmap(base_cmap)
    
    # Get the 'Greys' colormap
    greys_cmap = plt.get_cmap('Greys')
    
    # Convert the colormaps to lists of colors
    cmap_colors = cmap(np.arange(cmap.N))
    greys_colors = greys_cmap(np.linspace(0, 0.5, cmap.N // 5))
    
    # Replace the first 1/8 of the base colormap with the greyscale transition
    cmap_colors[:cmap.N // 5] = greys_colors
    
    # Create a new colormap from the modified list of colors
    custom_cmap = ListedColormap(cmap_colors)
    
    return custom_cmap

# Create the custom colormap
#custom_cmap = create_custom_colormap()
custom_cmap = 'jet'


def RSM_plot_sep(
    dat: Tuple[Any],
    threshold_filt: float = 0,
    grid_filt: float = 100,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal',
    show_key: bool = True,
    resolution_scaling: Optional[float] = None, # New parameter to scale the resolution of the meshgrid
    string_suppress: bool = False,  # New parameter to suppress string labels
    intensity_lim: Optional[List[float]] = None, # New parameter to set intensity limits
    title_on: bool = False, # New parameter to show/hide plot title
    color_bar: bool = True, # New parameter to show/hide colorbar
) -> Figure:
    """
    Plot multiple datasets in separate subplots.

    Parameters:
        dat: A sequence of objects containing data. Each object is expected to have
             a 'qxqz_df' attribute (a sequence of dictionaries with data) and may have
             a 'plot_string' attribute and 'lat_param_df' for lattice parameters.
        threshold_filt: Intensity threshold used to filter data.
        grid_filt: Filter the grid data based on threshold.
        x_lim: Optional tuple specifying the x-axis limits (min, max).
        y_lim: Optional tuple specifying the y-axis limits (min, max).
        plot_type: Type of plot to generate; one of 'scatter', 'contour', or 'mesh'.
        ax_unit: Unit type for the axes; either 'reciprocal' or 'lattice'.
        resolution_scaling: Optional parameter to scale the resolution of the meshgrid.
        intensity_lim: Optional intensity limits [min, max].

    Returns:
        The matplotlib Figure object containing the subplots.
    """
    total_plots = sum(len(class_obj.qxqz_df) for class_obj in dat)
    n_rows = total_plots // 2 + total_plots % 2

    fig = plt.figure(figsize=(15, 5.5 * n_rows))
    gs = fig.add_gridspec(n_rows, 2)
    
    subplot_counter = 0
   
    for class_obj in dat:
        first_index = True # New parameter to prevent double plotting of labels
        for count_d, d in enumerate(class_obj.qxqz_df):
            ax = fig.add_subplot(gs[subplot_counter // 2, subplot_counter % 2])
            
            # Only add the label string if string_suppress is False
            if string_suppress:
                label_str = None
            else:
                label_str = (class_obj.plot_string[count_d] if hasattr(class_obj, 'plot_string') else None)
                
            _plot_on_ax(ax, d, class_obj, count_d, threshold_filt = threshold_filt,
                        grid_filt = grid_filt,
                        plot_type = plot_type, ax_unit = ax_unit, colormap= custom_cmap , label_str=label_str,
                        x_lim=x_lim, y_lim=y_lim, resolution_scaling=resolution_scaling, first_index=first_index,
                        show_key=show_key, intensity_lim=intensity_lim, show_title=title_on, show_colorbar=color_bar)
            subplot_counter += 1
            first_index = False
            
    fig.tight_layout()
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    plt.show()
    return fig


def RSM_plot_single(
    dat: Any,
    threshold_filt: float = 0,
    grid_filt: float = 100,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal',
    show_key: bool = True,
    resolution_scaling: Optional[float] = None,
    string_suppress: bool = False,
    intensity_lim: Optional[List[float]] = None,
    title_on: bool = False,
    color_bar: bool = True,
    export_data: bool = False,
    fitted_peaks: Optional[pd.DataFrame] = None,
) -> Figure:
    """
    Plot a single RSM dataset with consistent figure styling.

    Parameters:
        dat: A single data object containing RSM data. Expected to have
             a 'qxqz_df' attribute (a sequence of dictionaries with data) and may have
             a 'plot_string' attribute and 'lat_param_df' for lattice parameters.
        threshold_filt: Intensity threshold used to filter data.
        grid_filt: Filter the grid data based on threshold.
        x_lim: Optional tuple specifying the x-axis limits (min, max).
        y_lim: Optional tuple specifying the y-axis limits (min, max).
        plot_type: Type of plot to generate; one of 'scatter', 'contour', or 'mesh'.
        ax_unit: Unit type for the axes; either 'reciprocal' or 'lattice'.
        resolution_scaling: Optional parameter to scale the resolution of the meshgrid.
        intensity_lim: Optional intensity limits [min, max].
        title_on: Whether to show the plot title (default: False).
        color_bar: Whether to show the colorbar (default: True).
        export_data: If True, use publication-ready figure size [2.625, 3.5].
                     If False, use visualisation size [6, 9].

    Returns:
        The matplotlib Figure object containing the plot.
    """
    # Get figure size from plot_style and invert x/y for RSM plots
    fig_size = set_plot_style(export_data=export_data)
    fig_size_inverted = [fig_size[1], fig_size[0]]  # Invert x and y
    
    fig, ax = plt.subplots(figsize=fig_size_inverted)
    
    # Handle both single object and list input
    if isinstance(dat, list):
        class_obj = dat[0]  # Take first element if list is passed
    else:
        class_obj = dat
    
    count_d = 0
    d = class_obj.qxqz_df[count_d]
    
    # Only add the label string if string_suppress is False
    if string_suppress:
        label_str = None
    else:
        label_str = (class_obj.plot_string[count_d] if hasattr(class_obj, 'plot_string') else None)
    
    _plot_on_ax(ax, d, class_obj, count_d, threshold_filt=threshold_filt,
                grid_filt=grid_filt,
                plot_type=plot_type, ax_unit=ax_unit, colormap=custom_cmap, label_str=label_str,
                x_lim=x_lim, y_lim=y_lim, resolution_scaling=resolution_scaling, first_index=True,
                show_key=show_key, intensity_lim=intensity_lim, show_title=title_on, show_colorbar=color_bar,
                fitted_peaks=fitted_peaks)
    
    fig.tight_layout()
    plt.show()
    return fig


def RSM_plot_fitting(
    dat: Any,
    fitted_peaks: pd.DataFrame,
    threshold_filt: float = 10,
    grid_filt: float = 90,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    plot_type: str = 'mesh',
    ax_unit: str = 'reciprocal',
    resolution_scaling: Optional[float] = None,
    intensity_lim: Optional[List[float]] = None,
    export_data: bool = False,
    residuals: bool = False,
) -> Figure:
    """
    Display RSM data alongside fitted simulation for visual validation.

    Left panel:  Measured RSM data with peak markers and STO exclusion zone.
    Right panel: Simulated intensity from the fitted Pseudo-Voigt peaks.

    Uses the same threshold_filt as peak_search's intensity_threshold so the
    visible data matches what the fitting algorithm operated on.

    Parameters:
        dat            : RSM data object (or single-element list).
        fitted_peaks   : DataFrame returned by peak_search().
        threshold_filt : Intensity threshold (should match intensity_threshold in peak_search).
        grid_filt      : Contour grid filter.
        x_lim          : X-axis limits.
        y_lim          : Y-axis limits.
        plot_type      : Plot type ('mesh', 'scatter', 'contour').
        ax_unit        : Axis unit ('reciprocal' or 'lattice').
        resolution_scaling : Grid resolution multiplier.
        intensity_lim  : Log10 intensity limits [min, max].
        export_data    : If True, use publication figure size.
        residuals      : If True, add a third panel showing measured minus simulated
                         to reveal any unmodelled peaks.

    Returns:
        The matplotlib Figure object.
    """
    fig_size = set_plot_style(export_data=export_data)
    fig_size_inverted = [fig_size[1], fig_size[0]]
    # Panel layout: 2 columns (+ optional residual row)
    n_cols = 2
    n_rows = 2 if residuals else 1
    fig_width = fig_size_inverted[0] * 2.1
    fig_height = fig_size_inverted[1] * (1.8 if residuals else 1.0)

    if residuals:
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height),
                                 gridspec_kw={'height_ratios': [1, 0.8]})
        ax_data, ax_sim = axes[0, 0], axes[0, 1]
        ax_resid = axes[1, 0]
        axes[1, 1].axis('off')  # empty quadrant
    else:
        fig, (ax_data, ax_sim) = plt.subplots(1, 2, figsize=(fig_width, fig_size_inverted[1]))

    if isinstance(dat, list):
        class_obj = dat[0]
    else:
        class_obj = dat

    count_d = 0
    d = class_obj.qxqz_df[count_d]

    # --- Left panel: measured RSM data ---
    _plot_on_ax(ax_data, d, class_obj, count_d, threshold_filt=threshold_filt,
                grid_filt=grid_filt,
                plot_type=plot_type, ax_unit=ax_unit, colormap=custom_cmap,
                label_str=None,
                x_lim=x_lim, y_lim=y_lim, resolution_scaling=resolution_scaling,
                first_index=False, show_key=False, intensity_lim=intensity_lim,
                show_title=False, show_colorbar=False,
                fitted_peaks=None)
    ax_data.set_title('Measured', fontsize=9)

    # --- Overlay fitted peak markers on left panel ---
    _overlay_peak_markers(ax_data, fitted_peaks, ax_unit)

    # --- Right panel: simulated intensity from fitted peaks ---
    grid_data = class_obj.qxqz_np[count_d]
    grid_qx = grid_data[:, :, 0]
    grid_qz = grid_data[:, :, 1]

    # Optionally interpolate to higher resolution (match left panel)
    if resolution_scaling is not None:
        nx = int(grid_qx.shape[1] * resolution_scaling)
        ny = int(grid_qx.shape[0] * resolution_scaling)
        xi = np.linspace(grid_qx.min(), grid_qx.max(), nx)
        yi = np.linspace(grid_qz.min(), grid_qz.max(), ny)
        sim_qx, sim_qz = np.meshgrid(xi, yi)
    else:
        sim_qx = grid_qx
        sim_qz = grid_qz

    # Build simulated intensity from all fitted peaks
    sim_int = np.zeros_like(sim_qx)
    coords_flat = np.array([sim_qx.flatten(), sim_qz.flatten()])

    if fitted_peaks is not None and not fitted_peaks.empty:
        # Use the film background (from any film peak row, they share the same value)
        film_rows = fitted_peaks[~fitted_peaks['is_substrate']]
        if not film_rows.empty and 'background' in film_rows.columns:
            bg = film_rows.iloc[0]['background']
        else:
            bg = 0.0
        sim_flat = np.full(coords_flat.shape[1], bg)

        for _, pk in fitted_peaks.iterrows():
            sim_flat += _pseudo_voigt_2d(
                coords_flat, pk['intensity'], pk['qx'], pk['qz'],
                pk['sigma_x'], pk['sigma_z'], pk['eta'],
            )
        sim_int = sim_flat.reshape(sim_qx.shape)

    # Apply same threshold as the data panel
    sim_int[sim_int < threshold_filt] = 0.0

    vmin = None
    vmax = None
    if intensity_lim is not None:
        vmin = 10**intensity_lim[0]
        vmax = 10**intensity_lim[1]

    # Plot simulated RSM
    sim_int_plot = np.ma.masked_where(sim_int <= 0, sim_int)
    p_sim = ax_sim.pcolormesh(sim_qx, sim_qz, sim_int_plot, cmap=custom_cmap,
                              shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax))
    ax_sim.set_title('Fitted simulation', fontsize=9)
    ax_sim.set_xlabel(r"$Q_{\parallel} \,(2\pi/\mathrm{\AA})$")
    ax_sim.set_ylabel(r"$Q_{\perp} \,(2\pi/\mathrm{\AA})$")

    if x_lim is not None:
        ax_sim.set_xlim(x_lim)
    if y_lim is not None:
        ax_sim.set_ylim(y_lim)

    # Overlay peak markers on simulation panel too
    _overlay_peak_markers(ax_sim, fitted_peaks, ax_unit)

    # --- Residual panel (measured - simulated) ---
    if residuals:
        # Interpolate the measured data onto the simulation grid
        meas_grid_data = class_obj.qxqz_np[count_d]
        meas_qx = meas_grid_data[:, :, 0]
        meas_qz = meas_grid_data[:, :, 1]
        meas_int = meas_grid_data[:, :, 2].astype(float).copy()

        # Threshold the measured data the same way
        meas_int[meas_int < threshold_filt] = 0.0

        # Interpolate measured data onto the simulation grid
        from scipy.interpolate import griddata as _griddata
        meas_on_sim = _griddata(
            (meas_qx.flatten(), meas_qz.flatten()),
            meas_int.flatten(),
            (sim_qx, sim_qz),
            method='nearest',
        )

        # Compute residual: measured - simulated (both on same grid, same threshold)
        resid = meas_on_sim - sim_int

        # Zero out the STO exclusion zone so it doesn't dominate the residual
        if fitted_peaks is not None and not fitted_peaks.empty:
            sto_rows = fitted_peaks[fitted_peaks['is_substrate'] == True]
            for _, pk in sto_rows.iterrows():
                er = pk.get('sto_exclusion_radius', 0.1)
                dist_sto = np.sqrt((sim_qx - pk['qx'])**2 + (sim_qz - pk['qz'])**2)
                resid[dist_sto < er] = 0.0

        # Use a diverging colormap so positive = unmodelled signal
        resid_abs_max = max(np.abs(resid).max(), 1e-6)
        from matplotlib.colors import TwoSlopeNorm
        resid_norm = TwoSlopeNorm(vmin=-resid_abs_max, vcenter=0, vmax=resid_abs_max)

        p_resid = ax_resid.pcolormesh(sim_qx, sim_qz, resid, cmap='RdBu_r',
                                      shading='auto', norm=resid_norm)
        plt.colorbar(p_resid, ax=ax_resid, pad=0.02, label='Measured \u2212 Simulated')
        ax_resid.set_title('Residual (Measured \u2212 Simulated)', fontsize=9)
        ax_resid.set_xlabel(r"$Q_{\parallel} \,(2\pi/\mathrm{\AA})$")
        ax_resid.set_ylabel(r"$Q_{\perp} \,(2\pi/\mathrm{\AA})$")

        if x_lim is not None:
            ax_resid.set_xlim(x_lim)
        if y_lim is not None:
            ax_resid.set_ylim(y_lim)

        _overlay_peak_markers(ax_resid, fitted_peaks, ax_unit)

    fig.tight_layout()
    plt.show()
    return fig


def _overlay_peak_markers(ax: plt.Axes, fitted_peaks: Optional[pd.DataFrame],
                          ax_unit: str = 'reciprocal') -> None:
    """Draw numbered peak markers and STO exclusion zone on an axis."""
    if fitted_peaks is None or fitted_peaks.empty:
        return

    for _, pk in fitted_peaks.iterrows():
        if ax_unit == 'reciprocal':
            px, pz = pk['qx'], pk['qz']
        elif ax_unit == 'lattice':
            px, pz = pk['a'], pk['c']
        else:
            continue

        is_sub = pk.get('is_substrate', False)

        # Peak position marker — black for all
        ax.scatter(px, pz, marker='+', s=120, linewidths=1.5,
                   color='black', zorder=10)

        # Numbered label (STO gets "STO", film peaks get 1, 2, 3...)
        if is_sub:
            label_text = 'STO'
        else:
            # peak_id 0 is STO when present, so film peaks start from 1
            label_text = str(int(pk['peak_id']))
        ax.annotate(label_text, (px, pz), textcoords='offset points',
                    xytext=(6, 6), fontsize=8, fontweight='bold', color='black',
                    path_effects=[pe.withStroke(linewidth=2, foreground='white')],
                    zorder=11)

        # STO exclusion zone (red dashed)
        if (ax_unit == 'reciprocal' and is_sub
                and 'sto_exclusion_radius' in pk.index):
            er = pk['sto_exclusion_radius']
            excl = plt.Circle((px, pz), er, fill=False,
                              edgecolor='red', linewidth=1.2,
                              linestyle='--', alpha=0.9, zorder=9)
            ax.add_patch(excl)


def RSM_plot_same(
    dat: Tuple[Any],
    threshold_filt: float = 0,
    grid_filt: float = 100,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal',
    show_key: bool = True,
    resolution_scaling: Optional[float] = None, # New parameter to scale the resolution of the meshgrid
    title_on: bool = False, # New parameter to show/hide plot title
    color_bar: bool = True, # New parameter to show/hide colorbar
) -> Figure:
    """
    Plot multiple datasets on a single set of axes.

    Parameters:
        dat: A sequence of objects containing data. Each object is expected to have
             a 'qxqz_df' attribute (a sequence of dictionaries with data) and may have
             a 'plot_string' attribute and 'lat_param_df' for lattice parameters.
        threshold_filt: Intensity threshold used to filter data.
        x_lim: Optional tuple specifying the x-axis limits (min, max).
        y_lim: Optional tuple specifying the y-axis limits (min, max).
        plot_type: Type of plot to generate; one of 'scatter', 'contour', or 'mesh'.
        ax_unit: Unit type for the axes; either 'reciprocal' or 'lattice'.
        resolution_scaling: Optional parameter to scale the resolution of the meshgrid.

    Returns:
        The matplotlib Figure object containing the plot.
    """
    fig, ax = plt.subplots(figsize=(15, 5.5))
    
    col_maps = ['Purples', 'Greens', 'Oranges', 'Blues', 'Reds', 'Greys']
    col_counter = 0
    first_index = True # New parameter to prevent double plotting of labels
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.qxqz_df):
            label_str = (class_obj.plot_string[count_d]
                         if hasattr(class_obj, 'plot_string') else None)
            cmap = col_maps[col_counter % len(col_maps)]
            _plot_on_ax(ax, d, class_obj, count_d, threshold_filt = threshold_filt,
                        grid_filt = grid_filt,
                        plot_type = plot_type, ax_unit = ax_unit, colormap=cmap, label_str=label_str,
                        x_lim=x_lim, y_lim=y_lim, contour_color=True, show_colorbar=color_bar, first_index=first_index, resolution_scaling=resolution_scaling,
                        show_key=show_key, show_title=title_on)
            col_counter += 1
            first_index = False
    fig.tight_layout()
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    plt.show()
    return fig


# =============================================================================
# Materials database for perovskite oxides (103 reflection)
# =============================================================================
PEROVSKITE_MATERIALS = {
    'STO':      {'name': 'SrTiO3',              'a_bulk': 3.905},
    'BSO':      {'name': 'BaSnO3',              'a_bulk': 4.116},
    'SSO':      {'name': 'SrSnO3',              'a_bulk': 4.036},
    'LSO':      {'name': 'LaScO3',              'a_bulk': 4.050},
    'BTO':      {'name': 'BaTiO3',              'a_bulk': 4.010},
    'BSSO_25':  {'name': 'Ba0.25Sr0.75SnO3',    'a_bulk': 4.056},
    'BSSO_40':  {'name': 'Ba0.40Sr0.60SnO3',    'a_bulk': 4.068},
    'BSSO_55':  {'name': 'Ba0.55Sr0.45SnO3',    'a_bulk': 4.080},
    'BSSO_70':  {'name': 'Ba0.70Sr0.30SnO3',    'a_bulk': 4.092},
}

# Material marker shapes for comparison plots
MATERIAL_MARKERS = {
    'STO': '*', 'BSO': 'o', 'SSO': 's', 'LSO': '^', 'BTO': 'd',
    'BSSO': 'X',  # thick X marker for all BSSO variants
}


def _get_material_marker(mat: str, default: str = 'x') -> str:
    """Return marker for a material, matching BSSO prefix for all BSSO_xx variants."""
    if mat in MATERIAL_MARKERS:
        return MATERIAL_MARKERS[mat]
    # Prefix match: BSSO_25, BSSO_40, BSSO_55, BSSO_70 etc. all get the BSSO marker
    for prefix in MATERIAL_MARKERS:
        if mat.startswith(prefix + '_'):
            return MATERIAL_MARKERS[prefix]
    return default


def _material_matches(mat_name: str, select_list: List[str]) -> bool:
    """Check if *mat_name* matches any entry in *select_list* (prefix-aware).

    'BSSO' in the select list will match 'BSSO', 'BSSO_25', 'BSSO_70', etc.
    """
    mat_upper = mat_name.upper()
    for s in select_list:
        s_upper = s.upper()
        if mat_upper == s_upper or mat_upper.startswith(s_upper + '_'):
            return True
    return False

# Known STO (103) position in reciprocal space
STO_QX_REF = 2 * np.pi / 3.905
STO_QZ_REF = 6 * np.pi / 3.905


def _Qx_to_a(Qx):
    return (2 * np.pi) / Qx

def _Qz_to_c(Qz):
    return (6 * np.pi) / Qz


# =============================================================================
# 2D Pseudo-Voigt peak model for XRD peak fitting
# =============================================================================
def _pseudo_voigt_2d(coords: np.ndarray, amplitude: float, x0: float, z0: float,
                     sigma_x: float, sigma_z: float, eta: float) -> np.ndarray:
    """Single 2D Pseudo-Voigt peak: eta*Lorentzian + (1-eta)*Gaussian."""
    x, z = coords
    dx = (x - x0) / sigma_x
    dz = (z - z0) / sigma_z
    r2 = dx**2 + dz**2
    gaussian = np.exp(-0.5 * r2)
    lorentzian = 1.0 / (1.0 + r2)
    eta_c = np.clip(eta, 0, 1)
    return amplitude * (eta_c * lorentzian + (1 - eta_c) * gaussian)


def _multi_peak_residuals(params: np.ndarray, coords: np.ndarray,
                          observed: np.ndarray, n_peaks: int) -> np.ndarray:
    """Residuals for simultaneous multi-peak 2D Pseudo-Voigt fit."""
    background = params[0]
    model = np.full_like(observed, background)
    for i in range(n_peaks):
        p = params[1 + i * 6: 1 + (i + 1) * 6]  # amp, x0, z0, sx, sz, eta
        model += _pseudo_voigt_2d(coords, *p)
    return model - observed


# =============================================================================
# Peak search and fitting
# =============================================================================
def peak_search(
    RSM_data,
    no_peaks: int = 4,
    smooth_sigma: float = 5.0,
    fit_radius: float = 0.1,
    search_index: int = 0,
    sto_search_tol: float = 0.06,
    sto_exclusion_radius: float = 0.1,
    intensity_threshold: float = 10.0,
    min_separation: float = 0.08,
    max_fit_points: int = 2000,
    ftol: float = 1e-6,
    xtol: float = 1e-6,
    max_nfev: int = 1000,
    refine: bool = True,
    qx_range: Optional[Tuple[float, float]] = None,
    qz_range: Optional[Tuple[float, float]] = None,
) -> pd.DataFrame:
    """
    Semi-automatic peak detection and fitting on a 2D reciprocal space map.

    Algorithm:
        1. Gaussian smoothing of the meshgrid data
        2. Auto-detect & individually fit STO substrate peak near known (103) position
        3. Exclude data within sto_exclusion_radius of fitted STO
        4. Apply absolute intensity threshold to remaining data
        5. Local maxima detection on thresholded, STO-excluded data
        6. Select top `no_peaks` film peaks with min_separation enforcement
        7. Simultaneous multi-peak 2D Pseudo-Voigt fit of all film peaks

    Parameters:
        RSM_data              : RSM class instance containing the data.
        no_peaks              : Number of film peaks to search for (excluding STO).
        smooth_sigma          : Sigma for Gaussian smoothing in the coarse search.
        fit_radius            : Radius in reciprocal space units around each peak for fitting.
        search_index          : Index of the dataset within RSM_data to use (default 0).
        sto_search_tol        : Tolerance in Qx/Qz around known STO position for auto-detection.
        sto_exclusion_radius  : Radius around STO to exclude from film peak search.
        intensity_threshold   : Absolute intensity threshold; data below this is excluded from
                                film peak search. Use same value as threshold_filt in plots.
        min_separation        : Minimum distance between distinct film peaks in reciprocal space.
        max_fit_points        : Maximum data points for fitting (subsampled if larger).
        ftol                  : Function tolerance for least_squares solver (lower = more accurate).
        xtol                  : Parameter tolerance for least_squares solver (lower = more accurate).
        max_nfev              : Maximum number of function evaluations (higher = more iterations).
        refine                : If True, do Pseudo-Voigt refinement. If False, use coarse only.
        qx_range              : Optional (min, max) bounds for Qx; peaks outside are ignored.
        qz_range              : Optional (min, max) bounds for Qz; peaks outside are ignored.

    Returns:
        DataFrame with columns: peak_id, qx, qz, a, c, intensity, sigma_x, sigma_z,
                                eta, is_substrate, material, sample_name,
                                fit_radius, sto_exclusion_radius, intensity_threshold
    """
    import time
    t0 = time.time()

    # Extract meshgrid data
    grid_data = RSM_data.qxqz_np[search_index]
    grid_qx = grid_data[:, :, 0]
    grid_qz = grid_data[:, :, 1]
    grid_int = grid_data[:, :, 2]

    sample_name = RSM_data.plot_string[search_index] if hasattr(RSM_data, 'plot_string') else 'unknown'
    sample_code = RSM_data.sample_code if hasattr(RSM_data, 'sample_code') and RSM_data.sample_code else sample_name
    print(f"\n  [peak_search] Sample: {sample_name} (code: {sample_code})")
    print(f"  [peak_search] Grid shape: {grid_int.shape} ({grid_int.size} points)")

    flat_qx = grid_qx.flatten()
    flat_qz = grid_qz.flatten()
    flat_int = grid_int.flatten()

    # --- Step 1: Coarse peak search on smoothed data ---
    print(f"  [Step 1] Gaussian smoothing (sigma={smooth_sigma})...", end=' ')
    smoothed = gaussian_filter(grid_int.astype(float), sigma=smooth_sigma)
    print(f"done ({time.time()-t0:.1f}s)")

    # --- Step 2: STO auto-detection and individual fit ---
    print(f"  [Step 2] STO auto-detection (tol={sto_search_tol})...", end=' ')

    # Find local maxima on full smoothed data for STO detection
    neighbourhood_size = max(3, int(smooth_sigma * 2.5)) | 1
    local_max_mask = (smoothed == maximum_filter(smoothed, size=neighbourhood_size))
    max_indices = np.argwhere(local_max_mask)
    max_intensities = smoothed[local_max_mask]
    sort_order = np.argsort(-max_intensities)
    max_indices = max_indices[sort_order]
    max_intensities = max_intensities[sort_order]
    max_qx = grid_qx[max_indices[:, 0], max_indices[:, 1]]
    max_qz = grid_qz[max_indices[:, 0], max_indices[:, 1]]

    sto_peak = None
    sto_mask = (np.abs(max_qx - STO_QX_REF) < sto_search_tol) & \
               (np.abs(max_qz - STO_QZ_REF) < sto_search_tol)

    if np.any(sto_mask):
        sto_idx = np.where(sto_mask)[0][0]
        sto_coarse = {
            'qx': max_qx[sto_idx], 'qz': max_qz[sto_idx],
            'intensity': max_intensities[sto_idx],
        }
        print(f"FOUND at Qx={sto_coarse['qx']:.5f}, Qz={sto_coarse['qz']:.5f}")

        # Individually fit STO with Pseudo-Voigt
        if refine:
            print(f"  [Step 2b] Fitting STO individually...", end=' ')
            t_sto = time.time()
            dist_sto = np.sqrt((flat_qx - sto_coarse['qx'])**2 +
                               (flat_qz - sto_coarse['qz'])**2)
            sto_fit_mask = dist_sto < fit_radius
            s_qx = flat_qx[sto_fit_mask]
            s_qz = flat_qz[sto_fit_mask]
            s_int = flat_int[sto_fit_mask]
            n_sto = len(s_qx)

            if n_sto > max_fit_points:
                rng = np.random.default_rng(42)
                sub_idx = rng.choice(n_sto, max_fit_points, replace=False)
                s_qx, s_qz, s_int = s_qx[sub_idx], s_qz[sub_idx], s_int[sub_idx]
                n_sto = max_fit_points

            if n_sto >= 10:
                coords = np.array([s_qx, s_qz])
                bg0 = float(np.percentile(s_int, 10))
                p0 = np.array([bg0, sto_coarse['intensity'], sto_coarse['qx'],
                               sto_coarse['qz'], 0.02, 0.02, 0.5])
                lb = np.array([0, 0, sto_coarse['qx'] - 0.05, sto_coarse['qz'] - 0.05,
                               0.001, 0.001, 0.0])
                ub = np.array([sto_coarse['intensity'], sto_coarse['intensity'] * 10,
                               sto_coarse['qx'] + 0.05, sto_coarse['qz'] + 0.05, 0.2, 0.2, 1.0])

                def _single_residual(params, coords, obs):
                    return params[0] + _pseudo_voigt_2d(coords, *params[1:]) - obs

                try:
                    result = least_squares(
                        _single_residual, p0, args=(coords, s_int),
                        bounds=(lb, ub), method='trf',
                        max_nfev=max_nfev, ftol=ftol, xtol=xtol,
                    )
                    fp = result.x
                    sto_peak = {
                        'qx': fp[2], 'qz': fp[3], 'intensity': fp[1],
                        'sigma_x': fp[4], 'sigma_z': fp[5], 'eta': fp[6],
                        'background': fp[0],
                    }
                    print(f"converged={result.success}, Qx={fp[2]:.5f}, "
                          f"Qz={fp[3]:.5f} ({time.time()-t_sto:.2f}s)")
                except Exception as e:
                    print(f"fit failed ({e}), using coarse position")
                    sto_peak = {**sto_coarse, 'sigma_x': 0.02, 'sigma_z': 0.02, 'eta': 0.5, 'background': 0.0}
            else:
                sto_peak = {**sto_coarse, 'sigma_x': 0.02, 'sigma_z': 0.02, 'eta': 0.5, 'background': 0.0}
        else:
            sto_peak = {**sto_coarse, 'sigma_x': 0.02, 'sigma_z': 0.02, 'eta': 0.5, 'background': 0.0}
    else:
        print("NOT FOUND - check sto_search_tol or data range")

    # --- Step 3: Exclude STO region and apply intensity threshold ---
    print(f"  [Step 3] Excluding STO region (r={sto_exclusion_radius}) "
          f"& thresholding (intensity>{intensity_threshold})...", end=' ')

    # Work on the smoothed grid for film peak search
    masked_smoothed = smoothed.copy()

    if sto_peak is not None:
        # Zero out everything within sto_exclusion_radius of the fitted STO position
        dist_to_sto = np.sqrt((grid_qx - sto_peak['qx'])**2 +
                              (grid_qz - sto_peak['qz'])**2)
        masked_smoothed[dist_to_sto < sto_exclusion_radius] = 0.0

    # Apply absolute intensity threshold (same value as threshold_filt in plots)
    masked_smoothed[masked_smoothed < intensity_threshold] = 0.0

    # Apply Qx/Qz search bounds (zero out regions outside specified range)
    if qx_range is not None:
        masked_smoothed[(grid_qx < qx_range[0]) | (grid_qx > qx_range[1])] = 0.0
    if qz_range is not None:
        masked_smoothed[(grid_qz < qz_range[0]) | (grid_qz > qz_range[1])] = 0.0

    print(f"done ({time.time()-t0:.1f}s)")

    # --- Step 4: Find local maxima on masked data ---
    print(f"  [Step 4] Finding film peak maxima (min_sep={min_separation})...", end=' ')
    local_max_film = (masked_smoothed == maximum_filter(masked_smoothed, size=neighbourhood_size))
    local_max_film &= (masked_smoothed > 0)  # exclude zeroed regions
    film_indices = np.argwhere(local_max_film)
    film_max_int = masked_smoothed[local_max_film]

    if len(film_indices) == 0:
        print("no film maxima found")
        film_peaks = []
    else:
        sort_f = np.argsort(-film_max_int)
        film_indices = film_indices[sort_f]
        film_max_int = film_max_int[sort_f]
        film_max_qx = grid_qx[film_indices[:, 0], film_indices[:, 1]]
        film_max_qz = grid_qz[film_indices[:, 0], film_indices[:, 1]]
        film_max_raw = grid_int[film_indices[:, 0], film_indices[:, 1]]

        # Enforce min_separation
        selected = []
        for i in range(len(film_indices)):
            too_close = False
            for j in selected:
                dist = np.sqrt((film_max_qx[i] - film_max_qx[j])**2 +
                               (film_max_qz[i] - film_max_qz[j])**2)
                if dist < min_separation:
                    too_close = True
                    break
            if not too_close:
                selected.append(i)
            if len(selected) >= no_peaks:
                break

        film_peaks = []
        for idx in selected:
            film_peaks.append({
                'qx': film_max_qx[idx], 'qz': film_max_qz[idx],
                'intensity': film_max_raw[idx],
            })
        print(f"found {len(film_peaks)} film peaks ({time.time()-t0:.1f}s)")
        for i, pk in enumerate(film_peaks):
            print(f"    Film peak {i}: Qx={pk['qx']:.5f}, Qz={pk['qz']:.5f}, "
                  f"I={pk['intensity']:.0f}")

    n_film = len(film_peaks)

    # --- Step 5: Simultaneous multi-peak fit of film peaks ---
    refined_film = []
    if refine and n_film > 0:
        print(f"  [Step 5] Simultaneous {n_film}-peak fit  (radius={fit_radius}, "
              f"max_pts={max_fit_points})...", end=' ')
        t_fit = time.time()

        # Collect data within fit_radius of ANY film peak (and outside STO exclusion)
        combined_mask = np.zeros(len(flat_qx), dtype=bool)
        for pk in film_peaks:
            dist = np.sqrt((flat_qx - pk['qx'])**2 + (flat_qz - pk['qz'])**2)
            combined_mask |= (dist < fit_radius)

        # Also exclude STO region from fitting data
        if sto_peak is not None:
            dist_sto = np.sqrt((flat_qx - sto_peak['qx'])**2 +
                               (flat_qz - sto_peak['qz'])**2)
            combined_mask &= (dist_sto >= sto_exclusion_radius)

        fit_qx = flat_qx[combined_mask]
        fit_qz = flat_qz[combined_mask]
        fit_int = flat_int[combined_mask]
        n_pts = len(fit_qx)

        # Subsample if needed
        if n_pts > max_fit_points:
            rng = np.random.default_rng(42)
            sub_idx = rng.choice(n_pts, max_fit_points, replace=False)
            fit_qx, fit_qz, fit_int = fit_qx[sub_idx], fit_qz[sub_idx], fit_int[sub_idx]
            n_pts = max_fit_points

        if n_pts < 10 * n_film:
            print(f"too few points ({n_pts}), using coarse positions")
            for pk in film_peaks:
                refined_film.append({**pk, 'sigma_x': 0.02, 'sigma_z': 0.02, 'eta': 0.5, 'background': 0.0})
        else:
            coords = np.array([fit_qx, fit_qz])
            # Build initial guess: [bg, amp1, x01, z01, sx1, sz1, eta1, amp2, ...]
            bg0 = float(np.percentile(fit_int, 10))
            p0 = [bg0]
            lb = [0.0]
            ub = [float(np.max(fit_int))]
            for pk in film_peaks:
                p0.extend([pk['intensity'], pk['qx'], pk['qz'], 0.02, 0.02, 0.5])
                lb.extend([0, pk['qx'] - 0.05, pk['qz'] - 0.05, 0.001, 0.001, 0.0])
                ub.extend([pk['intensity'] * 10, pk['qx'] + 0.05, pk['qz'] + 0.05,
                           0.2, 0.2, 1.0])
            p0 = np.array(p0)
            lb = np.array(lb)
            ub = np.array(ub)

            try:
                result = least_squares(
                    _multi_peak_residuals, p0, args=(coords, fit_int, n_film),
                    bounds=(lb, ub), method='trf',
                    max_nfev=max_nfev, ftol=ftol, xtol=xtol,
                )
                fp = result.x
                film_background = fp[0]
                for i in range(n_film):
                    pp = fp[1 + i * 6: 1 + (i + 1) * 6]
                    refined_film.append({
                        'qx': pp[1], 'qz': pp[2], 'intensity': pp[0],
                        'sigma_x': pp[3], 'sigma_z': pp[4], 'eta': pp[5],
                        'background': film_background,
                    })
                print(f"converged={result.success}, {n_pts} pts ({time.time()-t_fit:.2f}s)")
                for i, rpk in enumerate(refined_film):
                    print(f"    Film {i}: Qx={rpk['qx']:.5f}, Qz={rpk['qz']:.5f}")
            except Exception as e:
                print(f"fit failed ({e}), using coarse positions")
                for pk in film_peaks:
                    refined_film.append({**pk, 'sigma_x': 0.02, 'sigma_z': 0.02, 'eta': 0.5, 'background': 0.0})
    elif not refine:
        print(f"  [Step 5] Skipping refinement (refine=False)")
        for pk in film_peaks:
            refined_film.append({**pk, 'sigma_x': 0.02, 'sigma_z': 0.02, 'eta': 0.5, 'background': 0.0})

    # --- Step 6: Build output DataFrame ---
    rows = []
    # STO first
    if sto_peak is not None:
        rows.append({
            'peak_id': 0,
            'qx': sto_peak['qx'],
            'qz': sto_peak['qz'],
            'a': _Qx_to_a(sto_peak['qx']),
            'c': _Qz_to_c(sto_peak['qz']),
            'intensity': sto_peak['intensity'],
            'sigma_x': sto_peak.get('sigma_x', 0.02),
            'sigma_z': sto_peak.get('sigma_z', 0.02),
            'eta': sto_peak.get('eta', 0.5),
            'background': sto_peak.get('background', 0.0),
            'is_substrate': True,
            'material': 'STO',
            'sample_name': sample_name,
            'sample_code': sample_code,
            'Ba_doping': None,
            'La_doping': None,
            'fit_radius': fit_radius,
            'sto_exclusion_radius': sto_exclusion_radius,
            'intensity_threshold': intensity_threshold,
        })
    # Film peaks
    for i, pk in enumerate(refined_film):
        rows.append({
            'peak_id': len(rows),
            'qx': pk['qx'],
            'qz': pk['qz'],
            'a': _Qx_to_a(pk['qx']),
            'c': _Qz_to_c(pk['qz']),
            'intensity': pk['intensity'],
            'sigma_x': pk['sigma_x'],
            'sigma_z': pk['sigma_z'],
            'eta': pk['eta'],
            'background': pk.get('background', 0.0),
            'is_substrate': False,
            'material': None,
            'sample_name': sample_name,
            'sample_code': sample_code,
            'Ba_doping': None,
            'La_doping': None,
            'fit_radius': fit_radius,
            'sto_exclusion_radius': sto_exclusion_radius,
            'intensity_threshold': intensity_threshold,
        })

    peaks_df = pd.DataFrame(rows)

    # --- Step 7: Print formatted output ---
    print(f"\n{'='*80}")
    print(f"  Peak Search Results for: {sample_name}  (total time: {time.time()-t0:.1f}s)")
    print(f"{'='*80}")
    print(f"  {'ID':<4} {'Material':<12} {'Qx (2π/Å)':<12} {'Qz (2π/Å)':<12} "
          f"{'a (Å)':<10} {'c (Å)':<10} {'Intensity':<12} {'η':<6}")
    print(f"  {'-'*76}")

    for _, row in peaks_df.iterrows():
        mat_str = row['material'] if row['material'] else '---'
        sub_tag = ' [substrate]' if row['is_substrate'] else ''
        print(f"  {row['peak_id']:<4} {mat_str:<12} {row['qx']:<12.5f} {row['qz']:<12.5f} "
              f"{row['a']:<10.4f} {row['c']:<10.4f} {row['intensity']:<12.1f} {row['eta']:<6.3f}{sub_tag}")

    print(f"{'='*80}")

    # Print known bulk references
    print(f"\n  Bulk references (103):")
    for key, mat in PEROVSKITE_MATERIALS.items():
        a_bulk = mat['a_bulk']
        qx_bulk = 2 * np.pi / a_bulk
        qz_bulk = 6 * np.pi / a_bulk
        print(f"    {key:<10} {mat['name']:<25} a={a_bulk:.3f} Å   "
              f"Qx={qx_bulk:.4f}  Qz={qz_bulk:.4f}")
    print()

    return peaks_df


# =============================================================================
# Peak assignment (interactive)
# =============================================================================
def peak_assign(
    peaks_df: pd.DataFrame,
    materials: Optional[Dict] = None,
) -> pd.DataFrame:
    """
    Interactively assign material labels to fitted peaks.

    For each non-substrate peak, the user is prompted to select a material from a
    numbered list. The DataFrame is updated in-place and returned.

    Parameters:
        peaks_df  : DataFrame returned by peak_search().
        materials : Dict of material keys → info (default: PEROVSKITE_MATERIALS).

    Returns:
        Updated DataFrame with material column filled in.
    """
    if materials is None:
        materials = PEROVSKITE_MATERIALS

    # Build numbered option list (excluding STO which is auto-assigned)
    mat_keys = [k for k in materials.keys() if k != 'STO']
    print("\nAvailable materials:")
    for i, key in enumerate(mat_keys):
        print(f"  {i + 1}: {key} ({materials[key]['name']}, a={materials[key]['a_bulk']:.3f} Å)")
    print(f"  {len(mat_keys) + 1}: Other (enter custom label)")
    print(f"  {len(mat_keys) + 2}: Skip")
    print()

    for idx, row in peaks_df.iterrows():
        if row['is_substrate']:
            continue  # STO already assigned
        if row['material'] is not None:
            continue  # already assigned from a previous run

        print(f"\nPeak {int(row['peak_id'])}: Qx={row['qx']:.5f}, Qz={row['qz']:.5f} | "
              f"a={row['a']:.4f} Å, c={row['c']:.4f} Å | Intensity={row['intensity']:.1f}")

        while True:
            try:
                choice = int(input(f"  Peak {int(row['peak_id'])}: Select material (1-{len(mat_keys) + 2}): "))
                if 1 <= choice <= len(mat_keys):
                    peaks_df.at[idx, 'material'] = mat_keys[choice - 1]
                    break
                elif choice == len(mat_keys) + 1:
                    custom = input("  Enter custom label: ").strip()
                    if custom:
                        peaks_df.at[idx, 'material'] = custom
                    break
                elif choice == len(mat_keys) + 2:
                    break
                else:
                    print(f"  Invalid choice. Enter 1-{len(mat_keys) + 2}.")
            except ValueError:
                print("  Please enter a number.")

    print("\nAssignment complete:")
    for _, row in peaks_df.iterrows():
        mat_str = row['material'] if row['material'] else '---'
        print(f"  Peak {row['peak_id']}: {mat_str} (a={row['a']:.4f}, c={row['c']:.4f})")

    return peaks_df


# =============================================================================
# Peak storage (.xlsx)
# =============================================================================
def peak_store(
    peaks_df: pd.DataFrame,
    xlsx_path: str = 'RSM_peaks.xlsx',
) -> None:
    """
    Save fitted and assigned peaks to an Excel file.

    If the file exists, removes existing rows for the same sample_name before appending
    (so re-running is idempotent). Creates the file if it doesn't exist.

    Parameters:
        peaks_df : DataFrame returned by peak_search() / peak_assign().
        xlsx_path: Path to the .xlsx file.
    """
    xlsx_path = Path(xlsx_path)

    if xlsx_path.exists():
        existing = pd.read_excel(xlsx_path, engine='openpyxl')
        sample_names = peaks_df['sample_name'].unique()
        existing = existing[~existing['sample_name'].isin(sample_names)]
        combined = pd.concat([existing, peaks_df], ignore_index=True)
    else:
        combined = peaks_df.copy()

    combined.to_excel(xlsx_path, index=False, engine='openpyxl')
    print(f"Peaks saved to {xlsx_path} ({len(peaks_df)} rows for "
          f"{peaks_df['sample_name'].iloc[0]})")


def peak_load(
    xlsx_path: str = 'RSM_peaks.xlsx',
    sample_names: Optional[List[str]] = None,
    sample_codes: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load peak data from the shared Excel file.

    Parameters:
        xlsx_path    : Path to the .xlsx file.
        sample_names : Optional list of sample names (plot_string) to filter by.
        sample_codes : Optional list of sample codes (e.g. 'HC021') to filter by (case-insensitive).

    Returns:
        DataFrame of peak data.
    """
    xlsx_path = Path(xlsx_path)
    if not xlsx_path.exists():
        print(f"File not found: {xlsx_path}")
        return pd.DataFrame()

    df = pd.read_excel(xlsx_path, engine='openpyxl')

    if sample_names is not None:
        df = df[df['sample_name'].isin(sample_names)]
    elif sample_codes is not None:
        codes_upper = [c.upper() for c in sample_codes]
        if 'sample_code' in df.columns:
            df = df[df['sample_code'].str.upper().isin(codes_upper)]
        else:
            # Fallback: try matching against sample_name
            df = df[df['sample_name'].str.upper().isin(codes_upper)]

    return df


# =============================================================================
# Peak comparison plots
# =============================================================================
def _resolve_sample_codes(sample_list: list) -> List[str]:
    """
    Resolve a mixed list of RSM objects and/or sample code strings to a list
    of uppercase sample code strings.

    Accepts:
        - RSM class instances (uses .sample_code attribute)
        - Strings like 'HC021', 'em023' (case-insensitive)
    """
    codes = []
    for item in sample_list:
        if hasattr(item, 'sample_code'):
            codes.append(item.sample_code.upper())
        elif isinstance(item, str):
            codes.append(item.upper())
        else:
            codes.append(str(item).upper())
    return codes


def lattice_param_2D(
    sample_list: list,
    xlsx_path: str = 'RSM_peaks.xlsx',
    sto_correction: bool = False,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    export_data: bool = False,
    show_bulk_refs: bool = True,
    show_STO: bool = False,
    material_select: Union[bool, List[str]] = False,
    materials: Optional[Dict] = None,
) -> Figure:
    """
    Scatter plot of in-plane lattice parameter (a) vs out-of-plane (c).

    Each sample gets a unique colour (seaborn colorblind palette).
    Each material type gets a unique marker shape.

    Parameters:
        sample_list      : RSM objects or sample code strings (case-insensitive).
        xlsx_path        : Path to the .xlsx peak database.
        sto_correction   : Shift peaks so each sample's STO aligns with bulk reference.
        x_lim, y_lim     : Axis limits for a, c.
        export_data      : Publication-ready sizing from plot_style.
        show_bulk_refs   : Overlay bulk reference positions (grey crosses) for plotted materials only.
        show_STO         : If True, include STO data points with coloured markers and a bold
                           green vertical reference line.  If False (default), STO is hidden.
        material_select  : False to plot all fitted materials, or a list of material labels
                           (e.g. ['SSO', 'BSSO_40']) to restrict displayed peaks.
        materials        : Material dictionary (default: PEROVSKITE_MATERIALS).

    Returns:
        matplotlib Figure.
    """
    if materials is None:
        materials = PEROVSKITE_MATERIALS

    sample_codes = _resolve_sample_codes(sample_list)
    df = peak_load(xlsx_path, sample_codes=sample_codes)
    if df.empty:
        print("No peak data found for the specified samples.")
        return plt.figure()

    if 'sample_code' not in df.columns:
        df['sample_code'] = df['sample_name']

    # Filter by material_select (keep substrate rows for now — needed for STO correction)
    if material_select is not False and material_select:
        df = df[
            df['material'].apply(lambda m: _material_matches(m, material_select) if pd.notna(m) else False)
            | (df['is_substrate'] == True)
        ]

    # Remove STO rows from plot unless show_STO
    if not show_STO:
        df = df[df['is_substrate'] != True]

    if df.empty:
        print("No peaks remain after filtering.")
        return plt.figure()

    fig_size = set_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    # Colour per sample code
    sample_colors = sns.color_palette("colorblind", len(sample_codes))
    sample_color_map = {code: sample_colors[i] for i, code in enumerate(sample_codes)}

    # STO correction offsets (always from full dataset)
    offsets = {}
    if sto_correction:
        df_full = peak_load(xlsx_path, sample_codes=sample_codes)
        if 'sample_code' not in df_full.columns:
            df_full['sample_code'] = df_full['sample_name']
        for code in sample_codes:
            sto_rows = df_full[(df_full['sample_code'].str.upper() == code) & (df_full['is_substrate'] == True)]
            if not sto_rows.empty:
                sto_row = sto_rows.iloc[0]
                offsets[code] = (3.905 - sto_row['a'], 3.905 - sto_row['c'])
            else:
                offsets[code] = (0.0, 0.0)
                print(f"  Warning: No STO peak for {code}, no correction applied.")

    # Plot each peak
    for _, row in df.iterrows():
        scode = str(row.get('sample_code', row['sample_name'])).upper()
        mat = row['material'] if pd.notna(row['material']) else 'Unknown'
        color = sample_color_map.get(scode, 'grey')
        marker = _get_material_marker(mat)

        a_val, c_val = row['a'], row['c']
        if sto_correction and scode in offsets:
            a_val += offsets[scode][0]
            c_val += offsets[scode][1]

        ax.scatter(a_val, c_val, color=color, marker=marker, s=80,
                   edgecolors='black', linewidths=0.5, zorder=5)

    # Determine which film materials were actually plotted
    plotted_materials = df[df['is_substrate'] != True]['material'].dropna().unique()

    # Bulk reference positions (grey crosses — only for plotted materials)
    if show_bulk_refs:
        for key, mat_info in materials.items():
            if key == 'STO':
                continue
            if key not in plotted_materials:
                continue
            a_bulk = mat_info['a_bulk']
            ax.scatter(a_bulk, a_bulk, marker='x', s=40, color='grey',
                       alpha=0.5, linewidths=0.8, zorder=3)
            ax.annotate(key, (a_bulk, a_bulk), textcoords='offset points',
                        xytext=(4, 4), fontsize=6, color='grey', alpha=0.7)

    # STO vertical reference line (only when show_STO)
    if show_STO:
        ax.axvline(x=3.905, color='green', linestyle='--', linewidth=1.2, alpha=0.7)

    ax.set_xlabel(r"$a_{\parallel}$ (in-plane) $(\mathrm{\AA})$")
    ax.set_ylabel(r"$c_{\perp}$ (out-of-plane) $(\mathrm{\AA})$")

    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    # Dual legend: samples (colour patches) + materials (marker shapes)
    legend_elements = []
    for code in sample_codes:
        legend_elements.append(mpatches.Patch(facecolor=sample_color_map[code],
                                              edgecolor='black', linewidth=0.5, label=code))
    for mat in plotted_materials:
        marker = _get_material_marker(mat)
        legend_elements.append(Line2D([0], [0], marker=marker, color='grey',
                                      markerfacecolor='grey', markersize=6,
                                      linestyle='None', label=mat))
    ax.legend(handles=legend_elements, loc='best', framealpha=0.4, fontsize=7)

    fig.tight_layout()
    plt.show()
    return fig


def lattice_param_1D(
    sample_list: list,
    xlsx_path: str = 'RSM_peaks.xlsx',
    sto_correction: bool = False,
    export_data: bool = False,
    show_bulk_refs: bool = False,
    show_STO: bool = False,
    material_select: Union[bool, List[str]] = False,
    materials: Optional[Dict] = None,
    sample_order: Optional[List[str]] = None,
) -> Figure:
    """
    Strip plots: sample code on x-axis, in-plane (a) and out-of-plane (c) on
    separate y-axes in two subplots.  Each material type gets a unique marker
    and colour.

    Parameters:
        sample_list      : RSM objects or sample code strings (case-insensitive).
        xlsx_path        : Path to .xlsx peak database.
        sto_correction   : Shift peaks so each sample's STO aligns with bulk reference.
        export_data      : Publication-ready sizing from plot_style.
        show_bulk_refs   : Show horizontal lines for bulk lattice parameters (plotted materials only).
        show_STO         : If True, include STO data points and bold green STO reference line.
        material_select  : False to plot all fitted materials, or a list of material labels.
        materials        : Material dictionary (default: PEROVSKITE_MATERIALS).
        sample_order     : Optional explicit ordering for x-axis.

    Returns:
        matplotlib Figure.
    """
    if materials is None:
        materials = PEROVSKITE_MATERIALS

    sample_codes = _resolve_sample_codes(sample_list)
    if sample_order is None:
        sample_order = sample_codes

    df = peak_load(xlsx_path, sample_codes=sample_codes)
    if df.empty:
        print("No peak data found for the specified samples.")
        return plt.figure()

    if 'sample_code' not in df.columns:
        df['sample_code'] = df['sample_name']

    # Filter by material_select
    if material_select is not False and material_select:
        df = df[
            df['material'].apply(lambda m: _material_matches(m, material_select) if pd.notna(m) else False)
            | (df['is_substrate'] == True)
        ]

    # Exclude STO unless show_STO
    if not show_STO:
        df = df[df['is_substrate'] != True]

    if df.empty:
        print("No peaks remain after filtering.")
        return plt.figure()

    fig_size = set_plot_style(export_data=export_data)
    fig, (ax_a, ax_c) = plt.subplots(2, 1, figsize=[fig_size[0], fig_size[1] * 1.6],
                                      sharex=True)

    # Material colours
    mat_colors = sns.color_palette("colorblind", len(materials))
    mat_color_map = {key: mat_colors[i] for i, key in enumerate(materials.keys())}

    # STO correction offsets (from full dataset)
    offsets = {}
    if sto_correction:
        df_full = peak_load(xlsx_path, sample_codes=sample_codes)
        if 'sample_code' not in df_full.columns:
            df_full['sample_code'] = df_full['sample_name']
        for code in sample_codes:
            sto_rows = df_full[(df_full['sample_code'].str.upper() == code) & (df_full['is_substrate'] == True)]
            if not sto_rows.empty:
                sto_row = sto_rows.iloc[0]
                offsets[code] = (3.905 - sto_row['a'], 3.905 - sto_row['c'])
            else:
                offsets[code] = (0.0, 0.0)

    # Map sample codes to x positions
    x_map = {code: i for i, code in enumerate(sample_order)}

    for _, row in df.iterrows():
        scode = str(row.get('sample_code', row['sample_name'])).upper()
        if scode not in x_map:
            continue
        mat = row['material'] if pd.notna(row['material']) else 'Unknown'
        marker = _get_material_marker(mat)
        color = mat_color_map.get(mat, 'grey')
        x_pos = x_map[scode]

        a_val, c_val = row['a'], row['c']
        if sto_correction and scode in offsets:
            a_val += offsets[scode][0]
            c_val += offsets[scode][1]

        ax_a.scatter(x_pos, a_val, color=color, marker=marker, s=80,
                     edgecolors='black', linewidths=0.5, zorder=5)
        ax_c.scatter(x_pos, c_val, color=color, marker=marker, s=80,
                     edgecolors='black', linewidths=0.5, zorder=5)

    # Determine which materials were actually plotted
    plotted_materials = df[df['is_substrate'] != True]['material'].dropna().unique()

    # Bulk reference lines (only for plotted materials, thicker)
    if show_bulk_refs:
        for key, mat_info in materials.items():
            if key == 'STO':
                continue
            if key not in plotted_materials:
                continue
            a_bulk = mat_info['a_bulk']
            color = mat_color_map.get(key, 'grey')
            ax_a.axhline(y=a_bulk, color=color, linestyle='--', linewidth=1.2,
                         alpha=0.7, label=f'{key} bulk')
            ax_c.axhline(y=a_bulk, color=color, linestyle='--', linewidth=1.2,
                         alpha=0.7, label=f'{key} bulk')

    # STO reference line (only when show_STO)
    if show_STO:
        ax_a.axhline(y=3.905, color='green', linestyle='--', linewidth=1.2, alpha=0.7)
        ax_c.axhline(y=3.905, color='green', linestyle='--', linewidth=1.2, alpha=0.7)

    ax_a.set_ylabel(r"$a_{\parallel}$ (in-plane) $(\mathrm{\AA})$")
    ax_c.set_ylabel(r"$c_{\perp}$ (out-of-plane) $(\mathrm{\AA})$")
    ax_c.set_xticks(range(len(sample_order)))
    ax_c.set_xticklabels(sample_order, rotation=45, ha='right', fontsize=7)

    # Material legend (only plotted materials)
    legend_elements = []
    for mat in plotted_materials:
        marker = _get_material_marker(mat)
        color = mat_color_map.get(mat, 'grey')
        legend_elements.append(Line2D([0], [0], marker=marker, color=color,
                                      markerfacecolor=color, markersize=6,
                                      linestyle='None', label=mat))
    ax_a.legend(handles=legend_elements, loc='best', framealpha=0.4, fontsize=7)

    fig.tight_layout()
    plt.show()
    return fig


def lattice_param_doping(
    sample_list: list,
    xlsx_path: str = 'RSM_peaks.xlsx',
    doping_column: str = 'Ba_doping',
    sto_correction: bool = False,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    export_data: bool = False,
    show_bulk_refs: bool = True,
    material_select: Union[bool, List[str]] = False,
    materials: Optional[Dict] = None,
) -> Figure:
    """
    Doping-dependence plot: doping level on x-axis, lattice parameter on y-axis.
    In-plane (a) plotted in red, out-of-plane (c) plotted in blue.

    Parameters:
        sample_list      : RSM objects or sample code strings (case-insensitive).
        xlsx_path        : Path to .xlsx peak database.
        doping_column    : Column name for doping: 'Ba_doping' or 'La_doping'.
        sto_correction   : Shift peaks so each sample's STO aligns with bulk reference.
        x_lim, y_lim     : Axis limits.
        export_data      : Publication-ready sizing from plot_style.
        show_bulk_refs   : Show horizontal lines at bulk lattice parameters for plotted materials.
        material_select  : False to plot all fitted materials, or a list of material labels.
        materials        : Material dictionary (default: PEROVSKITE_MATERIALS).

    Returns:
        matplotlib Figure.
    """
    if materials is None:
        materials = PEROVSKITE_MATERIALS

    sample_codes = _resolve_sample_codes(sample_list)
    df = peak_load(xlsx_path, sample_codes=sample_codes)
    if df.empty:
        print("No peak data found for the specified samples.")
        return plt.figure()

    if 'sample_code' not in df.columns:
        df['sample_code'] = df['sample_name']

    # Ensure doping column exists
    if doping_column not in df.columns:
        print(f"  Column '{doping_column}' not found in the spreadsheet. "
              f"Available columns: {list(df.columns)}")
        return plt.figure()

    # Filter by material_select
    if material_select is not False and material_select:
        df = df[
            df['material'].apply(lambda m: _material_matches(m, material_select) if pd.notna(m) else False)
        ]
    else:
        # Exclude substrate
        df = df[df['is_substrate'] != True]

    # Drop rows without doping value
    df = df.dropna(subset=[doping_column])
    if df.empty:
        print(f"No peaks with '{doping_column}' values found. "
              "Please fill in doping values in the .xlsx file.")
        return plt.figure()

    # STO correction offsets
    offsets = {}
    if sto_correction:
        df_full = peak_load(xlsx_path, sample_codes=sample_codes)
        if 'sample_code' not in df_full.columns:
            df_full['sample_code'] = df_full['sample_name']
        for code in sample_codes:
            sto_rows = df_full[(df_full['sample_code'].str.upper() == code) & (df_full['is_substrate'] == True)]
            if not sto_rows.empty:
                sto_row = sto_rows.iloc[0]
                offsets[code] = (3.905 - sto_row['a'], 3.905 - sto_row['c'])
            else:
                offsets[code] = (0.0, 0.0)

    fig_size = set_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    plotted_materials = df['material'].dropna().unique()

    for _, row in df.iterrows():
        scode = str(row.get('sample_code', row['sample_name'])).upper()
        mat = row['material'] if pd.notna(row['material']) else 'Unknown'
        marker = _get_material_marker(mat, 'o')
        doping = float(row[doping_column])

        a_val, c_val = row['a'], row['c']
        if sto_correction and scode in offsets:
            a_val += offsets[scode][0]
            c_val += offsets[scode][1]

        ax.scatter(doping, a_val, color='red', marker=marker, s=80,
                   edgecolors='black', linewidths=0.5, zorder=5)
        ax.scatter(doping, c_val, color='blue', marker=marker, s=80,
                   edgecolors='black', linewidths=0.5, zorder=5)

    # Bulk reference lines (only plotted materials)
    if show_bulk_refs:
        for key, mat_info in materials.items():
            if key == 'STO':
                continue
            if key not in plotted_materials:
                continue
            a_bulk = mat_info['a_bulk']
            ax.axhline(y=a_bulk, color='grey', linestyle='--', linewidth=1.0,
                       alpha=0.5, label=f'{key} bulk')

    doping_label = doping_column.replace('_', ' ').title()
    ax.set_xlabel(f"{doping_label}")
    ax.set_ylabel(r"Lattice parameter $(\mathrm{\AA})$")

    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    # Legend: in-plane / out-of-plane + material markers
    legend_elements = [
        mpatches.Patch(facecolor='red', edgecolor='black', linewidth=0.5,
                       label=r'$a_{\parallel}$ (in-plane)'),
        mpatches.Patch(facecolor='blue', edgecolor='black', linewidth=0.5,
                       label=r'$c_{\perp}$ (out-of-plane)'),
    ]
    for mat in plotted_materials:
        marker = _get_material_marker(mat, 'o')
        legend_elements.append(Line2D([0], [0], marker=marker, color='grey',
                                      markerfacecolor='grey', markersize=6,
                                      linestyle='None', label=mat))
    ax.legend(handles=legend_elements, loc='best', framealpha=0.4, fontsize=7)

    fig.tight_layout()
    plt.show()
    return fig


def lattice_param_tetragonality(
    sample_list: list,
    xlsx_path: str = 'RSM_peaks.xlsx',
    sto_correction: bool = False,
    export_data: bool = False,
    show_bulk_refs: bool = True,
    material_select: Union[bool, List[str]] = False,
    materials: Optional[Dict] = None,
    sample_order: Optional[List[str]] = None,
    show_unity: bool = True,
    scatter_line: bool = True,
) -> Figure:
    """
    Scatter / line chart of tetragonality (c/a) across samples, one trace per material.

    Tetragonality = c_perp / a_parallel.  A value of 1.0 indicates cubic symmetry;
    >1 means the unit cell is elongated out-of-plane, <1 means compressed.

    Parameters:
        sample_list      : RSM objects or sample code strings (case-insensitive).
        xlsx_path        : Path to .xlsx peak database.
        sto_correction   : Shift peaks so each sample's STO aligns with bulk reference.
        export_data      : Publication-ready sizing from plot_style.
        show_bulk_refs   : Show horizontal lines at bulk c/a for plotted materials.
        material_select  : False to plot all fitted materials, or a list of material labels.
        materials        : Material dictionary (default: PEROVSKITE_MATERIALS).
        sample_order     : Optional explicit ordering for x-axis.
        show_unity       : Draw a horizontal line at c/a = 1 (cubic reference).
        scatter_line     : If True (default) connect markers with a line per material.

    Returns:
        matplotlib Figure.
    """
    if materials is None:
        materials = PEROVSKITE_MATERIALS

    sample_codes = _resolve_sample_codes(sample_list)
    if sample_order is None:
        sample_order = sample_codes

    df = peak_load(xlsx_path, sample_codes=sample_codes)
    if df.empty:
        print("No peak data found for the specified samples.")
        return plt.figure()

    if 'sample_code' not in df.columns:
        df['sample_code'] = df['sample_name']

    # Filter by material_select
    if material_select is not False and material_select:
        df = df[
            df['material'].apply(lambda m: _material_matches(m, material_select) if pd.notna(m) else False)
            | (df['is_substrate'] == True)
        ]

    # Exclude substrate
    df = df[df['is_substrate'] != True]

    if df.empty:
        print("No film peaks remain after filtering.")
        return plt.figure()

    # STO correction offsets (from full dataset)
    offsets = {}
    if sto_correction:
        df_full = peak_load(xlsx_path, sample_codes=sample_codes)
        if 'sample_code' not in df_full.columns:
            df_full['sample_code'] = df_full['sample_name']
        for code in sample_codes:
            sto_rows = df_full[(df_full['sample_code'].str.upper() == code) & (df_full['is_substrate'] == True)]
            if not sto_rows.empty:
                sto_row = sto_rows.iloc[0]
                offsets[code] = (3.905 - sto_row['a'], 3.905 - sto_row['c'])
            else:
                offsets[code] = (0.0, 0.0)

    # Compute tetragonality
    df = df.copy()
    df['scode'] = df['sample_code'].str.upper()
    a_vals = df['a'].values.copy()
    c_vals = df['c'].values.copy()
    if sto_correction:
        for i, row in df.iterrows():
            sc = row['scode']
            if sc in offsets:
                a_vals[df.index.get_loc(i)] += offsets[sc][0]
                c_vals[df.index.get_loc(i)] += offsets[sc][1]
    df['a_corr'] = a_vals
    df['c_corr'] = c_vals
    df['tetragonality'] = df['c_corr'] / df['a_corr']

    # Unique materials present
    plotted_materials = df['material'].dropna().unique()

    # Material colours
    mat_colors = sns.color_palette("colorblind", len(materials))
    mat_color_map = {key: mat_colors[i] for i, key in enumerate(materials.keys())}

    fig_size = set_plot_style(export_data=export_data)
    fig, ax = plt.subplots(figsize=fig_size)

    x_map = {code: i for i, code in enumerate(sample_order)}
    linestyle = '-' if scatter_line else 'None'

    for mat in plotted_materials:
        mat_df = df[df['material'] == mat]
        color = mat_color_map.get(mat, 'grey')
        if color == 'grey':
            for prefix in mat_color_map:
                if mat.startswith(prefix + '_'):
                    color = mat_color_map[prefix]
                    break

        marker = _get_material_marker(mat, 'o')

        # Build ordered x, y arrays following sample_order
        xs, ys = [], []
        for sc in sample_order:
            rows = mat_df[mat_df['scode'] == sc]
            for _, row in rows.iterrows():
                xs.append(x_map[sc])
                ys.append(row['tetragonality'])

        ax.plot(xs, ys, marker=marker, color=color, linestyle=linestyle,
                linewidth=1.2, markersize=7, markeredgecolor='black',
                markeredgewidth=0.5, label=mat, zorder=5)

    # Cubic reference line
    if show_unity:
        ax.axhline(y=1.0, color='black', linestyle='-', linewidth=0.8, alpha=0.5,
                   label='Cubic (c/a = 1)')

    # Bulk tetragonality reference lines (bulk perovskites are cubic, c/a = 1)
    if show_bulk_refs:
        for key, mat_info in materials.items():
            if key == 'STO':
                continue
            if key not in plotted_materials:
                if not any(m.startswith(key + '_') for m in plotted_materials):
                    continue
            ax.axhline(y=1.0, color=mat_color_map.get(key, 'grey'),
                       linestyle='--', linewidth=1.0, alpha=0.5)

    ax.set_ylabel(r"Tetragonality ($c / a$)")
    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=45, ha='right', fontsize=7)

    # Build legend from plotted traces
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, loc='best', framealpha=0.4, fontsize=7)

    fig.tight_layout()
    return fig

