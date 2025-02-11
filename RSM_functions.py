#import all the libraries needed
from import_dep import *
from typing import Optional, Tuple, Any, Dict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.interpolate import griddata

def _plot_on_ax(
    ax: plt.Axes,
    d: pd.DataFrame,
    class_obj: Any,
    count_d: int,
    threshold_filt: float = 0,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal',
    colormap: str = 'plasma',
    label_str: Optional[str] = None,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    contour_color: bool = False,
    show_colorbar: bool = True,
) -> None:
    """
    A helper function to plot a single dataset on the provided axis.

    Parameters:
        ax           : The matplotlib Axes instance to draw on.
        d            : Dictionary containing the data (e.g., 'qx', 'qz', 'Intensity').
        class_obj    : An object that may hold additional data (like lattice parameters).
        count_d      : Index to reference specific dataset properties from class_obj.
        threshold_filt: Intensity threshold for filtering data.
        plot_type    : Type of plot to generate: 'scatter', 'contour', or 'mesh'.
        ax_unit      : Type of axis unit to use, e.g., 'reciprocal' or 'lattice'.
        colormap     : The colormap to use for plotting.
        label_str    : Label for the plot legend.
        x_lim        : Optional x-axis limits.
        y_lim        : Optional y-axis limits.
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
        
        ax.set_xlabel(r"$Q_{\parallel} \,(2\pi/\mathrm{\AA})$", fontsize=16)
        ax.set_ylabel(r"$Q_{\perp} \,(2\pi/\mathrm{\AA})$", fontsize=16)
        ax.set_title(r"Reciprocal Space Map ($Q_x$ vs $Q_z$)", fontsize=14)
    elif ax_unit == 'lattice':
        # df format
        x = class_obj.lat_param_df[count_d]['a']
        y = class_obj.lat_param_df[count_d]['c']
        z = class_obj.lat_param_df[count_d]['Intensity']
        # Numpy format where it is reshaped into meshgrid
        grid_x = class_obj.lat_param_np[count_d][:,:,0]
        grid_y = class_obj.lat_param_np[count_d][:,:,1]
        grid_intensity = class_obj.lat_param_np[count_d][:,:,2]
        
        ax.set_xlabel(r"$a \,(\mathrm{\AA})$", fontsize=16)
        ax.set_ylabel(r"$c \,(\mathrm{\AA})$", fontsize=16)
        ax.set_title(r"Lattice Parameters ($a$ vs $c$)", fontsize=14)
    else:
        raise ValueError("Unknown ax_unit: choose 'reciprocal' or 'lattice'")

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
    grid_intensity[grid_intensity < 100] = 0
    
    # Apply logarithmic scaling to plot_intensity for alpha mapping
    log_plot_intensity = np.log10(plot_intensity + 1e-99)

    # Normalize log_plot_intensity to range [0, 1] for alpha mapping
    norm_log_intensity = (log_plot_intensity - log_plot_intensity.min()) / (log_plot_intensity.max() - log_plot_intensity.min())    
    
    norm_log_intensity[norm_log_intensity > 0.1] = 1
    norm_log_intensity[norm_log_intensity < 0.1] += 0.06
    
    # Invert the normalized intensity to get higher transparency for smaller values
    alpha_values = norm_log_intensity
    

    # Generate plot based on plot_type
    if plot_type == 'scatter':
        p_out = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=colormap, s=2.5, edgecolors='none', label=label_str, norm=LogNorm(),alpha=alpha_values)
    
    elif plot_type == 'contour':
        
        p_out = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=colormap, s=2.5, edgecolors='none', label=label_str, norm=LogNorm(),alpha=alpha_values)
        
        contour_levels = np.logspace(np.log10(1 + 1e-99), np.log10(plot_intensity.max()), 10)
        
        # generate contours from the meshgrid data to prevent errors/streaking
        ax.contour(grid_x, grid_y, grid_intensity, levels=contour_levels, colors=c_color, linewidths=0.3)
        
    elif plot_type == 'mesh':

        
        p_out = ax.pcolormesh(grid_x, grid_y, grid_intensity, cmap=colormap, shading='auto', norm=LogNorm())
        scatter_for_legend = ax.scatter([], [], c=c_color, label=label_str)
        
        # Apply Gaussian filter to the interpolated data for contour plot
        #smoothed_intensity = gaussian_filter(grid_intensity, sigma=0.00001)
        
        contour_levels = np.logspace(np.log10(1 + 1e-99), np.log10(plot_intensity.max()), 10) 
        ax.contour(grid_x, grid_y, grid_intensity, levels=contour_levels, colors=c_color, linewidths=0.2)
        
    elif plot_type == 'surf':
        # Set up the 3D plot
        
        p_out = ax.plot_surface(grid_x, grid_y, grid_intensity, cmap=colormap, linewidth=0, antialiased=False)
        
        
    elif plot_type == 'triangle':
        
        ax.plot_trisurf(x_filt, y_filt, plot_intensity, cmap=colormap, shading='auto', norm=LogNorm())
        # Add a color bar
        p_out = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=colormap)
        
    
    else:
        raise ValueError("Unknown plot_type: choose 'scatter', 'contour', or 'mesh'")        

    # Generate a colorbar for the plot
    if show_colorbar: 
            cbar = plt.colorbar(p_out, ax=ax, pad=0.02, format=LogFormatterSciNotation())  
            cbar.set_label("Intensity", fontsize=12)
            cbar.ax.tick_params(labelsize=10)
            
    # Add guidelines for lattice plots
    if ax_unit == 'lattice':
        ax.scatter(4.14, 4.14, color='red', marker='+', s=100, label='BaSnO3 bulk')
        ax.scatter(4.035, 4.033, color='blue', marker='+', s=100, label='SrSnO3 bulk')
        ax.axvline(x=3.905, color='green', linestyle='--', linewidth=1)
        ax.axhline(y=3.905, color='green', linestyle='--', linewidth=1)

    # Set axis limits if provided
    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    # Tweak tick parameters
    ax.tick_params(axis='both', which='major', labelsize=14, length=4, width=0.6)
    ax.tick_params(axis='both', which='minor', labelsize=14, length=2, width=0.4)

    # Add a legend if a label is provided
    if label_str is not None:
        ax.legend(loc='lower right', framealpha=0.4)


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


def XRR_plot_sep(
    dat: Tuple[Any],
    threshold_filt: float = 0,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal'
) -> Figure:
    """
    Plot multiple datasets in separate subplots.

    Parameters:
        dat: A sequence of objects containing data. Each object is expected to have
             a 'qxqz_df' attribute (a sequence of dictionaries with data) and may have
             a 'plot_string' attribute and 'lat_param_df' for lattice parameters.
        threshold_filt: Intensity threshold used to filter data.
        x_lim: Optional tuple specifying the x-axis limits (min, max).
        y_lim: Optional tuple specifying the y-axis limits (min, max).
        plot_type: Type of plot to generate; one of 'scatter', 'contour', or 'mesh'.
        ax_unit: Unit type for the axes; either 'reciprocal' or 'lattice'.

    Returns:
        The matplotlib Figure object containing the subplots.
    """
    total_plots = sum(len(class_obj.qxqz_df) for class_obj in dat)
    n_rows = total_plots // 2 + total_plots % 2

    fig = plt.figure(figsize=(15, 5.5 * n_rows))
    gs = fig.add_gridspec(n_rows, 2)
    
    subplot_counter = 0
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.qxqz_df):
            ax = fig.add_subplot(gs[subplot_counter // 2, subplot_counter % 2])
            label_str = (class_obj.plot_string[count_d]
                         if hasattr(class_obj, 'plot_string') else None)
            _plot_on_ax(ax, d, class_obj, count_d, threshold_filt,
                        plot_type, ax_unit, colormap= custom_cmap , label_str=label_str,
                        x_lim=x_lim, y_lim=y_lim)
            subplot_counter += 1

    fig.tight_layout()
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    plt.show()
    return fig


def XRR_plot_same(
    dat: Tuple[Any],
    threshold_filt: float = 0,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    plot_type: str = 'scatter',
    ax_unit: str = 'reciprocal',
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

    Returns:
        The matplotlib Figure object containing the plot.
    """
    fig, ax = plt.subplots(figsize=(15, 5.5))
    
    col_maps = ['Purples', 'Greens', 'Oranges', 'Blues', 'Reds', 'Greys']
    col_counter = 0
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.qxqz_df):
            label_str = (class_obj.plot_string[count_d]
                         if hasattr(class_obj, 'plot_string') else None)
            cmap = col_maps[col_counter % len(col_maps)]
            _plot_on_ax(ax, d, class_obj, count_d, threshold_filt,
                        plot_type, ax_unit, colormap=cmap, label_str=label_str,
                        x_lim=x_lim, y_lim=y_lim, contour_color=True, show_colorbar=False)
            col_counter += 1

    fig.tight_layout()
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    plt.show()
    return fig

