#import all the libraries needed
from import_dep import *

def _plot_ttw_on_ax(
    ax: plt.Axes,
    d: pd.DataFrame,
    label_str: Optional[str] = None,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    linear_subtraction: bool = False,
    log_plot: bool = True,
    med_filt: int = 0,
) -> None:
    """
    A helper function to plot two theta vs intensity on the provided axis.

    Parameters:
        ax                : The matplotlib Axes instance to draw on.
        d                 : DataFrame containing the data (e.g., 'Angle', 'Intensity').
        label_str         : Label for the plot legend.
        x_lim             : Optional x-axis limits.
        y_lim             : Optional x-axis limits.
        linear_subtraction: Whether to subtract a linear fit from intensity data.
        log_plot          : Whether to use logarithmic scale for y-axis.
        med_filt          : Kernel size for median filter. If 0, no filtering is applied.
    """
    x = d['Angle']
    y = d['Intensity']

    if linear_subtraction:
        # Fit a linear model to the data
        coeffs = np.polyfit(x, y, 1)
        linear_fit = np.poly1d(coeffs)
        
        # Subtract the linear fit from the intensity data
        y = y - linear_fit(x)
    
    # Apply median filter if kernel size > 0
    if med_filt > 0:
        from scipy.signal import medfilt
        y = medfilt(y, kernel_size=med_filt)

    ax.plot(x, y, label=label_str)
    ax.set_xlabel(r"$2\theta \,(\mathrm{degrees})$")
    ax.set_ylabel(r"Intensity")
    #ax.set_title(r"Two Theta vs Intensity")
    
    ## Add reference dotted lines for the bulk phase in the positions they would be expected
    ax.axvline(x=24.37, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(x=23.58, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(x=23.11, color='grey', linestyle='--', alpha=0.5)
    

    if log_plot:
        ax.set_yscale('log')  # Set y-axis to logarithmic scale
    
    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    if label_str is not None:
        ax.legend()


def ttw_plot_sep(
    dat: Tuple[Any],
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    linear_subtraction: bool = False,
    log_plot: bool = True,
    med_filt: int = 0,
) -> Figure:
    """
    Plot multiple datasets in separate subplots.

    Parameters:
        dat               : A sequence of objects containing data. Each object is expected to have
                           a 'ttw_df' attribute (a sequence of dictionaries with data) and may have
                           a 'plot_string' attribute.
        x_lim             : Optional tuple specifying the x-axis limits (min, max).
        y_lim             : Optional tuple specifying the y-axis limits (min, max).
        linear_subtraction: Whether to subtract a linear fit from intensity data.
        log_plot          : Whether to use logarithmic scale for y-axis.
        med_filt          : Kernel size for median filter. If 0, no filtering is applied.

    Returns:
        The matplotlib Figure object containing the subplots.
    """

    total_plots = sum(len(class_obj.ttw_df) for class_obj in dat)
    n_rows = total_plots // 2 + total_plots % 2

    fig = plt.figure(figsize=(15, 5.5 * n_rows))
    gs = fig.add_gridspec(n_rows, 2)
    
    subplot_counter = 0
   
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.ttw_df):
            ax = fig.add_subplot(gs[subplot_counter // 2, subplot_counter % 2])
            label_str = (class_obj.plot_string[count_d]
                         if hasattr(class_obj, 'plot_string') else None)
            _plot_ttw_on_ax(
                ax, 
                d, 
                label_str=label_str, 
                x_lim=x_lim, 
                y_lim=y_lim,
                linear_subtraction=linear_subtraction,
                log_plot=log_plot,
                med_filt=med_filt
            )
            subplot_counter += 1
            
    plt.show()
    return fig


def ttw_plot_same(
    dat: Tuple[Any],
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    norm_plot: bool = False,
    linear_subtraction: bool = False,
    log_plot: bool = True,
    med_filt: int = 0,
) -> Figure:
    """
    Plot multiple datasets on a single set of axes, with an option to normalize to their own maximum intensity.

    Parameters:
        dat               : A sequence of objects containing data. Each object is expected to have
                           a 'ttw_df' attribute (a sequence of dictionaries with data) and may have
                           a 'plot_string' attribute.
        x_lim             : Optional tuple specifying the x-axis limits (min, max).
        y_lim             : Optional tuple specifying the y-axis limits (min, max).
        norm_plot         : Boolean indicating whether to normalize the data to their own maximum intensity.
        linear_subtraction: Whether to subtract a linear fit from intensity data.
        log_plot          : Whether to use logarithmic scale for y-axis.
        med_filt          : Kernel size for median filter. If 0, no filtering is applied.

    Returns:
        The matplotlib Figure object containing the plot.
    """
    fig, ax = plt.subplots()
    
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.ttw_df):
            if norm_plot:
                # Normalize intensity by the maximum intensity of the current dataset
                max_intensity = d['Intensity'].max()
                d['Normalized_Intensity'] = d['Intensity'] / max_intensity
                plot_data = d.assign(Intensity=d['Normalized_Intensity'])  # Use normalized data
            else:
                plot_data = d  # Use raw data
            
            label_str = (class_obj.plot_string[count_d]
                         if hasattr(class_obj, 'plot_string') else None)
            _plot_ttw_on_ax(
                ax, 
                plot_data, 
                label_str=label_str, 
                x_lim=x_lim, 
                y_lim=y_lim,
                linear_subtraction=linear_subtraction,
                log_plot=log_plot,
                med_filt=med_filt
            )
    
    ax.set_ylabel("Intensity (a.u.)")  # Update y-axis label
    ax.set_yticks([])  # Remove y-axis ticks
    ax.yaxis.set_ticklabels([])  # Remove y-axis numbers
    plt.show()
    return fig


def update_plot_string(ttw):
    '''Update the plot_str variable for each file in the tuple.'''
    # Loop over the samples in the multiple ttw objects
    for t in ttw: 
        # Loop over the files within each sample folder
        for i in range(len(t.file_name)):
            new_plot_str = input(f"Enter new plot string for {t.file_name[i]} (current: {t.plot_string[i]}): ")
            t.plot_string[i] = new_plot_str if new_plot_str else t.plot_string[i] #ppms.material + ' - ' + 
            print(f"New plot string for {t.file_name[i]}: {t.plot_string[i]}")
    return ttw
