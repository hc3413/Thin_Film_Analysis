#import all the libraries needed
from import_dep import *

def _plot_ttw_on_ax(
    ax: plt.Axes,
    d: pd.DataFrame,
    label_str: Optional[str] = None,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
) -> None:
    """
    A helper function to plot two theta vs intensity on the provided axis.

    Parameters:
        ax       : The matplotlib Axes instance to draw on.
        d        : DataFrame containing the data (e.g., 'Angle', 'Intensity').
        label_str: Label for the plot legend.
        x_lim    : Optional x-axis limits.
        y_lim    : Optional y-axis limits.
    """
    x = d['Angle']
    y = d['Intensity']

    ax.plot(x, y, label=label_str)
    ax.set_xlabel(r"$2\theta \,(\mathrm{degrees})$", fontsize=16)
    ax.set_ylabel(r"Intensity", fontsize=16)
    ax.set_title(r"Two Theta vs Intensity", fontsize=14)

    ax.set_yscale('log')  # Set y-axis to logarithmic scale
    
    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    ax.tick_params(axis='both', which='major', labelsize=14, length=4, width=0.6)
    ax.tick_params(axis='both', which='minor', labelsize=14, length=2, width=0.4)

    if label_str is not None:
        ax.legend(loc='upper right', framealpha=0.4)


def ttw_plot_sep(
    dat: Tuple[Any],
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
) -> Figure:
    """
    Plot multiple datasets in separate subplots.

    Parameters:
        dat  : A sequence of objects containing data. Each object is expected to have
               a 'ttw_df' attribute (a sequence of dictionaries with data) and may have
               a 'plot_string' attribute.
        x_lim: Optional tuple specifying the x-axis limits (min, max).
        y_lim: Optional tuple specifying the y-axis limits (min, max).

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
            _plot_ttw_on_ax(ax, d, label_str=label_str, x_lim=x_lim, y_lim=y_lim)
            subplot_counter += 1
            
    fig.tight_layout()
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    plt.show()
    return fig


def ttw_plot_same(
    dat: Tuple[Any],
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
) -> Figure:
    """
    Plot multiple datasets on a single set of axes.

    Parameters:
        dat  : A sequence of objects containing data. Each object is expected to have
               a 'ttw_df' attribute (a sequence of dictionaries with data) and may have
               a 'plot_string' attribute.
        x_lim: Optional tuple specifying the x-axis limits (min, max).
        y_lim: Optional tuple specifying the y-axis limits (min, max).

    Returns:
        The matplotlib Figure object containing the plot.
    """
    fig, ax = plt.subplots(figsize=(15, 5.5))
    
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.ttw_df):
            label_str = (class_obj.plot_string[count_d]
                         if hasattr(class_obj, 'plot_string') else None)
            _plot_ttw_on_ax(ax, d, label_str=label_str, x_lim=x_lim, y_lim=y_lim)
    
    fig.tight_layout()
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    plt.show()
    return fig

