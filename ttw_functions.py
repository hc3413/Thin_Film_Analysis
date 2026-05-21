#import all the libraries needed
from import_dep import *

def calculate_002_peak_position(lattice_param, wavelength=1.5406):
    """
    Calculate the 2θ position for the 002 peak given a lattice parameter.
    
    Parameters:
        lattice_param : float
            Lattice parameter in Angstroms
        wavelength : float
            X-ray wavelength in Angstroms (default: Cu Kα = 1.5406 Å)
    
    Returns:
        float : 2θ position in degrees
    """
    d_spacing = lattice_param / 2  # For 002 reflection: d = a/2
    theta_rad = np.arcsin(wavelength / (2 * d_spacing))  # Bragg's law: λ = 2d sinθ
    two_theta_deg = 2 * np.degrees(theta_rad)
    return two_theta_deg

def _plot_ttw_on_ax(
    ax: plt.Axes,
    d: pd.DataFrame,
    label_str: Optional[str] = None,
    x_lim: Optional[Tuple[float, float]] = None,
    y_lim: Optional[Tuple[float, float]] = None,
    linear_subtraction: bool = False,
    log_plot: bool = True,
    med_filt: int = 0,
    offset_factor: float = 1.0,
    STO_align: bool = False,
    STO_align_showfit: bool = False,
    bulk_params: Union[bool, List[str]] = False,
    first_index: bool = True,
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
        offset_factor     : Multiplicative factor to offset the plot. For log scale, this creates
                           uniform visual spacing (e.g., 10 = 1 decade offset).
        STO_align         : If True, searches for the STO peak (002, 001, 003, 004) and aligns
                           the x-axis to the literature value.
        STO_align_showfit : If True, shows the fitted peak as a dotted line.
    """
    x = d['Angle'].copy()
    y = d['Intensity'].copy()

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
        
    if STO_align:
        from scipy.optimize import curve_fit
        
        def pseudo_voigt(x_val, amp, cen, wid, eta, bkg):
            # Pseudo-Voigt profile (linear combination of Lorentzian and Gaussian)
            # wid is the Half-Width at Half-Maximum (HWHM)
            # eta is the mixing parameter (1 = pure Lorentzian, 0 = pure Gaussian)
            L = 1 / (1 + ((x_val - cen) / wid)**2)
            G = np.exp(-np.log(2) * ((x_val - cen) / wid)**2)
            return amp * (eta * L + (1 - eta) * G) + bkg

        lattice_STO = 3.905
        wavelength = 1.5406
        
        def get_STO_pos(l_idx):
            d_spacing = lattice_STO / l_idx
            val = wavelength / (2 * d_spacing)
            if val > 1: return None
            return 2 * np.degrees(np.arcsin(val))
            
        fit_plot_data = None
        for reflection in [2, 1, 3, 4]:
            theo_pos = get_STO_pos(reflection)
            if theo_pos is None: continue
            
            mask = (x >= theo_pos - 0.5) & (x <= theo_pos + 0.5)
            x_window = x[mask].values
            y_window = y[mask].values
            
            if len(x_window) < 10:
                continue
                
            amp_guess = np.max(y_window) - np.min(y_window)
            cen_guess = x_window[np.argmax(y_window)]
            wid_guess = 0.05
            eta_guess = 0.5
            bkg_guess = np.median(y_window)  # better than min against noise
            
            try:
                popt, pcov = curve_fit(
                    pseudo_voigt, x_window, y_window, 
                    p0=[amp_guess, cen_guess, wid_guess, eta_guess, bkg_guess],
                    bounds=([0, theo_pos-0.5, 0.0001, 0.0, -np.inf], 
                            [np.inf, theo_pos+0.5, 2.0, 1.0, np.inf]),
                    maxfev=10000
                )
                
                fitted_cen = popt[1]
                shift = theo_pos - fitted_cen
                x = x + shift
                
                if STO_align_showfit:
                    x_fit = np.linspace(np.min(x_window), np.max(x_window), 500)
                    y_fit = pseudo_voigt(x_fit, *popt)
                    # We store the shifted x_fit and the y_fit so we can plot it later
                    fit_plot_data = (x_fit + shift, y_fit)
                
                break # Break if we find and fit a peak
            except Exception as e:
                print(f"Pseudo-Voigt failed: {e}. Retrying with Lorentzian fallback.")
                try:
                    def lorentzian(x_val, amp, cen, wid, bkg):
                        return amp * wid**2 / ((x_val - cen)**2 + wid**2) + bkg
                        
                    popt_fb, pcov_fb = curve_fit(
                        lorentzian, x_window, y_window,
                        p0=[amp_guess, cen_guess, wid_guess, bkg_guess],
                        bounds=([0, theo_pos-0.5, 0.0001, -np.inf], 
                                [np.inf, theo_pos+0.5, 2.0, np.inf]),
                        maxfev=10000
                    )
                    
                    popt = [popt_fb[0], popt_fb[1], popt_fb[2], 1.0, popt_fb[3]]
                    
                    fitted_cen = popt[1]
                    shift = theo_pos - fitted_cen
                    x = x + shift
                    
                    if STO_align_showfit:
                        x_fit = np.linspace(np.min(x_window), np.max(x_window), 500)
                        y_fit = pseudo_voigt(x_fit, *popt)
                        fit_plot_data = (x_fit + shift, y_fit)
                    
                    break
                except Exception as e2:
                    print(f"Fallback fit failed: {e2}")
                    pass
    
    # Apply offset (multiplicative for log scale compatibility)
    y = y * offset_factor

    ax.plot(x, y, label=label_str)
    
    if STO_align and STO_align_showfit and 'fit_plot_data' in locals() and fit_plot_data is not None:
        ax.plot(fit_plot_data[0], fit_plot_data[1] * offset_factor, color='black', linestyle=':', label='STO Fit')

    ax.set_xlabel(r"$2\theta \,(\mathrm{degrees})$")
    # Default y-axis label (may be overridden in ttw_plot_same if offset is used)
    ax.set_ylabel(r"Intensity (a.u.)")
    
    if bulk_params is not False and bulk_params:
        # Define materials database for 002 peaks
        PEROVSKITE_MATERIALS = {
            'STO': {'name': 'SrTiO3', 'a_bulk': 3.905, 'color': 'green'},
            'BSO': {'name': 'BaSnO3', 'a_bulk': 4.116, 'color': 'purple'},
            'SSO': {'name': 'SrSnO3', 'a_bulk': 4.036, 'color': 'orange'},
            'LSO': {'name': 'LaScO3', 'a_bulk': 4.052, 'color': 'gold'},
            'BTO': {'name': 'BaTiO3', 'a_bulk': 4.010, 'color': 'cyan'},
            'BTO_TETRAGONAL': {'name': 'BaTiO3 (Tetragonal)', 'a_bulk': 3.994, 'c_bulk': 4.0335, 'color': 'cyan'},
        }
        
        mats_to_plot = bulk_params if isinstance(bulk_params, list) else ['STO', 'BSO', 'SSO']
        for mat in mats_to_plot:
            for key in PEROVSKITE_MATERIALS:
                if mat.upper() == key.upper():
                    info = PEROVSKITE_MATERIALS[key]
                    
                    # If it's a tetragonal material, it has an explicit c_bulk we might want bounds for.
                    # Since this is a 002 scan, c_bulk is typically what's being probed.
                    # But we can plot both the typical a_bulk and c_bulk positions as reference bounds.
                    for param_key in ['a_bulk', 'c_bulk']:
                        if param_key in info:
                            peak_pos = calculate_002_peak_position(info[param_key])
                            # Do not add labels to the legend key for these reference lines
                            ax.axvline(
                                x=peak_pos, color=info['color'], linestyle='--', alpha=0.5
                            )
                    break


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
    STO_align: bool = False,
    STO_align_showfit: bool = False,
    bulk_params: Union[bool, List[str]] = False,
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
                med_filt=med_filt,
                offset_factor=1.0,  # No offset for separate plots
                STO_align=STO_align,
                STO_align_showfit=STO_align_showfit,
                bulk_params=bulk_params,
                first_index=(count_d == 0)
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
    offset: Optional[float] = None,
    y_label: bool = False,
    STO_align: bool = False,
    STO_align_showfit: bool = False,
    bulk_params: Union[bool, List[str]] = False,
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
        offset            : Multiplicative offset factor between successive plots. In log scale, each
                           successive plot is multiplied by offset^n (e.g., offset=10 gives 1 decade
                           spacing, offset=100 gives 2 decades). If None, no offset is applied.

    Returns:
        The matplotlib Figure object containing the plot.
    """
    fig, ax = plt.subplots()
    
    plot_counter = 0  # Track which plot we're on for offsetting
    
    for class_obj in dat:
        for count_d, d in enumerate(class_obj.ttw_df):
            if norm_plot:
                # Normalize intensity by the maximum intensity of the current dataset
                max_intensity = d['Intensity'].max()
                d['Normalized_Intensity'] = d['Intensity'] / max_intensity
                plot_data = d.assign(Intensity=d['Normalized_Intensity'])  # Use normalized data
            else:
                plot_data = d  # Use raw data
            
            # Calculate offset factor for this plot
            if offset is not None:
                offset_factor = offset ** plot_counter
            else:
                offset_factor = 1.0
            
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
                med_filt=med_filt,
                offset_factor=offset_factor,
                STO_align=STO_align,
                STO_align_showfit=STO_align_showfit,
                bulk_params=bulk_params,
                first_index=(plot_counter == 0)
            )
            
            plot_counter += 1  # Increment for next plot
    
    # Override y-axis label when offset is applied
    if offset is not None:
        # With offset, absolute values are meaningless
        # Add "Log" prefix for log scale since tick marks are removed (not visually obvious)
        if y_label == True:
            if log_plot:
                ax.set_ylabel(r"Log Intensity (a.u., offset)")
            else:
                ax.set_ylabel(r"Intensity (a.u., offset)")
        ax.set_yticks([])  # Remove y-axis ticks
        ax.yaxis.set_ticklabels([])  # Remove y-axis numbers
    
    plt.show()
    return fig


def update_plot_string(ttw):
    '''Update the plot_str variable for each file in the tuple.'''
    # Loop over the samples in the multiple ttw objects
    for t in ttw: 
        # Convert tuple to list for modification
        plot_string_list = list(t.plot_string)
        
        # Loop over the files within each sample folder
        for i in range(len(t.file_name)):
            new_plot_str = input(f"Enter new plot string for {t.file_name[i]} (current: {t.plot_string[i]}): ")
            plot_string_list[i] = new_plot_str if new_plot_str else t.plot_string[i] #ppms.material + ' - ' + 
            print(f"New plot string for {t.file_name[i]}: {plot_string_list[i]}")
        
        # Convert back to tuple and assign
        t.plot_string = tuple(plot_string_list)
        
    return ttw
