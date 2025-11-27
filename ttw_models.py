import numpy as np
import time
from scipy.signal import convolve
from scipy.optimize import minimize, differential_evolution

# ========== Core Functions: system generation, XRD calc, smoothing ==========

def smooth(y, window=5):
    """Simple moving-average smoothing (edge-preserving via 'same' conv)."""
    if window is None or window <= 1:
        return y
    kernel = np.ones(window) / window
    return convolve(y, kernel, mode='same')


def create_crystal_system(nx, ny, nz,
                          lattice_perfect,
                          lattice_gradient,
                          mode='perfect',
                          gradient_type='linear',
                          spatial_var_percent=0,
                          roughness_layers=0,
                          mosaicity_fwhm=0,
                          relaxation_length=15):
    """
    Generate lattice (per-layer c), height_map and tilt array for a set of crystallites.
    - nx,ny: lateral sampling of columns
    - nz: number of unit-cell layers in z
    - lattice_perfect: base c (Å) used for simplest 'perfect' mode
    - lattice_gradient: [c_interface, c_bulk] used by gradient/exponential profiles
    """
    # --- Lattice ---
    if mode == 'perfect':
        lattice_array = np.full((nx, ny, nz), lattice_perfect)
    elif mode == 'gradient' or mode == 'full':
        layer_indices = np.arange(nz)
        if gradient_type == 'linear':
            gradient_per_layer = (lattice_gradient[1] - lattice_gradient[0]) / float(max(nz - 1, 1))
            c_values = lattice_gradient[0] + gradient_per_layer * layer_indices
        elif gradient_type == 'exponential':
            total_change = lattice_gradient[1] - lattice_gradient[0]
            decay = 1 - np.exp(-layer_indices / float(relaxation_length))
            c_values = lattice_gradient[0] + total_change * decay
        else:
            raise ValueError(f"Unknown gradient_type: {gradient_type}")

        if mode == 'gradient':
            lattice_array = np.broadcast_to(c_values[np.newaxis, np.newaxis, :], (nx, ny, nz)).copy()
        else:
            # mode == 'full' -> allow spatial variation at z=0
            # FIX: Don't add absolute values! Add the VARIATION.
            base_c = lattice_gradient[0]
            sigma = base_c * (spatial_var_percent / 100.0)
            
            # Generate random starting values for each column
            c_start = np.random.normal(base_c, sigma, (nx, ny))
            
            # Calculate the shift for each column relative to the base profile
            c_shift = c_start - base_c
            
            # Apply shift to the gradient profile
            lattice_array = c_shift[:, :, np.newaxis] + c_values[np.newaxis, np.newaxis, :]
    else:
        raise ValueError(f"Unknown mode: {mode}")

    # --- Height / Roughness ---
    # Decoupled from 'mode': Apply if parameter is non-zero
    if roughness_layers > 0:
        heights = np.random.normal(nz, roughness_layers, (nx, ny))
        heights = np.clip(heights, 1, nz * 1.5)
        height_map = np.round(heights).astype(int)
    else:
        height_map = np.full((nx, ny), nz, dtype=int)

    # --- Mosaicity / Tilt ---
    # Decoupled from 'mode': Apply if parameter is non-zero
    if mosaicity_fwhm > 0:
        mosaicity_sigma = mosaicity_fwhm / 2.355
        tilt_map = np.random.normal(0.0, mosaicity_sigma, (nx, ny))
        tilt_array = np.broadcast_to(tilt_map[:, :, np.newaxis], (nx, ny, nz))
    else:
        tilt_array = np.zeros((nx, ny, nz))

    return lattice_array, height_map, tilt_array


def calculate_xrd_vectorized(lattice_array, height_map, tilt_array,
                             two_theta_range=(40,46), wavelength=1.5406,
                             attenuation_length_nm=0, instrumental_broadening=0,
                             n_points=3000, chunk_size=2000, normalize=True):
    """
    Compute XRD pattern from columns with optional absorption and instrumental broadening.
    Returns (two_theta, intensity_normalized).
    """
    start_time = time.time()
    nx, ny, nz = lattice_array.shape
    n_columns = nx * ny

    lattice_flat = lattice_array.reshape(n_columns, nz)
    tilt_flat = tilt_array.reshape(n_columns, nz)
    height_flat = height_map.reshape(n_columns)

    # z positions: accumulate layer thicknesses (approx c per layer -> z in Å)
    z_positions = np.cumsum(lattice_flat, axis=1)
    z_positions = np.concatenate([np.zeros((n_columns, 1)), z_positions[:, :-1]], axis=1)

    # mask of present layers per column
    idx_layer = np.arange(nz)[np.newaxis, :]
    mask = idx_layer < height_flat[:, np.newaxis]

    two_theta = np.linspace(two_theta_range[0], two_theta_range[1], n_points)
    theta = np.radians(two_theta / 2.0)
    Q_base = (4.0 * np.pi / wavelength) * np.sin(theta)

    # mosaic tilt factor (approx cos(tilt))
    tilt_rad = np.radians(tilt_flat)
    mosaicity_factor = np.cos(tilt_rad)

    # absorption: compute using a representative theta (use mean of range)
    theta_mean = np.radians((two_theta_range[0] + two_theta_range[1]) / 4.0)
    if attenuation_length_nm and attenuation_length_nm > 0:
        attenuation_A = attenuation_length_nm * 10.0  # nm -> Å
        path_length = 2.0 * z_positions / np.sin(theta_mean)
        absorption_factor = np.exp(-path_length / attenuation_A)
    else:
        absorption_factor = np.ones_like(z_positions)

    intensity_total = np.zeros(n_points, dtype=float)

    for start in range(0, n_columns, chunk_size):
        end = min(start + chunk_size, n_columns)
        z_chunk = z_positions[start:end, :]
        mask_chunk = mask[start:end, :]
        mos_chunk = mosaicity_factor[start:end, :]
        abs_chunk = absorption_factor[start:end, :]

        # phases: shape (n_points, n_columns_chunk, nz)
        phases = Q_base[:, None, None] * mos_chunk[None, :, :] * z_chunk[None, :, :]
        waves = np.exp(1j * phases) * abs_chunk[None, :, :]
        waves *= mask_chunk[None, :, :]

        F_col = np.sum(waves, axis=2)  # (n_points, n_columns_chunk)
        I_col = np.abs(F_col) ** 2
        intensity_total += np.sum(I_col, axis=1)

    # instrumental broadening
    if instrumental_broadening and instrumental_broadening > 0:
        # FIX: Calculate sigma in points correctly
        step_size = two_theta[1] - two_theta[0]
        sigma_pts = (instrumental_broadening / step_size) / 2.355
        
        # guard against extremely small sigma
        sigma_pts = max(sigma_pts, 0.5)
        
        kernel_radius = int(min(n_points // 2 - 1, max(1, int(4 * sigma_pts))))
        x_kernel = np.arange(-kernel_radius, kernel_radius + 1)
        kernel = np.exp(-x_kernel ** 2 / (2.0 * sigma_pts ** 2))
        kernel /= kernel.sum()
        intensity_total = convolve(intensity_total, kernel, mode='same')

    # normalize
    if normalize and intensity_total.max() > 0:
        intensity_total = intensity_total / intensity_total.max()

    elapsed = time.time() - start_time
    # print(f"XRD calc finished in {elapsed:.3f}s (columns={n_columns})")
    return two_theta, intensity_total


def run_xrd_fitting(ttw_object, peak_intensity, wavelength, 
                    fit_range=(42.0, 45.5), 
                    nx_fit=8, ny_fit=8, 
                    use_log_residual=True,
                    p0=None, bounds=None,
                    tol=0.01, maxiter=None,
                    background_intensity=None,
                    film_scale=1.0,
                    use_global_search=True,
                    de_maxiter=50,
                    de_popsize=10):
    """
    Run the XRD fitting optimization with optional global search.
    
    Parameters:
    -----------
    use_global_search : bool, default=True
        If True, use Differential Evolution (global) followed by Powell (local refinement).
        If False, use only Powell starting from p0.
    de_maxiter : int, default=50
        Maximum iterations for Differential Evolution (global search phase).
        Typical: 30-100. Higher = more thorough but slower.
    de_popsize : int, default=10
        Population size multiplier for DE (actual pop = de_popsize * n_params).
        Typical: 10-15. Higher = better exploration but slower.
    tol : float, default=0.01
        Tolerance for local refinement (Powell). Lower = more precise.
    maxiter : int or None
        Maximum iterations for local refinement (Powell). None = unlimited.
    background_intensity : array or None
        Optional background to add to simulation.
    film_scale : float, default=1.0
        Scaling factor for simulated film intensity.
    """
    
    print(f"Setting up optimization in range {fit_range} with grid {nx_fit}x{ny_fit}...")
    bg_msg = "Enabled" if background_intensity is not None else "None"
    search_mode = "Global (DE) + Local (Powell)" if use_global_search else "Local (Powell only)"
    print(f"Optimization mode: {search_mode}")
    print(f"  Settings: film_scale={film_scale:.4f}, background={bg_msg}")
    if use_global_search:
        print(f"  DE params: maxiter={de_maxiter}, popsize={de_popsize} (total evals ~{de_maxiter * de_popsize * len(bounds)})")
    print(f"  Powell params: tol={tol}, maxiter={maxiter}")
    if background_intensity is not None:
        print("  -> Fitting target: (Simulated_Film * Scale) + Background_Intensity")

    # Extract and normalize data from ttw object
    df = ttw_object.ttw_df[0]
    x_exp = df['Angle'].values
    y_exp = df['Intensity'].values
    y_exp_norm = y_exp / peak_intensity

    # Prepare Experimental Data for Fitting (Subset)
    mask_fit = (x_exp >= fit_range[0]) & (x_exp <= fit_range[1])
    x_fit = x_exp[mask_fit]
    y_fit = y_exp_norm[mask_fit]

    # Define Objective Function
    iteration_count = [0]  # Use list to make it mutable in nested function
    def xrd_objective(params):
        iteration_count[0] += 1
        """
        Objective function for optimization.
        Returns the sum of squared errors between simulation and experiment.
        """
        # Unpack parameters
        c_start, c_end, t_nm, s_var, r_nm, mos, rel_len, att_len, instr = params
        
        # Derived parameters
        nz_fit_layer = int((t_nm * 10) / c_end)
        r_layers = (r_nm * 10) / c_end
        
        # Generate System
        lat_arr, h_map, t_arr = create_crystal_system(
            nx_fit, ny_fit, nz_fit_layer,
            c_end,
            [c_start, c_end],
            mode='full',
            gradient_type='exponential',
            spatial_var_percent=s_var,
            roughness_layers=r_layers,
            mosaicity_fwhm=mos,
            relaxation_length=rel_len
        )
        
        # Calculate XRD
        tt_sim, int_sim = calculate_xrd_vectorized(
            lat_arr, h_map, t_arr,
            two_theta_range=fit_range,
            wavelength=wavelength,
            attenuation_length_nm=att_len,
            instrumental_broadening=instr,
            n_points=1000,
            normalize=True
        )
        
        # Scale the film simulation
        int_sim = int_sim * film_scale

        # Add background if provided
        if background_intensity is not None:
            int_sim = int_sim + background_intensity
        
        # Interpolate to experimental x-axis
        int_sim_interp = np.interp(x_fit, tt_sim, int_sim)
        
        # Calculate Residual
        if use_log_residual:
            epsilon = 1e-6
            residual = np.sum((np.log10(y_fit + epsilon) - np.log10(int_sim_interp + epsilon))**2)
        else:
            residual = np.sum((y_fit - int_sim_interp)**2)
        
        print(f"Eval {iteration_count[0]}: Residual = {residual:.6f}")
        return residual

    # Run optimization
    if use_global_search:
        print("\n=== STAGE 1: Global Search (Differential Evolution) ===")
        iteration_count[0] = 0
        
        res_global = differential_evolution(
            xrd_objective,
            bounds=bounds,
            strategy='best1bin',  # Good for physics problems
            maxiter=de_maxiter,
            popsize=de_popsize,
            tol=0.1,  # Coarse tolerance for global search
            polish=False,  # Don't use local polish yet
            init='latinhypercube',  # Better initial sampling
            atol=0,
            updating='deferred',  # Faster parallelization (if available)
            workers=1,  # Single thread (can be increased)
            disp=False
        )
        
        print(f"\nGlobal search complete: {res_global.nfev} evaluations")
        print(f"Best global residual: {res_global.fun:.6f}")
        
        print("\n=== STAGE 2: Local Refinement (Powell) ===")
        iteration_count[0] = 0
        p0_refined = res_global.x
        
    else:
        print("\n=== Local Optimization (Powell only) ===")
        iteration_count[0] = 0
        p0_refined = p0

    # Local refinement with Powell
    options = {}
    if maxiter is not None:
        options['maxiter'] = maxiter
        
    res = minimize(
        xrd_objective, 
        p0_refined, 
        bounds=bounds, 
        method='Powell', 
        tol=tol, 
        options=options
    )

    print("\n=== Optimization Complete ===")
    print(f"Success: {res.success}")
    print(f"Message: {res.message}")
    print(f"Final Residual: {res.fun:.6f}")
    
    if use_global_search:
        print(f"Total evaluations: {res_global.nfev + res.nfev}")
    
    return res


def find_peak(ttw_object, x_range):
    """
    Find the maximum peak in a given range from a ttw object.
    Returns (peak_pos, peak_height).
    """
    # Extract data from the first measurement
    df = ttw_object.ttw_df[0]
    x = df['Angle'].values
    y = df['Intensity'].values

    mask = (x >= x_range[0]) & (x <= x_range[1])
    x_subset = x[mask]
    y_subset = y[mask]
    
    if len(x_subset) == 0:
        raise ValueError(f"No data points found in range {x_range}")
        
    max_idx = np.argmax(y_subset)
    return x_subset[max_idx], y_subset[max_idx]


def run_multi_peak_xrd_fitting(ttw_object, peak_labels, peak_intensities, wavelength,
                                 fit_ranges, nx_fits, ny_fits,
                                 use_log_residual=True,
                                 p0_list=None, bounds_list=None,
                                 use_global_search=True,
                                 de_maxiter=50, de_popsize=10,
                                 tol=0.01, maxiter=None,
                                 background_intensity=None,
                                 normalization_peak_intensity=None):
    """
    Multi-peak XRD fitting function. Fits multiple overlapping peaks simultaneously.
    
    Parameters
    ----------
    ttw_object : ttw object
        Experimental data containing Angle and Intensity.
    peak_labels : list of str
        Labels for each peak (e.g., ['BSO', 'SSO']).
    peak_intensities : list of float
        Peak intensities for each film peak (for scaling).
    wavelength : float
        X-ray wavelength in Angstroms.
    fit_ranges : list of tuple
        Fitting range for each peak: [(min1, max1), (min2, max2), ...].
    nx_fits, ny_fits : list of int
        Grid sizes for each peak simulation.
    use_log_residual : bool, default=True
        If True, minimize in log-space (better for Laue oscillations).
    p0_list : list of array-like, optional
        Initial guesses for each peak's parameters (9 params each).
    bounds_list : list of list of tuple, optional
        Bounds for each peak's parameters.
    use_global_search : bool, default=True
        Use Differential Evolution before Powell refinement.
    de_maxiter : int, default=50
        Max iterations for DE.
    de_popsize : int, default=10
        Population size multiplier for DE.
    tol : float, default=0.01
        Tolerance for Powell minimization.
    maxiter : int, optional
        Max iterations for Powell.
    background_intensity : array-like, optional
        Pre-calculated background signal (e.g., from STO substrate).
    normalization_peak_intensity : float, optional
        Peak intensity for normalization (typically STO peak). If None, uses first peak intensity.
    
    Returns
    -------
    results : dict
        Contains optimized parameters for each peak and optimization details.
    """
    
    n_peaks = len(peak_labels)
    
    if normalization_peak_intensity is None:
        normalization_peak_intensity = peak_intensities[0]
    
    print(f"\n{'='*70}")
    print(f"Multi-Peak Fitting: {n_peaks} peaks ({', '.join(peak_labels)})")
    print(f"{'='*70}")
    
    # Calculate overall fitting range
    all_ranges = [r for r in fit_ranges]
    overall_range = (min([r[0] for r in all_ranges]), max([r[1] for r in all_ranges]))
    
    print(f"Overall fitting range: {overall_range}")
    for i, label in enumerate(peak_labels):
        print(f"  {label}: {fit_ranges[i]}, grid {nx_fits[i]}x{ny_fits[i]}")
    
    search_mode = "Global (DE) + Local (Powell)" if use_global_search else "Local (Powell only)"
    print(f"Optimization mode: {search_mode}")
    bg_msg = "Enabled" if background_intensity is not None else "None"
    print(f"  Background: {bg_msg}")
    
    # Extract and normalize experimental data
    df = ttw_object.ttw_df[0]
    x_exp = df['Angle'].values
    y_exp = df['Intensity'].values
    y_exp_norm = y_exp / normalization_peak_intensity
    
    # Prepare data for fitting (use overall range)
    mask_fit = (x_exp >= overall_range[0]) & (x_exp <= overall_range[1])
    x_fit = x_exp[mask_fit]
    y_fit = y_exp_norm[mask_fit]
    
    # Calculate film scale factors
    film_scales = [peak_int / normalization_peak_intensity for peak_int in peak_intensities]
    print(f"\nFilm scale factors: {[f'{s:.4f}' for s in film_scales]}")
    
    # Flatten parameters: [peak1_params (9), peak2_params (9), ...]
    total_params = n_peaks * 9
    
    # Flatten bounds
    if bounds_list is None:
        raise ValueError("bounds_list must be provided for multi-peak fitting")
    bounds_flat = [b for bounds in bounds_list for b in bounds]
    
    # Flatten initial guess
    if p0_list is None:
        p0_flat = np.array([np.mean(b) for b in bounds_flat])
    else:
        p0_flat = np.concatenate([np.array(p0) for p0 in p0_list])
    
    print(f"\nTotal parameters to fit: {total_params} (9 per peak × {n_peaks} peaks)")
    
    # Define objective function
    iteration_count = [0]
    
    # Pre-interpolate background to fit range if provided
    if background_intensity is not None:
        # background_intensity is on full experimental grid (x_exp)
        # We need to interpolate it to the fit range (x_fit)
        bg_interp = np.interp(x_fit, x_exp, background_intensity)
    else:
        bg_interp = None
    
    def multi_peak_objective(params_flat):
        iteration_count[0] += 1
        """
        Objective function for multi-peak fitting.
        Simulates all peaks, sums them, and compares to experiment.
        """
        
        # Split flat parameters into per-peak arrays
        params_per_peak = [params_flat[i*9:(i+1)*9] for i in range(n_peaks)]
        
        # Initialize total simulated intensity
        int_sim_total = np.zeros_like(x_fit)
        
        # Add background if provided (already interpolated to fit range)
        if bg_interp is not None:
            int_sim_total += bg_interp
        
        # Simulate and sum each peak
        for i, (label, params) in enumerate(zip(peak_labels, params_per_peak)):
            c_start, c_end, t_nm, s_var, r_nm, mos, rel_len, att_len, instr = params
            
            # Derived parameters
            nz_fit_layer = int((t_nm * 10) / c_end)
            r_layers = (r_nm * 10) / c_end
            
            # Generate system for this peak
            lat_arr, h_map, t_arr = create_crystal_system(
                nx_fits[i], ny_fits[i], nz_fit_layer,
                c_end,
                [c_start, c_end],
                mode='full',
                gradient_type='exponential',
                spatial_var_percent=s_var,
                roughness_layers=r_layers,
                mosaicity_fwhm=mos,
                relaxation_length=rel_len
            )
            
            # Calculate XRD for this peak
            tt_sim, int_sim = calculate_xrd_vectorized(
                lat_arr, h_map, t_arr,
                two_theta_range=fit_ranges[i],
                wavelength=wavelength,
                attenuation_length_nm=att_len,
                instrumental_broadening=instr,
                n_points=800,
                normalize=True
            )
            
            # Scale this peak's intensity
            int_sim = int_sim * film_scales[i]
            
            # Interpolate to experimental grid and add to total
            int_sim_interp = np.interp(x_fit, tt_sim, int_sim)
            int_sim_total += int_sim_interp
        
        # Calculate residual
        if use_log_residual:
            epsilon = 1e-6
            residual = np.sum((np.log10(y_fit + epsilon) - np.log10(int_sim_total + epsilon))**2)
        else:
            residual = np.sum((y_fit - int_sim_total)**2)
        
        if iteration_count[0] % 50 == 0:
            print(f"  Iteration {iteration_count[0]}: Residual = {residual:.6f}")
        
        return residual
    
    # Run optimization
    print("\n=== Starting Optimization ===")
    
    if use_global_search:
        print(f"Stage 1: Differential Evolution (exploring parameter space)")
        print(f"  maxiter={de_maxiter}, popsize={de_popsize}")
        
        res_global = differential_evolution(
            multi_peak_objective,
            bounds=bounds_flat,
            strategy='best1bin',
            maxiter=de_maxiter,
            popsize=de_popsize,
            tol=0.01,
            mutation=(0.5, 1.0),
            recombination=0.7,
            seed=None,
            init='latinhypercube',
            atol=0,
            updating='immediate',
            workers=1
        )
        
        print(f"  DE Complete: Residual = {res_global.fun:.6f}, Evals = {res_global.nfev}")
        p0_refined = res_global.x
    else:
        p0_refined = p0_flat
    
    # Stage 2: Local refinement with Powell
    print(f"\nStage 2: Powell Local Refinement (polishing solution)")
    print(f"  tol={tol}, maxiter={maxiter if maxiter else 'unlimited'}")
    
    iteration_count[0] = 0  # Reset counter
    
    options = {}
    if maxiter is not None:
        options['maxiter'] = maxiter
    
    res = minimize(
        multi_peak_objective,
        p0_refined,
        bounds=bounds_flat,
        method='Powell',
        tol=tol,
        options=options
    )
    
    print("\n=== Optimization Complete ===")
    print(f"Success: {res.success}")
    print(f"Message: {res.message}")
    print(f"Final Residual: {res.fun:.6f}")
    
    if use_global_search:
        print(f"Total evaluations: {res_global.nfev + res.nfev}")
    
    # Unpack results
    params_optimized = [res.x[i*9:(i+1)*9] for i in range(n_peaks)]
    
    results = {
        'optimization_result': res,
        'params_per_peak': params_optimized,
        'peak_labels': peak_labels,
        'bounds_list': bounds_list,
        'film_scales': film_scales,
        'use_global_search': use_global_search
    }
    
    if use_global_search:
        results['global_result'] = res_global
    
    return results
