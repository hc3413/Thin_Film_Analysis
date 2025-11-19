import numpy as np
import time
from scipy.signal import convolve
from scipy.optimize import minimize

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
                             n_points=3000, chunk_size=2000):
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
    if intensity_total.max() > 0:
        intensity_total = intensity_total / intensity_total.max()

    elapsed = time.time() - start_time
    print(f"XRD calc finished in {elapsed:.3f}s (columns={n_columns})")
    return two_theta, intensity_total


def run_xrd_fitting(ttw_object, peak_intensity, wavelength, 
                    fit_range=(42.0, 45.5), 
                    nx_fit=8, ny_fit=8, 
                    use_log_residual=True,
                    p0=None, bounds=None,
                    tol=0.01, maxiter=None):
    """
    Run the XRD fitting optimization.
    tol: Tolerance for termination (lower = more precise, slower). Default 0.01.
    maxiter: Maximum number of iterations.
    """
    
    print(f"Setting up optimization in range {fit_range} with grid {nx_fit}x{ny_fit}...")
    print(f"Optimization settings: tol={tol}, maxiter={maxiter}")

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
    def xrd_objective(params):
        """
        Objective function for optimization.
        Returns the sum of squared errors between simulation and experiment.
        """
        # Unpack parameters
        # c_start: Lattice parameter at interface
        # c_end:   Lattice parameter at surface (relaxed/perfect)
        c_start, c_end, t_nm, s_var, r_nm, mos, rel_len, att_len, instr = params
        
        # Derived parameters
        # Use c_end (relaxed) as the reference for layer count
        nz_fit_layer = int((t_nm * 10) / c_end)
        r_layers = (r_nm * 10) / c_end
        
        # Generate System
        lat_arr, h_map, t_arr = create_crystal_system(
            nx_fit, ny_fit, nz_fit_layer,
            c_end,                  # lattice_perfect (reference for spatial var)
            [c_start, c_end],       # lattice_gradient [interface, surface]
            mode='full',
            gradient_type='exponential',
            spatial_var_percent=s_var,
            roughness_layers=r_layers,
            mosaicity_fwhm=mos,
            relaxation_length=rel_len
        )
        
        # Calculate XRD
        # Use a moderate resolution for fitting speed (1000 points)
        tt_sim, int_sim = calculate_xrd_vectorized(
            lat_arr, h_map, t_arr,
            two_theta_range=fit_range, # Use the fitting range!
            wavelength=wavelength,
            attenuation_length_nm=att_len,
            instrumental_broadening=instr,
            n_points=1000 
        )
        
        # Interpolate Simulation to Experimental X-axis to calculate residual
        int_sim_interp = np.interp(x_fit, tt_sim, int_sim)
        
        # Calculate Residual (Least Squares)
        if use_log_residual:
            # Log scale fitting: Good for Laue oscillations
            epsilon = 1e-6 # Avoid log(0)
            residual = np.sum((np.log10(y_fit + epsilon) - np.log10(int_sim_interp + epsilon))**2)
        else:
            # Linear scale fitting: Good for main Bragg peak position/width
            residual = np.sum((y_fit - int_sim_interp)**2)
        
        return residual

    print(f"Initial parameters: {p0}")
    print("Optimizing... (This may take 1-2 minutes)")

    # We use 'Powell' method because it handles the non-smooth nature of the 
    # integer thickness (nz) better than gradient-based methods like BFGS.
    options = {}
    if maxiter is not None:
        options['maxiter'] = maxiter
        
    res = minimize(xrd_objective, p0, bounds=bounds, method='Powell', tol=tol, options=options)

    print("\nOptimization Result:")
    print(f"Success: {res.success}")
    print(f"Message: {res.message}")
    print(f"Final Residual: {res.fun:.6f}")
    
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
