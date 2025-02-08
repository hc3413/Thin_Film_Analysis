#import all the libraries needed
from import_dep import *

    
def plot_2D_sep(dat: tuple = None, threshold_filt: float = 0, log_scale: bool = True, x_lim = None, y_lim = None, plot_type: str = 'scatter', ax_unit: str = 'reciprocal'):
    '''
    Plot the 2D data 
    Compare different datasets which are fed into a tuple 'dat'
    '''
    col_maps = ['Purples', 'Greens', 'Oranges', 'Blues', 'Reds', 'Greys']
    # Calculate the number of rows needed based on the number of RSMs
    total_plots = 0
    for count_c, class_obj in enumerate(dat):
        for count_d, d in enumerate(class_obj.qxqz_df):  
            total_plots += 1
    n_rows = total_plots // 2 + total_plots % 2
    
    # Intialize the figure with the number of subplots and thus aspect ratio based on the number of RSMs
    fig = plt.figure(figsize=(15,5.5*n_rows))
    gs = fig.add_gridspec(n_rows, 2)
    
    # Initialize the counter for the subplot to cycle through the colors
    col_counter = 0
    
    for count_c, class_obj in enumerate(dat):
        for count_d, d in enumerate(class_obj.qxqz_df): 
        
            ax = fig.add_subplot(gs[col_counter//2, col_counter % 2])
            
            # Extract the data
            if ax_unit == 'reciprocal':
                
                # Set plot values for qxqz
                x = d['qx']
                y = d['qz']
                z = d['Intensity']
                
                # Set axis labels and title
                ax.set_xlabel(r"$Q_x \,(2\pi/\mathrm{\AA})$", fontsize=12)
                ax.set_ylabel(r"$Q_z \,(2\pi/\mathrm{\AA})$", fontsize=12)
                ax.set_title(r"Reciprocal Space Map ($Q_x$ vs $Q_z$)", fontsize=14)
            
            elif ax_unit == 'lattice':
                
                # Set plot values
                x = class_obj.lat_param_df[count_d]['a']
                y = class_obj.lat_param_df[count_d]['c']
                z = class_obj.lat_param_df[count_d]['Intensity']
                
                # Set axis labels and title
                ax.set_xlabel(r"$a \,(\mathrm{\AA})$", fontsize=12)
                ax.set_ylabel(r"$c \,(\mathrm{\AA})$", fontsize=12)
                ax.set_title(r"Lattice Parameters ($a$ vs $c$)", fontsize=14)
            
            # Filter out non-positive intensities (logarithm can't handle zero or negative values)
            filt_indices = z > threshold_filt 
            
            # Filter the data accordingly
            x_filt = x[filt_indices]
            y_filt = y[filt_indices]
            z_filt = z[filt_indices]
            
            if log_scale:
                plot_intensity = np.log10(z_filt + 1e-99) # Apply a small offset to avoid log(0) for the intensity values 
            else:
                plot_intensity = z_filt
                
    
            
            if plot_type == 'scatter':   
                # Set up the plot
                sc = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap='plasma', s=2.5, edgecolors='none', label=f'{class_obj.plot_string[count_d]}')

                # Add a color bar
                cbar = plt.colorbar(sc, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.ax.tick_params(labelsize=12)
                
            elif plot_type == 'contour':
                
                # Manually specify contour levels to exclude lower intensity values
                contour_levels = np.linspace(1*np.log10(threshold_filt), plot_intensity.max(), 10)

                contour = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap='plasma', s=2.5, edgecolors='none', label=f'{class_obj.plot_string[count_d]}')
                contour_lines = ax.tricontour(x_filt, y_filt, plot_intensity, levels=contour_levels, colors='black', linewidths=0.5)

                # Add a color bar
                cbar = plt.colorbar(contour, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
                cbar.ax.tick_params(labelsize=10)
                
    
                
            elif plot_type == 'mesh':
                from scipy.interpolate import griddata

                # Create a grid for interpolation
                grid_x, grid_y = np.mgrid[min(x_filt):max(x_filt):1000j, min(y_filt):max(y_filt):1000j]

                # Interpolate the data
                grid_intensity = griddata((x_filt, y_filt), plot_intensity, (grid_x, grid_y), method='linear')

                #contour = ax.contourf(grid_x, grid_y, grid_intensity, levels=10, colors='black', linewidths=0.5)
                mesh = ax.pcolormesh(grid_x, grid_y, grid_intensity, cmap='plasma', shading='auto', label=f'{class_obj.plot_string[count_d]}') 
                contour_lines = ax.contour(grid_x, grid_y, grid_intensity, levels=10, colors='black', linewidths=0.5)

                # Add a color bar
                cbar = plt.colorbar(mesh, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
                cbar.ax.tick_params(labelsize=10)
            
            if log_scale:
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
            else:
                cbar.set_label(r"Intensity", fontsize=12)
            

            # Add red lines at x=4.1 and y=4.1
            ax.axvline(x=4.1, color='red', linestyle='--', linewidth=1)
            ax.axhline(y=4.1, color='red', linestyle='--', linewidth=1)
            
            # Add green lines for STO
            ax.axvline(x=3.905, color='green', linestyle='--', linewidth=1)
            ax.axhline(y=3.905, color='green', linestyle='--', linewidth=1)
            
            # Modify tick parameters for better visibility
            ax.tick_params(axis='both', which='major', labelsize=12, length=4, width=0.6)
            ax.tick_params(axis='both', which='minor', labelsize=12, length=2, width=0.4)

            # Set axis limits (optional, adjust to your data)
            ax.set_xlim(x_lim)
            ax.set_ylim(y_lim)

            # Improve layout
            fig.tight_layout()

            # Adjust spacing so that the plotted region aligns with the figure without a colorbar
            fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
            
            # Add a legend
            ax.legend(loc='best', framealpha=0.4)

    # Show the plot
    plt.show()
    return fig
    
    
def plot_2D_same(dat: tuple = None, threshold_filt: float = 0, log_scale: bool = True, x_lim = None, y_lim = None, plot_type: str = 'scatter', ax_unit: str = 'reciprocal'):
    '''
    Plot the 2D data 
    Compare different datasets which are fed into a tuple 'dat'
    '''
    col_maps = ['Purples', 'Greens', 'Oranges', 'Blues', 'Reds', 'Greys']

    
    # Intialize the figure 
    fig = plt.figure(figsize=(15,5.5))
    
    # Initialize the counter for the subplot to cycle through the colors
    col_counter = 0
    
    for count_c, class_obj in enumerate(dat):
        for count_d, d in enumerate(class_obj.qxqz_df): 
        
            ax = fig.add_subplot(111)
            
            # Extract the data
            if ax_unit == 'reciprocal':
                
                # Set plot values for qxqz
                x = d['qx']
                y = d['qz']
                z = d['Intensity']
                
                # Set axis labels and title
                ax.set_xlabel(r"$Q_x \,(2\pi/\mathrm{\AA})$", fontsize=12)
                ax.set_ylabel(r"$Q_z \,(2\pi/\mathrm{\AA})$", fontsize=12)
                ax.set_title(r"Reciprocal Space Map ($Q_x$ vs $Q_z$)", fontsize=14)
            
            elif ax_unit == 'lattice':
                
                # Set plot values
                x = class_obj.lat_param_df[count_d]['a']
                y = class_obj.lat_param_df[count_d]['c']
                z = class_obj.lat_param_df[count_d]['Intensity']
                
                # Set axis labels and title
                ax.set_xlabel(r"$a \,(\mathrm{\AA})$", fontsize=12)
                ax.set_ylabel(r"$c \,(\mathrm{\AA})$", fontsize=12)
                ax.set_title(r"Lattice Parameters ($a$ vs $c$)", fontsize=14)
            
            # Filter out non-positive intensities (logarithm can't handle zero or negative values)
            filt_indices = z > threshold_filt 
            
            # Filter the data accordingly
            x_filt = x[filt_indices]
            y_filt = y[filt_indices]
            z_filt = z[filt_indices]
            
            if log_scale:
                plot_intensity = np.log10(z_filt + 1e-99) # Apply a small offset to avoid log(0) for the intensity values 
            else:
                plot_intensity = z_filt
                
    
            
            if plot_type == 'scatter':   
                # Set up the plot
                sc = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=col_maps[col_counter], s=2.5, edgecolors='none', label=f'{class_obj.plot_string[count_d]}')

                # Add a color bar
                cbar = plt.colorbar(sc, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.ax.tick_params(labelsize=12)
                
            elif plot_type == 'contour':
                
                # Manually specify contour levels to exclude lower intensity values
                contour_levels = np.linspace(1*np.log10(threshold_filt), plot_intensity.max(), 10)

                contour = ax.scatter(x_filt, y_filt, c=plot_intensity, cmap=col_maps[col_counter], s=2.5, edgecolors='none', label=f'{class_obj.plot_string[count_d]}')
                contour_lines = ax.tricontour(x_filt, y_filt, plot_intensity, levels=contour_levels, colors='black', linewidths=0.5)

                # Add a color bar
                cbar = plt.colorbar(contour, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
                cbar.ax.tick_params(labelsize=10)
                
    
                
            elif plot_type == 'mesh':
                from scipy.interpolate import griddata

                # Create a grid for interpolation
                grid_x, grid_y = np.mgrid[min(x_filt):max(x_filt):1000j, min(y_filt):max(y_filt):1000j]

                # Interpolate the data
                grid_intensity = griddata((x_filt, y_filt), plot_intensity, (grid_x, grid_y), method='linear')

                #contour = ax.contourf(grid_x, grid_y, grid_intensity, levels=10, colors='black', linewidths=0.5)
                mesh = ax.pcolormesh(grid_x, grid_y, grid_intensity, cmap=col_maps[col_counter], shading='auto', label=f'{class_obj.plot_string[count_d]}')
                contour_lines = ax.contour(grid_x, grid_y, grid_intensity, levels=10, colors='black', linewidths=0.5)

                # Add a color bar
                cbar = plt.colorbar(mesh, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
                cbar.ax.tick_params(labelsize=10)
            
                # Increment the color counter
                col_counter +=1 
                
    if log_scale:
        cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
    else:
        cbar.set_label(r"Intensity", fontsize=12)
    

    # Add red lines at x=4.1 and y=4.1
    ax.axvline(x=4.1, color='red', linestyle='--', linewidth=1)
    ax.axhline(y=4.1, color='red', linestyle='--', linewidth=1)
    
    # Add green lines for STO
    ax.axvline(x=3.905, color='green', linestyle='--', linewidth=1)
    ax.axhline(y=3.905, color='green', linestyle='--', linewidth=1)
    
    # Modify tick parameters for better visibility
    ax.tick_params(axis='both', which='major', labelsize=12, length=4, width=0.6)
    ax.tick_params(axis='both', which='minor', labelsize=12, length=2, width=0.4)

    # Set axis limits (optional, adjust to your data)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    # Improve layout
    fig.tight_layout()

    # Adjust spacing so that the plotted region aligns with the figure without a colorbar
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    
    # Add a legend
    ax.legend(loc='best', framealpha=0.4)

            
    # Show the plot
    plt.show()
    return fig
# class XRD_2theta:
#     def __init__(self, file_path: str):
#         self.file_path: str = file_path
#         self.data_import_np: np.ndarray = None  # raw data from the file
#         self.data_import_df: pd.DataFrame = None  # data in a pandas dataframe
#         self.data: np.ndarray = self._load_data()

#     def _load_data(self) -> np.ndarray:
#         # Logic to load data from the file
#         pass

#     def analyze_data(self):
#         # Logic to analyze the data
#         pass

#     def plot_data(self):
#         # Logic to plot the data
#         pass

# class XRD_XRR:
#     def __init__(self, file_path: str):
#         self.file_path: str = file_path
#         self.data_import_np: np.ndarray = None  # raw data from the file
#         self.data_import_df: pd.DataFrame = None  # data in a pandas dataframe
#         self.data: np.ndarray = self._load_data()

#     def _load_data(self) -> np.ndarray:
#         # Logic to load data from the file
#         pass

#     def analyze_data(self):
#         # Logic to analyze the data
#         pass

#     def plot_data(self):
#         # Logic to plot the data
#         pass

    
    