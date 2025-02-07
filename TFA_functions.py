#import all the libraries needed
from import_dep import *


class XRD_RSM:
    '''Class to store, process and plot XRD RSM data
    '''
    def __init__(self, root_XRD: str, folder_name: str, wavelength: float = 1.5405980):
        # Initialize the class with the root folder and the data folder
        self.root_XRD: Path = Path(root_XRD)
        self.folder_name: Path = Path(folder_name)
        # Construct the full path to the folder from the root directory and folder name
        self.folder_path: Path = self.root_XRD / self.folder_name
           
        # Tuple storing tuples for each RSM in case multiple exist in the folder
        self.RSM_df: tuple = None
        self.RSM_np: tuple = None
        self.file_name: tuple = None
        self.plot_string: tuple = None
        self.qxqz_df: tuple = None # tuple to store the processed qx and qz vs intensity values
        
        # Wavelength of the X-ray source in Ångstroms is 1.5405980 by default but can be changed in the class initialization
        self.wavelength: float = wavelength
        
        self.points_per_scan: int = None
          
        # Load the data as you initialize the class
        self._load_data()
        self._shape_data()
        self._extract_data()

    def _load_data(self) -> None:
        '''Load data from all the csv/txt files in the folder and store it in the class
        '''
        RSM_df = [] # list to store the dataframes
        RSM_np = [] # list to store the numpy arrays
        
        # Get a list of all the files in the folder as an iterator of path objects
        files = [f for f in self.folder_path.iterdir() if f.suffix in ['.csv', '.txt'] and 'RSM' in f.stem]
        # Check if no files were found
        if not files:
            print("Error: No files with .csv or .txt extension containing 'RSM' in the filename were found.")
            return
    
        # Sort the files alphabetically
        files.sort()
        
        print(files)
        
        for i, fi in enumerate(files):
            # Load the data from the file
            try:
                with open(fi, 'r') as file_check:
                    # Read lines until we find the header line
                    header_line = None
                    for line_number, line in enumerate(file_check):
                        # Try to find start of data in the file
                        if "2Theta position" in line:
                            header_line = line_number
                            print(header_line)
                        
                        # Check if the number of points per scan is specified in the file
                        if "No. of points per scan" in line:
                            self.points_per_scan = int(line.split(',')[-1].strip())
                            print(f"Points per scan: {self.points_per_scan}")
                        
                        if header_line is not None and self.points_per_scan is not None:
                            break
                        
                    if header_line is None:
                        print(f"Error: can't find '2Theta position' in file {fi}")
                        continue
                
                # Read the data starting from the header line
                df = pd.read_csv(fi, sep=',', skiprows=header_line, header=0)
                
            except Exception as e:
                print(f"Error with file: {fi}, {e}")
                continue
            
            # Strip any leading/trailing spaces in the column names
            df.columns = df.columns.astype(str).str.strip()
            
            print(df)
            
            RSM_df.append(df)
            RSM_np.append(df.to_numpy())
            
        self.RSM_df = RSM_df
        self.RSM_np = RSM_np
        self.file_name = [str(p) for p in files] # store the file names as strings not path objects
        self.plot_string = [p.stem for p in files ]
        
    def _shape_data(self):
        if not self.RSM_np:
            print("No data to shape")
            return
        
        RSM_np_reshaped = []
        
        for n in self.RSM_np:
            print(n.shape)
            # Extract the unique 2Theta and Omega values ready to reshape the data
            twotheta_unique = np.unique(n[:,0]).shape[0]
            omega_unique = np.unique(n[:,1]).shape[0]
            print("2Theta unique values:", twotheta_unique, "Omega unique values:", omega_unique)
            
            # Check if the number of unique 2Theta values matches the number of points per scan extracted from the file header
            if twotheta_unique != self.points_per_scan:
                print("Error: Number of 2Theta unique values does not match the number of points per scan extracted from the file.")
            
            # Reshape and transpose the data to a 3D array with dimensions (2Theta, Omega, 3)
            n_reshaped = np.reshape(n, (omega_unique,twotheta_unique,  3))
            n_reshaped_transposed = np.transpose(n_reshaped, (1,0,2))
            RSM_np_reshaped.append(n_reshaped_transposed)

        
        # Update self.RSM_np with the reshaped arrays
        self.RSM_np = RSM_np_reshaped


    def _extract_data(self):
        if not self.RSM_df:
            print("No data to extract")
            return
        
        qxqz_df = []
        
        for d in self.RSM_df:
                          
            # Extract columns
            theta = np.radians(d['2Theta position'].values)  # Convert to radians
            omega = np.radians(d['Omega position'].values)  # Convert to radians
            intensity = d['Intensity'].values

            # Step 2: Calculate q_x and q_z for reciprocal space coordinates
            q_x = ((4 * np.pi) / self.wavelength) * (np.sin(theta / 2) * np.sin((theta / 2) - omega))
            q_z = ((4 * np.pi) / self.wavelength) * (np.sin(theta / 2) * np.cos((theta / 2) - omega))
      
            qxqz_df.append(pd.DataFrame({'qx': q_x, 'qz': q_z, 'Intensity': intensity}))
            
        self.qxqz_df = qxqz_df      



    def plot_2D(self, threshold_filt: float = 0, log_scale: bool = True, x_lim = None, y_lim = None, plot_type: str = 'scatter'):
        '''
        Plot the 2D data in a 3D scatter plot
        threshold_filt lets you filter out low intensity values
        '''
        # Calculate the number of rows needed based on the number of RSMs
        n_rows = len(self.qxqz_df) // 2 + len(self.qxqz_df) % 2
        
        # Intialize the figure with the number of subplots and thus aspect ratio based on the number of RSMs
        fig = plt.figure(figsize=(15,5.5*n_rows))
        gs = fig.add_gridspec(n_rows, 2)
        
        for count, q in enumerate(self.qxqz_df):
            
            ax = fig.add_subplot(gs[count//2, count % 2])
            
            # Extract the data
            qx = q['qx']
            qz = q['qz']
            intensity = q['Intensity']
            
            # Filter out non-positive intensities (logarithm can't handle zero or negative values)
            filt_indices = intensity > threshold_filt 
            
            # Filter the data accordingly
            q_x_filt = qx[filt_indices]
            q_z_filt = qz[filt_indices]
            intensity_filt = intensity[filt_indices]
            
            if log_scale:
                plot_intensity = np.log10(intensity_filt + 1e-99) # Apply a small offset to avoid log(0) for the intensity values 
            else:
                plot_intensity = intensity_filt
            
            if plot_type == 'scatter':   
                # Set up the plot
                sc = ax.scatter(q_x_filt, q_z_filt, c=plot_intensity, cmap='viridis', s=7, edgecolors='none')

                # Add a color bar
                cbar = plt.colorbar(sc, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.ax.tick_params(labelsize=12)
                
            elif plot_type == 'contour':
                contour = ax.tricontourf(q_x_filt, q_z_filt, plot_intensity, levels=1000, cmap='plasma')
                contour_lines = ax.tricontour(q_x_filt, q_z_filt, plot_intensity, levels=10, colors='black', linewidths=0.5)

                # Add a color bar
                cbar = plt.colorbar(contour, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
                cbar.ax.tick_params(labelsize=10)
            
            if log_scale:
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
            else:
                cbar.set_label(r"Intensity", fontsize=12)
                
            

            # Set axis labels and title
            ax.set_xlabel(r"$Q_x \,(2\pi/\mathrm{\AA})$", fontsize=12)
            ax.set_ylabel(r"$Q_z \,(2\pi/\mathrm{\AA})$", fontsize=12)
            ax.set_title(r"Reciprocal Space Map ($Q_x$ vs $Q_z$)", fontsize=14)

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

        # Show the plot
        plt.show()
        return fig
    
    
    def plot_3D(self, threshold_filt: float = 0, log_scale: bool = True, x_lim = None, y_lim = None):
        '''Plot the 3D data in a 3D surface plot
        '''
        # Calculate the number of rows needed based on the number of RSMs
        n_rows = len(self.qxqz_df) // 2 + len(self.qxqz_df) % 2
        
        # Intialize the figure with the number of subplots and thus aspect ratio based on the number of RSMs
        fig = plt.figure(figsize=(15,5.5*n_rows))
        gs = fig.add_gridspec(n_rows, 2)
        
        for count, q in enumerate(self.qxqz_df):
            
            ax = fig.add_subplot(gs[count//2, count % 2])
            
            # Extract the data
            qx = q['qx']
            qz = q['qz']
            intensity = q['Intensity']
            
            # Filter out non-positive intensities (logarithm can't handle zero or negative values)
            filt_indices = intensity > threshold_filt 
            
            # Filter the data accordingly
            q_x_filt = qx[filt_indices]
            q_z_filt = qz[filt_indices]
            intensity_filt = intensity[filt_indices]
            
            if log_scale:
                plot_intensity = np.log10(intensity_filt + 1e-99) # Apply a small offset to avoid log(0) for the intensity values 
            else:
                plot_intensity = intensity_filt
            

            # Set up the 3D plot
            ax = fig.add_subplot(gs[count//2, count % 2], projection='3d')
            ax.plot_trisurf(q_x_filt, q_z_filt, plot_intensity, cmap='viridis', edgecolor='none')

            # Add a color bar
            sc = ax.scatter(q_x_filt, q_z_filt, c=plot_intensity, cmap='viridis')
            cbar = plt.colorbar(sc, ax=ax, pad=0.02, format=ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
            cbar.ax.tick_params(labelsize=12)
                
 
            
            if log_scale:
                cbar.set_label(r"$\log_{10}(\mathrm{Intensity})$", fontsize=12)
            else:
                cbar.set_label(r"Intensity", fontsize=12)
                
            

            # Set axis labels and title
            ax.set_xlabel(r"$Q_x \,(2\pi/\mathrm{\AA})$", fontsize=12)
            ax.set_ylabel(r"$Q_z \,(2\pi/\mathrm{\AA})$", fontsize=12)
            ax.set_title(r"Reciprocal Space Map ($Q_x$ vs $Q_z$)", fontsize=14)

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

        # Show the plot
        plt.show()
        return fig

    
def plot_3D_comparison(data1, data2):
    # Logic to plot the data
    pass
    
    
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

    
    