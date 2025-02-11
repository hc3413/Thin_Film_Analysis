#import all the libraries needed
from import_dep import *
from RSM_functions import *

class RSM:
    '''Class to store, process and plot XRD RSM data
    '''
    def __init__(self, root_XRD: str, folder_name: str, wavelength: float = 1.5405980):
        # Initialize the class with the root folder and the data folder
        self.root_XRD: Path = Path(root_XRD)
        self.folder_name: Path = Path(folder_name)
        # Construct the full path to the folder from the root directory and folder name
        self.folder_path: Path = self.root_XRD / self.folder_name
           
        # Tuple storing tuples of the raw data for each RSM in case multiple exist in the folder
        self.RSM_df: tuple = None
        self.RSM_np: tuple = None
        self.file_name: tuple = None
        self.plot_string: tuple = None
        
        # tuple to store the processed qx and qz vs intensity values
        self.qxqz_df: tuple = None 
        self.qxqz_np: tuple = None 
        
        # tuple to store the lattice parameters vs intensity values
        self.lat_param_df: tuple = None
        self.lat_param_np: tuple = None
        
        # Wavelength of the X-ray source in Ångstroms is 1.5405980 by default but can be changed in the class initialization
        self.wavelength: float = wavelength
        
        self.points_per_scan: int = None
          
        # Load the data as you initialize the class
        self._load_data()
        self._extract_data()
        self._shape_data()

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
        
    # Functions to convert between reciprocal space coordinates and lattice parameters
    @staticmethod
    def Qx_to_a(Qx):
        return (2 * np.pi) / Qx  # a = 2π/Qx
    @staticmethod
    def a_to_Qx(a):
        return (2 * np.pi) / a
    @staticmethod
    def Qz_to_c(Qz):
        return (6 * np.pi) / Qz  # c = 6π/Qz (from l=3)
    @staticmethod
    def c_to_Qz(c):
        return (6 * np.pi) / c
    #####

    def _extract_data(self):
        if not self.RSM_df:
            print("No data to extract")
            return
        
        lat_param_df = []
        lat_param_np = []
        q_df = []
        q_np = []
        
        for d in self.RSM_df:
                          
            # Step 1: Extract the values from the dataframe
            two_theta = np.radians(d['2Theta position'].values)  # Convert to radians
            omega = np.radians(d['Omega position'].values)  # Convert to radians
            intensity = d['Intensity'].values

            # Step 2: Calculate q_x and q_z for reciprocal space coordinates
            q_x = ((4 * np.pi) / self.wavelength) * (np.sin(two_theta / 2) * np.sin((two_theta / 2) - omega))
            q_z = ((4 * np.pi) / self.wavelength) * (np.sin(two_theta / 2) * np.cos((two_theta / 2) - omega))
            #q_x = ((2 * np.pi) / self.wavelength) * (np.cos(omega) - np.cos(two_theta-omega)) # For the other convention
            #q_z = ((2 * np.pi) / self.wavelength) * (np.sin(omega) + np.sin(two_theta-omega)) # For the other convention
            
            # Step 3: Calculate a and c lattice parameters from q_x and q_z based off this being the (103 peak)
            a = self.Qx_to_a(q_x)
            c = self.Qz_to_c(q_z)
            
            # Step 4: Store the q data in a dataframe and append to the lists
            df_q_params = pd.DataFrame({'qx': q_x, 'qz': q_z, 'Intensity': intensity})
            q_df.append(df_q_params)
            q_np.append(df_q_params.to_numpy())
            
            # Step 5: Store the lattice parameters in a dataframe and append to the lists
            df_lat = pd.DataFrame({'a': a, 'c': c, 'Intensity': intensity})
            lat_param_df.append(df_lat)
            lat_param_np.append(df_lat.to_numpy())
            
            
        self.qxqz_df = q_df
        self.qxqz_np = q_np
        self.lat_param_df = lat_param_df
        self.lat_param_np = lat_param_np      

    def _shape_data(self):
        '''Reshape the data to a 3D array with dimensions (2Theta, Omega, 3) or (Qx, Qz, Intensity)
        Intention is that the data should now be in the format of a meshgrid ready for rendering'''
        if not self.RSM_np:
            print("No data to shape")
            return
        
        RSM_np_reshaped = []
        qxqz_np_reshaped = []
        lat_param_np_reshaped = []
        
        for count, ttw in enumerate(self.RSM_np):
            print(ttw.shape)
            # Extract the unique 2Theta and Omega values ready to reshape the data
            twotheta_unique = np.unique(ttw[:,0]).shape[0]
            omega_unique = np.unique(ttw[:,1]).shape[0]
            print("2Theta unique values:", twotheta_unique, "Omega unique values:", omega_unique)
            
            # Check if the number of unique 2Theta values matches the number of points per scan extracted from the file header
            if twotheta_unique != self.points_per_scan:
                print("Error: Number of 2Theta unique values does not match the number of points per scan extracted from the file.")
            
            # Reshape and transpose the data to a 3D array with dimensions (2Theta, Omega, 3) or (Qx, Qz, Intensity)
            n_reshaped = np.reshape(ttw, (omega_unique,twotheta_unique,  3))
            q_reshaped = np.reshape(self.qxqz_np, (omega_unique,twotheta_unique,  3))
            l_reshaped = np.reshape(self.lat_param_np, (omega_unique,twotheta_unique,  3))
            
            # Transpose the reshaped arrays to move 2Theta/Qx/a to the first dimension and Omega/Qz/c to the second dimension
            n_reshaped_transposed = np.transpose(n_reshaped, (1,0,2))
            q_reshaped_transposed = np.transpose(q_reshaped, (1,0,2))
            l_reshaped_transposed = np.transpose(l_reshaped, (1,0,2))
            
            # Append the reshaped arrays to the list
            RSM_np_reshaped.append(n_reshaped_transposed)
            qxqz_np_reshaped.append(q_reshaped_transposed)
            lat_param_np_reshaped.append(l_reshaped_transposed)
            print(f'{self.plot_string} shape = {n_reshaped_transposed.shape}')

        
        # Update self.RSM_np with the reshaped arrays
        self.RSM_np = RSM_np_reshaped
        self.qxqz_np = qxqz_np_reshaped
        self.lat_param_np = lat_param_np_reshaped
