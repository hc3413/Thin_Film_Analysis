#import all the libraries needed
from import_dep import *
from RSM_functions import *

class ttw:
    '''Class to store, process and plot two theta omega data
    '''
    def __init__(self, root_XRD: str, folder_name: str, wavelength: float = 1.5405980):
        # Initialize the class with the root folder and the data folder
        self.root_XRD: Path = Path(root_XRD)
        self.folder_name: Path = Path(folder_name)
        # Construct the full path to the folder from the root directory and folder name
        self.folder_path: Path = self.root_XRD / self.folder_name
           
        # Tuple storing tuples of the raw data for each RSM in case multiple exist in the folder
        self.ttw_df: tuple = None
        self.ttw_np: tuple = None
        self.file_name: tuple = None
        self.plot_string: tuple = None
        
        # tuple to store the lattice parameters vs intensity values
        self.lat_param_df: tuple = None
        self.lat_param_np: tuple = None
        
        # Wavelength of the X-ray source in Ångstroms is 1.5405980 by default but can be changed in the class initialization
        self.wavelength: float = wavelength
        
        self.points_per_scan: int = None
          
        # Load the data as you initialize the class
        self._load_data()
        if self.ttw_df:  # Ensure data is loaded before extracting
            self._extract_data()


    def _load_data(self) -> None:
        '''Load data from all the csv/txt files in the folder and store it in the class
        '''
        ttw_df = [] # list to store the dataframes
        ttw_np = [] # list to store the numpy arrays
        
        # Get a list of all the files in the folder as an iterator of path objects
        files = [f for f in self.folder_path.iterdir() if f.suffix in ['.csv', '.txt'] and ('2tw' in f.stem or 'GIXRD' in f.stem)]
        # Check if no files were found
        if not files:
            print("Error: No files with .csv or .txt extension were found with 2tw in their name.")
            return
    
        # Sort the files alphabetically
        files.sort()
        
        print(files)
        
        for i, fi in enumerate(files):
            # Load the data from the file
            try:
                with open(fi, 'r', encoding='utf-8', errors='ignore') as file_check:
                    # Read lines until we find the header line
                    header_line = None
                    for line_number, line in enumerate(file_check):
                        # Try to find start of data in the file
                        if "Angle, TimePerStep, Intensity, ESD" in line:
                            header_line = line_number
                            print(header_line)
                            break
                        
                    if header_line is None:
                        print(f"Error: can't find 'Angle, TimePerStep, Intensity, ESD' in file {fi}")
                        continue
                
                # Read the data starting from the header line
                df = pd.read_csv(fi, sep=',', skiprows=header_line, header=0, encoding='utf-8')  # Removed 'errors' argument
                
            except Exception as e:
                print(f"Error with file: {fi}, {e}")
                continue
            
            # Strip any leading/trailing spaces in the column names
            df.columns = df.columns.astype(str).str.strip()
            
            print(df)
            
            ttw_df.append(df)
            ttw_np.append(df.to_numpy())
            
        self.ttw_df = ttw_df
        self.ttw_np = ttw_np
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
        if not self.ttw_df:
            print("No data to extract")
            return
        
        lat_param_df = []
        lat_param_np = []
        
        for d in self.ttw_df:
                          
            # Step 1: Extract the values from the dataframe
            angle = np.radians(d['Angle'].values)  # Convert to radians
            intensity = d['Intensity'].values

            # Step 2: Calculate q_x and q_z for reciprocal space coordinates
            q_z = ((4 * np.pi) / self.wavelength) * np.sin(angle / 2)
            
            # Step 3: Calculate a and c lattice parameters from q_x and q_z
            c = self.Qz_to_c(q_z)
            
            # Step 4: Store the lattice parameters in a dataframe and append to the lists
            df_lat = pd.DataFrame({'c': c, 'Intensity': intensity})
            lat_param_df.append(df_lat)
            lat_param_np.append(df_lat.to_numpy())
            
        self.lat_param_df = lat_param_df
        self.lat_param_np = lat_param_np

    def plot_data(self):
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots()
        
        for df, label_str in zip(self.lat_param_df, self.plot_string):
            x = df['c']
            y = df['Intensity']
            ax.plot(x, y, label=label_str)
        
        ax.set_yscale('log')  # Set y-axis to logarithmic scale
        ax.set_xlabel('Lattice Parameter c')
        ax.set_ylabel('Intensity')
        ax.legend()
        plt.show()
    
    def select_measurements(self, indices):
        """
        Select specific measurements by index and create a new ttw object with only those measurements.
        
        Parameters:
            indices : list or int
                Index or list of indices of the measurements to keep (0-based indexing)
        
        Returns:
            ttw : A new ttw object containing only the selected measurements
        """
        # Convert single index to list
        if isinstance(indices, int):
            indices = [indices]
        
        # Validate indices
        if not self.ttw_df:
            print("No data loaded to select from")
            return self
        
        max_idx = len(self.ttw_df) - 1
        invalid_indices = [i for i in indices if i < 0 or i > max_idx]
        if invalid_indices:
            print(f"Invalid indices: {invalid_indices}. Available indices: 0-{max_idx}")
            return self
        
        # Create a new ttw object (copy of current one)
        new_ttw = ttw.__new__(ttw)  # Create without calling __init__
        
        # Copy basic attributes
        new_ttw.root_XRD = self.root_XRD
        new_ttw.folder_name = self.folder_name
        new_ttw.folder_path = self.folder_path
        new_ttw.wavelength = self.wavelength
        new_ttw.points_per_scan = self.points_per_scan
        
        # Select only the specified indices
        new_ttw.ttw_df = tuple([self.ttw_df[i] for i in indices])
        new_ttw.ttw_np = tuple([self.ttw_np[i] for i in indices])
        new_ttw.file_name = tuple([self.file_name[i] for i in indices])
        new_ttw.plot_string = tuple([self.plot_string[i] for i in indices])
        
        # Select corresponding lattice parameter data if it exists
        if self.lat_param_df:
            new_ttw.lat_param_df = tuple([self.lat_param_df[i] for i in indices])
            new_ttw.lat_param_np = tuple([self.lat_param_np[i] for i in indices])
        else:
            new_ttw.lat_param_df = tuple()
            new_ttw.lat_param_np = tuple()
        
        return new_ttw
    
    def list_measurements(self):
        """
        Print a list of all available measurements with their indices and filenames.
        """
        if not self.ttw_df:
            print("No data loaded")
            return
        
        print(f"Available measurements in {self.folder_name}:")
        print("-" * 50)
        for i, (filename, plot_str) in enumerate(zip(self.file_name, self.plot_string)):
            print(f"Index {i}: {filename}")
            print(f"  Plot string: {plot_str}")
            if self.ttw_df[i] is not None:
                print(f"  Data points: {len(self.ttw_df[i])}")
            print()
