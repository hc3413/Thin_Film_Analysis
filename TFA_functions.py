#import all the libraries needed
from import_dep import *


class XRD_RSM:
    '''Class to store, process and plot XRD RSM data
    '''
    def __init__(self, root_XRD: str, folder_name: str):
        # Initialize the class with the root folder and the data folder
        self.root_XRD: Path = Path(root_XRD)
        self.folder_name: Path = Path(folder_name)
        # Construct the full path to the folder from the root directory and folder name
        self.folder_path: Path = self.root_XRD / self.folder_name
           
        # Tuple storing tuples for each RSM in case multiple exist in the folder
        self.RSM_data: tuple = None
        self.file_name: tuple = None
        self.plot_string: tuple = None
          
        # Load the data as you initialize the class
        self._load_data()

    def _load_data(self) -> None:
        RSM_data = []
        # Get a list of all the files in the folder as an iterator of path objects
        files = list(self.folder.iterdir())
        # Sort the files alphabetically
        files.sort()
        
        print(files)
        
        for i, fi in enumerate(files):
            # Load the data from the file
            try:
                RSM_df = pd.read_csv(self.folder_path / fi, sep='\t', skiprows=6, header=None, comment='N')
                
            except Exception as e:
                print(f"Error with file: {fi}, {e}")
                continue
        
            # Strip any leading/trailing spaces in the column names
            RSM_df.columns = RSM_df.columns.str.strip()
            
            print(RSM_df)
            
            RSM_data.append(RSM_df)
            
        self.RSM_data = RSM_data
        self.file_name = files
        #self.plot_string = tuple()
        


    def analyze_data(self):
        # Logic to analyze the data
        pass

    def plot_2D(self):
        # Logic to plot the data
        pass
    
    def plot_3D(self):
        # Logic to plot the data
        pass
    
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

    
    