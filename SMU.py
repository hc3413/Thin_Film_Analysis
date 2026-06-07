#import all the libraries needed
from import_dep import *
import re

class SMU:
    '''Class to store, process and plot SMU (Source Measure Unit) data
    '''
    def __init__(self, root_SMU: str, folder_name: Optional[str] = None):
        # Initialize the class with the root folder 
        self.root_SMU: Path = Path(root_SMU)
        
        # If folder_name is provided, use it; otherwise use root_SMU directly
        if folder_name:
            self.folder_name: Path = Path(folder_name)
            # Construct the full path to the folder from the root directory and folder name
            self.folder_path: Path = self.root_SMU / self.folder_name
        else:
            self.folder_path: Path = self.root_SMU
            self.folder_name: Path = Path('')
           
        # Lists to store data for each measurement type
        self.cv_data: list = []
        self.iv_data: list = []
        self.iv3T_data: list = []
        self.cv_run_nums: list = []
        self.iv_run_nums: list = []
        self.iv3T_run_nums: list = []
        self.cv_file_names: list = []
        self.iv_file_names: list = []
        self.iv3T_file_names: list = []
        self.cv_plot_strings: list = []
        self.iv_plot_strings: list = []
        self.iv3T_plot_strings: list = []
        
        # Load the data as you initialize the class
        self._load_data()

    def _load_data(self) -> None:
        '''Load data from all the xlsx files in the folder and categorize them
        '''
        # Get a list of all the xlsx files in the folder
        files = [f for f in self.folder_path.iterdir() if f.suffix == '.xlsx' and not f.name.startswith('~$')]
        
        # Check if no files were found
        if not files:
            print("Error: No .xlsx files were found.")
            return
    
        # Sort the files alphabetically
        files.sort()
        
        print(f"Found {len(files)} Excel files")
        
        for fi in files:
            try:
                # Extract run number from filename
                run_match = re.search(r'Run(\d+)', fi.name)
                if run_match:
                    run_num = int(run_match.group(1))
                else:
                    run_num = None
                    print(f"Warning: Could not extract run number from {fi.name}")
                
                # Check file type and process accordingly
                if 'cv-cap' in fi.name.lower() :
                    self._process_cv_file(fi, run_num)
                elif 'res2t' in fi.name.lower() or 'iv_sweep' in fi.name.lower():
                    self._process_iv_file(fi, run_num)
                elif 'fet' in fi.name.lower() or 'id_vs_vg' in fi.name.lower():
                    self._process_iv3T_file(fi, run_num)
                else:
                    print(f"Warning: Unknown file type for {fi.name}")
                    
            except Exception as e:
                print(f"Error processing file {fi.name}: {e}")
                continue
        
        print(f"Loaded {len(self.cv_data)} CV, {len(self.iv_data)} IV, {len(self.iv3T_data)} IV3T files")

    def _process_cv_file(self, file_path: Path, run_num: Optional[int]) -> None:
        '''Process CV (capacitance-voltage) files'''
        try:
            # Read the Excel file
            df = pd.read_excel(file_path, header=0)  # Use first row as header
            
            # Extract the required columns based on the specification
            # Column 1 = Cp, Column 2 = Gp, Column 4 = V_DC, Column 5 = frequency
            cv_df = pd.DataFrame({
                'Cp': df.iloc[:, 0],          # Column 1 (0-indexed)
                'Gp': df.iloc[:, 1],          # Column 2 (0-indexed)  
                'V_DC': df.iloc[:, 2],        # Column 3 (0-indexed) - this is DCV_AB
                'Frequency_Hz': df.iloc[:, 3] # Column 4 (0-indexed) - this is F_AB
            })
            
            # Remove the header row if it contains text
            cv_df = cv_df[pd.to_numeric(cv_df['Cp'], errors='coerce').notna()]
            
            # Convert to numeric
            for col in cv_df.columns:
                cv_df[col] = pd.to_numeric(cv_df[col], errors='coerce')
            
            # Store the data
            self.cv_data.append(cv_df)
            self.cv_run_nums.append(run_num)
            self.cv_file_names.append(str(file_path))
            # Clean up plot string to avoid LaTeX issues
            clean_name = self._clean_plot_string(file_path.stem)
            self.cv_plot_strings.append(f"Run{run_num}" if run_num is not None else clean_name)
            
        except Exception as e:
            print(f"Error processing CV file {file_path.name}: {e}")

    def _process_iv_file(self, file_path: Path, run_num: Optional[int]) -> None:
        '''Process IV (current-voltage) files'''
        try:
            # Read the Excel file
            df = pd.read_excel(file_path, header=0)  # Use first row as header
            
            # Extract the required columns based on the specification
            # Column 1 = I (A), Column 2 = V (V)
            iv_df = pd.DataFrame({
                'I_A': df.iloc[:, 0],    # Column 1 (0-indexed) - current in Amperes
                'V_V': df.iloc[:, 1]     # Column 2 (0-indexed) - voltage in Volts
            })
            
            # Remove the header row if it contains text
            iv_df = iv_df[pd.to_numeric(iv_df['I_A'], errors='coerce').notna()]
            
            # Convert to numeric
            for col in iv_df.columns:
                iv_df[col] = pd.to_numeric(iv_df[col], errors='coerce')
            
            # Calculate resistance R = V/I (avoiding division by zero)
            iv_df['R_Ohm'] = iv_df['V_V'] / iv_df['I_A']
            iv_df['R_Ohm'] = iv_df['R_Ohm'].replace([np.inf, -np.inf], np.nan)
            
            # Store the data
            self.iv_data.append(iv_df)
            self.iv_run_nums.append(run_num)
            self.iv_file_names.append(str(file_path))
            # Clean up plot string to avoid LaTeX issues
            clean_name = self._clean_plot_string(file_path.stem)
            self.iv_plot_strings.append(f"Run{run_num}" if run_num is not None else clean_name)
            
        except Exception as e:
            print(f"Error processing IV file {file_path.name}: {e}")
            
            
    def _process_iv3T_file(self, file_path: Path, run_num: Optional[int]) -> None:
        '''Process IV 3-terminal (current-voltage) files for FETs with gate, drain, source'''
        try:
            # Read the Excel file
            df = pd.read_excel(file_path, header=0)  # Use first row as header
            
            # Extract the required columns based on the specification
            # Column 1 = I (A), Column 2 = V (V)
            iv3T_df = pd.DataFrame({
                'D_I_A': df.iloc[:, 0],    # Column 1 (0-indexed) - drain current in Amperes
                'D_V_V': df.iloc[:, 1],     # Column 2 (0-indexed) - drain voltage in Volts
                'G_I_A': df.iloc[:, 2],    # Column 3 (0-indexed) - gate current in Amperes
                'G_V_V': df.iloc[:, 3],     # Column 4 (0-indexed) - gate voltage in Volts
                'S_I_A': df.iloc[:, 4],    # Column 5 (0-indexed) - source current in Amperes
                'S_V_V': df.iloc[:, 5]     # Column 6 (0-indexed) - source voltage in Volts
            })
            
            # Remove the header row if it contains text
            iv3T_df = iv3T_df[pd.to_numeric(iv3T_df['D_I_A'], errors='coerce').notna()]
            
            # Convert to numeric
            for col in iv3T_df.columns:
                iv3T_df[col] = pd.to_numeric(iv3T_df[col], errors='coerce')
            
            # Calculate channel resistance R = V/I (avoiding division by zero)
            iv3T_df['D_R_Ohm'] = iv3T_df['D_V_V'] / iv3T_df['D_I_A']
            iv3T_df['D_R_Ohm'] = iv3T_df['D_R_Ohm'].replace([np.inf, -np.inf], np.nan)
            
            # Calculate gate resistance Rg = Vg/Ig (avoiding division by zero)
            iv3T_df['G_R_Ohm'] = iv3T_df['G_V_V'] / iv3T_df['G_I_A']
            iv3T_df['G_R_Ohm'] = iv3T_df['G_R_Ohm'].replace([np.inf, -np.inf], np.nan)
            
            # Store the data
            self.iv3T_data.append(iv3T_df)
            self.iv3T_run_nums.append(run_num)
            self.iv3T_file_names.append(str(file_path))
            # Clean up plot string to avoid LaTeX issues
            clean_name = self._clean_plot_string(file_path.stem)
            self.iv3T_plot_strings.append(f"Run{run_num}" if run_num is not None else clean_name)
            
        except Exception as e:
            print(f"Error processing IV3T file {file_path.name}: {e}")

    def _clean_plot_string(self, text: str) -> str:
        """Clean text string to avoid LaTeX rendering issues"""
        # Replace problematic characters for LaTeX
        text = text.replace('#', 'num')
        text = text.replace('@', 'at')
        text = text.replace('&', 'and')
        text = text.replace('%', 'pct')
        text = text.replace('$', 'dollar')
        text = text.replace('_', ' ')
        text = text.replace('~', '-')
        return text

    def get_cv_data(self, run_nums: Optional[list] = None) -> tuple:
        '''Get CV data for specified run numbers or all if none specified'''
        if run_nums is None:
            return self.cv_data, self.cv_plot_strings
        
        selected_data = []
        selected_strings = []
        for i, run_num in enumerate(self.cv_run_nums):
            if run_num in run_nums:
                selected_data.append(self.cv_data[i])
                selected_strings.append(self.cv_plot_strings[i])
        
        return selected_data, selected_strings

    def get_iv_data(self, run_nums: Optional[list] = None) -> tuple:
        '''Get IV data for specified run numbers or all if none specified'''
        if run_nums is None:
            return self.iv_data, self.iv_plot_strings
        
        selected_data = []
        selected_strings = []
        for i, run_num in enumerate(self.iv_run_nums):
            if run_num in run_nums:
                selected_data.append(self.iv_data[i])
                selected_strings.append(self.iv_plot_strings[i])
        
        return selected_data, selected_strings
    
    def get_iv3T_data(self, run_nums: Optional[list] = None) -> tuple:
        '''Get IV3T data for specified run numbers or all if none specified'''
        if run_nums is None:
            return self.iv3T_data, self.iv3T_plot_strings
        
        selected_data = []
        selected_strings = []
        for i, run_num in enumerate(self.iv3T_run_nums):
            if run_num in run_nums:
                selected_data.append(self.iv3T_data[i])
                selected_strings.append(self.iv3T_plot_strings[i])
        
        return selected_data, selected_strings

    def __repr__(self):
        return f"SMU(folder_path='{self.folder_path}', CV_files={len(self.cv_data)}, IV_files={len(self.iv_data)}, IV3T_files={len(self.iv3T_data)})"
