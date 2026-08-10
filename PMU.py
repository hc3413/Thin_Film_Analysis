"""Import classes for pulsed / ferroelectric measurements (PMU).

This module is the pulsed-measurement counterpart to ``SMU.py``.  Where ``SMU``
handles DC source-measure data (IV, CV, FET transfer curves), ``PMU`` handles
measurements made with a pulse-measure unit or a dedicated ferroelectric
tester: PUND, DHM (dynamic hysteresis) and fatigue.

Structure follows the convention used in the sibling ``IS_Analysis`` project
(``IS_Import.py``): a dataclass holds one measurement, the dataclasses live in a
dict inside a per-sample container object, and there is one concrete container
class per (instrument x measurement type):

    PUND      Keithley_PUND, Keithley_Haoran_PUND, AixACCT_PUND
    DHM       AixACCT_DHM
    Fatigue   AixACCT_Fatigue
    Custom    Keithley_Custom

Each container scans one sample folder and picks up only the files belonging to
its own measurement type, so a folder holding a PundFile, a HysterFile and a
FatigueFile is opened as three separate objects.

Constructing an object both reads the files **and** processes them: each class
calls ``_load_data()`` then ``_transform_data()``, so by the time it returns,
every run already carries its processed ``.pund`` / ``.dhm`` loop and its
``.extracted_params``.  The plotting functions in ``PMU_functions.py`` therefore
only read stored data, and the same numbers back the plot, the summary table and
the CSV export.  Pass ``transform=False`` to import without processing.

The transform algorithms themselves live in ``PMU_functions.py``; this module
calls into them.
"""

# import all the libraries needed
from import_dep import *
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import copy
import re


# =============================================================================
#  Dataclasses - one per measurement
# =============================================================================

def _fmt_time(seconds: float) -> str:
    """Human-readable time with an SI prefix, e.g. 0.001 -> '1 ms'.

    Plain text rather than LaTeX: these go into legend labels that pass through
    :func:`_clean_plot_string`, which strips '$'.
    """
    for scale, unit in ((1.0, 's'), (1e-3, 'ms'), (1e-6, 'us'), (1e-9, 'ns')):
        if abs(seconds) >= scale:
            return f"{seconds / scale:.3g} {unit}"
    return f"{seconds:.3g} s"


def _clean_plot_string(text: str) -> str:
    """Strip characters that break LaTeX rendering in matplotlib labels."""
    for bad, good in (('#', 'num'), ('@', 'at'), ('&', 'and'), ('%', 'pct'),
                      ('$', 'dollar'), ('_', ' '), ('~', '-'), ('^', ' ')):
        text = text.replace(bad, good)
    return text


@dataclass
class PMUdata:
    """A single pulsed measurement: one PUND train, one DHM loop set, or one
    custom waveform capture.

    ``raw_data`` always holds the as-measured trace.  ``pund`` / ``dhm`` are
    filled in later by the transform functions in ``PMU_functions.py`` and hold
    the processed loop; ``extracted_params`` holds the scalars derived from it.
    ``instrument_params`` keeps whatever the instrument itself reported, so our
    numbers and the instrument's can be compared without either overwriting the
    other.
    """

    # --- identity ---
    file_name: str
    run_number: int
    folder_path: Path
    instrument: str                 # 'keithley' | 'keithley_haoran' | 'aixacct'
    meas_type: str                  # 'PUND' | 'DHM' | 'CUSTOM'

    # --- as measured ---
    raw_data: pd.DataFrame
    metadata: dict = field(default_factory=dict)

    # --- device geometry ---
    electrode_area: Optional[float] = None      # m^2
    film_thickness_nm: Optional[float] = None   # nm

    # --- transformed (filled by PMU_functions) ---
    pund: Optional[pd.DataFrame] = None         # V, I_sw, P, branch
    dhm: Optional[pd.DataFrame] = None          # V, I, P, loop
    extracted_params: dict = field(default_factory=dict)
    instrument_params: dict = field(default_factory=dict)

    # --- labelling ---
    sample: Optional[str] = None
    _plot_string_override: str = ""

    @property
    def plot_string(self) -> str:
        """Legend label, built from the metadata unless explicitly overridden.

        The instrument's own ``SampleName`` is deliberately left out: it is a
        long free-text field ('SM04_vacann_40umetched_isolated_d1_100kwakeup')
        that swamps a legend. It stays in ``metadata['sample']``. Pass
        ``sample='SM04'`` to the importer for a short prefix instead.
        """
        if self._plot_string_override:
            return self._plot_string_override

        parts = []
        if self.sample:
            parts.append(str(self.sample))
        parts.append(f"Run{self.run_number}")
        # A few metadata fields that meaningfully distinguish runs
        for key, fmt in (('table_no', 'T{}'),
                         ('cycles', 'N={:.3g}'),
                         ('amplitude_V', '{:.3g} V'),
                         ('frequency_Hz', '{:.4g} Hz')):
            val = self.metadata.get(key)
            if val is not None:
                try:
                    parts.append(fmt.format(val))
                except (ValueError, TypeError):
                    parts.append(str(val))

        # Delay between the switching and non-switching pulse. It sets how much
        # polarisation relaxes back before the read pulse arrives, so two PUND
        # runs that differ only in delay are not comparable without it.
        delay = self.metadata.get('delay_s')
        if delay is not None:
            parts.append(f"td={_fmt_time(delay)}")

        return _clean_plot_string(', '.join(parts))

    @plot_string.setter
    def plot_string(self, value: str):
        self._plot_string_override = value if value else self._plot_string_override

    @property
    def thickness_cm(self) -> Optional[float]:
        """Film thickness in cm, for electric-field conversion."""
        if self.film_thickness_nm is None:
            return None
        return float(self.film_thickness_nm) * 1e-7

    def __repr__(self):
        return (f"PMUdata({self.meas_type}, run={self.run_number}, "
                f"instrument='{self.instrument}', n={len(self.raw_data)}, "
                f"'{self.plot_string}')")


@dataclass
class FatigueRun:
    """One fatigue experiment - an aixACCT 'Result Table'.

    A fatigue measurement interrupts the cycling at a series of cycle counts and
    records a full ferroelectric loop at each.  ``summary`` is the instrument's
    own table of extracted values against cycle number; ``sub_runs`` holds the
    loops themselves as ordinary :class:`PMUdata` objects, so every DHM/PUND
    plotting function works on them unmodified.
    """

    # --- identity ---
    file_name: str
    run_number: int
    folder_path: Path
    instrument: str
    loop_type: str                  # 'DHM' | 'PUND' - what was measured each step

    # --- as measured ---
    summary: pd.DataFrame           # cycles + the instrument's Pr/Vc/... values
    sub_runs: list = field(default_factory=list)     # list[PMUdata], one per cycle point
    metadata: dict = field(default_factory=dict)

    # --- device geometry ---
    electrode_area: Optional[float] = None
    film_thickness_nm: Optional[float] = None

    # --- transformed ---
    summary_extracted: Optional[pd.DataFrame] = None
    extracted_params: dict = field(default_factory=dict)

    # --- labelling ---
    sample: Optional[str] = None
    _plot_string_override: str = ""

    @property
    def plot_string(self) -> str:
        if self._plot_string_override:
            return self._plot_string_override
        parts = []
        if self.sample:
            parts.append(str(self.sample))
        parts.append(f"Run{self.run_number}")
        parts.append(self.loop_type)
        # The cycle count actually reached, not the one that was configured -
        # a run that was stopped early would otherwise be mislabelled.
        if self.summary is not None and 'cycles' in self.summary \
                and len(self.summary):
            parts.append(f"to {self.summary['cycles'].max():.3g} cycles")
        elif self.metadata.get('total_cycles') is not None:
            parts.append(f"{self.metadata['total_cycles']:.3g} cycles")
        return _clean_plot_string(', '.join(parts))

    @plot_string.setter
    def plot_string(self, value: str):
        self._plot_string_override = value if value else self._plot_string_override

    def __repr__(self):
        return (f"FatigueRun(run={self.run_number}, loop_type='{self.loop_type}', "
                f"points={len(self.summary)}, sub_runs={len(self.sub_runs)}, "
                f"'{self.plot_string}')")


# =============================================================================
#  Container base class
# =============================================================================

class PMUContainer(ABC):
    """Abstract base for all PMU import classes.

    Subclasses implement :meth:`_load_data`, which populates
    ``self.measurements``: a dict of run number -> :class:`PMUdata` (or
    :class:`FatigueRun`).  Run numbers are assigned sequentially from 0 as files
    are read; the originating file and table are always recorded in the run's
    ``metadata`` so the mapping stays traceable.
    """

    #: measurement type handled by the subclass, e.g. 'PUND'
    meas_type: str = ''

    def __init__(self,
                 root_folder: str,
                 folder_name: Optional[str] = None,
                 electrode_area: Optional[float] = None,
                 film_thickness_nm: Optional[float] = None,
                 sample: Optional[str] = None,
                 verbose: bool = True):
        self.root_folder: Path = Path(root_folder)
        self.folder_name: Optional[str] = folder_name
        self.folder_path: Path = (self.root_folder / folder_name
                                  if folder_name else self.root_folder)
        # Geometry given here always wins over anything found in the files
        self.electrode_area = electrode_area
        self.film_thickness_nm = film_thickness_nm
        self.sample = sample
        self.verbose = verbose

        self.measurements: dict = {}
        self._run_counter: int = 0

    # ---- container protocol -------------------------------------------------

    def __iter__(self):
        """Iterate directly over the measurement objects."""
        return iter(self.measurements.values())

    def __getitem__(self, run_number: int):
        """Retrieve a measurement by run number, with a helpful error."""
        if run_number not in self.measurements:
            raise KeyError(
                f"Run {run_number} not found. Available: "
                f"{sorted(self.measurements)}")
        return self.measurements[run_number]

    def __len__(self):
        return len(self.measurements)

    def __repr__(self):
        return (f"{type(self).__name__}(folder_path='{self.folder_path}', "
                f"runs={len(self.measurements)})")

    # ---- index-matched views (convenient in notebooks) ----------------------

    @property
    def run_nums(self) -> list:
        return list(self.measurements.keys())

    @property
    def runs(self) -> list:
        return list(self.measurements.values())

    @property
    def metadata(self) -> list:
        return [m.metadata for m in self]

    @property
    def raw_data(self) -> list:
        return [m.raw_data for m in self]

    @property
    def extracted_params(self) -> list:
        return [m.extracted_params for m in self]

    @property
    def plot_strings(self) -> list:
        return [m.plot_string for m in self]

    # ---- selection ----------------------------------------------------------

    def get(self, run_nums: Optional[list] = None) -> tuple:
        """Return ``(runs, plot_strings)`` for the given run numbers, or all."""
        if run_nums is None:
            selected = list(self.measurements.values())
        else:
            selected = [self.measurements[r] for r in run_nums
                        if r in self.measurements]
            missing = [r for r in run_nums if r not in self.measurements]
            if missing:
                print(f"Warning: runs {missing} not found in {self!r}")
        return selected, [m.plot_string for m in selected]

    def find(self, **criteria) -> list:
        """Select runs whose metadata matches every ``key=value`` given.

        Values may be a scalar or a container of accepted values, e.g.
        ``pund.find(sample='SM04_vacann_40umetch', table_no=[3, 4])``.
        """
        out = []
        for m in self:
            ok = True
            for key, want in criteria.items():
                have = m.metadata.get(key, getattr(m, key, None))
                if isinstance(want, (list, tuple, set)):
                    ok = ok and (have in want)
                else:
                    ok = ok and (have == want)
            if ok:
                out.append(m)
        return out

    def summary_table(self) -> pd.DataFrame:
        """One row per run: identity plus the headline metadata. For orientation."""
        rows = []
        for m in self:
            row = {'run': m.run_number,
                   'file': Path(m.file_name).name,
                   'label': m.plot_string}
            row.update({k: v for k, v in m.metadata.items()
                        if isinstance(v, (int, float, str))})
            rows.append(row)
        return pd.DataFrame(rows)

    # ---- shared helpers -----------------------------------------------------

    @abstractmethod
    def _load_data(self):
        """Populate ``self.measurements``. Implemented by each subclass."""

    def _transform_data(self):
        """Process every run into its loop and derived parameters.

        Called by each subclass's ``__init__`` straight after ``_load_data``, so
        that by the time an object exists its runs already carry ``.pund`` /
        ``.dhm`` and ``.extracted_params``. Plotting functions then only read
        stored data rather than computing it, and the same numbers back the
        plot, the summary table and the CSV export.

        The import from ``PMU_functions`` is deliberately deferred to call time:
        that module imports the dataclasses from here, so a module-level import
        would be circular.
        """
        return

    def _next_run(self) -> int:
        run = self._run_counter
        self._run_counter += 1
        return run

    def _resolve_geometry(self,
                          file_area: Optional[float],
                          file_thickness: Optional[float],
                          context: str = '') -> tuple:
        """Decide the geometry for one run.

        Constructor arguments win over anything found in the file, because the
        file value describes the pad the instrument was told about, which is not
        always the pad that was actually probed.  If neither is available the
        caller is warned once - polarisation scales inversely with area and the
        electric field inversely with thickness, so silently defaulting would
        produce plausible-looking but wrong numbers.
        """
        area = self.electrode_area if self.electrode_area is not None else file_area
        thick = (self.film_thickness_nm if self.film_thickness_nm is not None
                 else file_thickness)

        if area is None and not getattr(self, '_warned_area', False):
            print(f"Warning: no electrode area for {context or self.folder_path.name}. "
                  f"Polarisation cannot be computed - pass electrode_area (m^2) "
                  f"to {type(self).__name__}().")
            self._warned_area = True
        if thick is None and not getattr(self, '_warned_thickness', False):
            print(f"Note: no film thickness for {context or self.folder_path.name}; "
                  f"electric-field axes will be unavailable. Pass "
                  f"film_thickness_nm to {type(self).__name__}() to enable them.")
            self._warned_thickness = True

        return area, thick

    @staticmethod
    def _extract_value(filename: str, pattern: str, default=None):
        """Pull a numeric value out of a filename via a regex with one group."""
        match = re.search(pattern, filename)
        return float(match.group(1)) if match else default

    @staticmethod
    def _extract_int(text: str, pattern: str, default=None):
        """Pull an integer out of a string via a regex with one group."""
        match = re.search(pattern, text)
        return int(match.group(1)) if match else default

    @staticmethod
    def _extract_string(filename: str, patterns: tuple, default=None):
        """Return the first of ``patterns`` that appears in ``filename``."""
        low = filename.lower()
        return next((p for p in patterns if p.lower() in low), default)

    def _report(self):
        if not self.verbose:
            return
        if not self.measurements:
            print(f"{type(self).__name__}: no {self.meas_type} data found in "
                  f"{self.folder_path}")
            return
        print(f"{type(self).__name__}: loaded {len(self.measurements)} "
              f"{self.meas_type} run(s) from {self.folder_path}")


# =============================================================================
#  aixACCT TFAnalyzer .dat parsing
# =============================================================================
#
# All three aixACCT exports share one block grammar:
#
#     <block name>                 - a bare line, no colon
#     key: value                   - metadata, repeated
#     col1 <tab> col2 <tab> ...    - a column header (non-numeric fields)
#     1.0e+000 <tab> ...           - numeric rows, repeated
#
# The file opens with a summary section (one row per measurement), then a block
# carrying ``TfaModule:`` which marks the start of the per-measurement data, then
# one block per measurement.  Everything before the TfaModule block is summary;
# everything after is data.

_AIX_KIND_BY_HEADER = {
    'pulseresult':              'PM',    # PUND
    'dynamichysteresisresult':  'DHM',   # dynamic hysteresis
    'fatigue':                  'FM',    # fatigue
}


def _aix_sniff(path: Path) -> Optional[str]:
    """Return the aixACCT module for a .dat file: 'PM', 'DHM', 'FM' or None.

    Detection is by file content, not filename: the exports are not consistently
    named (``NewPund.dat`` holds the same thing as ``PundFile1.dat``).  The first
    line names the file type; ``TfaModule:`` in the header confirms it.
    """
    try:
        with open(path, encoding='latin-1') as fh:
            first = fh.readline().strip().lower()
            kind = _AIX_KIND_BY_HEADER.get(first)
            if kind is not None:
                return kind
            # Fall back to TfaModule if the banner is unfamiliar
            for _ in range(2000):
                line = fh.readline()
                if not line:
                    break
                if line.startswith('TfaModule:'):
                    return line.split(':', 1)[1].strip().upper()
    except OSError as exc:
        print(f"Warning: could not read {path.name}: {exc}")
    return None


def _is_numeric_row(fields: list) -> bool:
    """True if every non-empty field parses as a float."""
    seen = False
    for f in fields:
        f = f.strip()
        if not f:
            continue
        try:
            float(f)
        except ValueError:
            return False
        seen = True
    return seen


def _aix_read_blocks(path: Path) -> tuple:
    """Parse an aixACCT .dat file into its blocks.

    Returns ``(file_header, summary_blocks, data_blocks)`` where each block is a
    dict with keys ``name``, ``meta``, ``columns`` and ``data`` (an ndarray, or
    None when the block is metadata only).  ``file_header`` is the metadata of
    the block containing ``TfaModule``.
    """
    blocks = []
    cur = None

    def _flush():
        nonlocal cur
        if cur is not None:
            cur['data'] = (np.asarray(cur.pop('rows'), dtype=float)
                           if cur['rows'] else None)
            blocks.append(cur)
        cur = None

    def _new(name):
        return {'name': name, 'meta': {}, 'columns': None, 'rows': []}

    with open(path, encoding='latin-1') as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            if not line.strip():
                continue

            if '\t' in line:
                fields = line.rstrip('\t').split('\t')
                if cur is None:
                    cur = _new('')
                if _is_numeric_row(fields):
                    try:
                        cur['rows'].append([float(f) if f.strip() else np.nan
                                            for f in fields])
                    except ValueError:
                        pass          # malformed row - skip rather than abort
                else:
                    cur['columns'] = [f.strip() for f in fields]
                continue

            if ':' in line:
                key, _, val = line.partition(':')
                if cur is None:
                    cur = _new('')
                cur['meta'][key.strip()] = val.strip()
                continue

            # A bare line starts a new block
            _flush()
            cur = _new(line.strip())

    _flush()

    # Split at the block that carries TfaModule
    file_header, split_at = {}, len(blocks)
    for i, b in enumerate(blocks):
        if 'TfaModule' in b['meta']:
            file_header, split_at = b['meta'], i
            break

    summary = [b for b in blocks[:split_at] if b['data'] is not None]
    data = [b for b in blocks[split_at + 1:]]
    return file_header, summary, data


_AIX_NUMERIC_KEYS = (
    'Area [mm2]', 'Thickness [nm]', 'Pund Frequency [Hz]', 'Pund Amplitude [V]',
    'Hysteresis Frequency [Hz]', 'Hysteresis Amplitude [V]',
    'Fatigue Frequency [Hz]', 'Fatigue Amplitude [V]', 'Fatigue Offset [V]',
    'Total Cycles', 'Pulse Points', 'Number of pulses',
    'Write Pulse Time [s]', 'Write Pulse Rise Time [s]',
    'Write Pulse Amplitude [V]', 'Write Pulse Delay [s]', 'Read Pulse Delay [s]',
    'Pr+ [uC/cm2]', 'Pr- [uC/cm2]', 'Prrel+ [uC/cm2]', 'Prrel- [uC/cm2]',
    'Psw [uC/cm2]', 'Pnsw [uC/cm2]', 'dPsw [uC/cm2]', 'Px [uC/cm2]',
    'Pvmax+ [uC/cm2]', 'Pvmax- [uC/cm2]', 'Vc+ [V]', 'Vc- [V]', 'VcShift [V]',
    'Vmax+ [V]', 'Vmax- [V]', 'Ipk+ [A]', 'Ipk- [A]', 'Wloss [uJ/cm2]',
    'Rav [Ohm]', 'Cls [F]', 'Epsls [1]', 'Measurement Status',
)

#: instrument-reported scalars worth keeping, mapped to short names
_AIX_PARAM_MAP = {
    'Pr+ [uC/cm2]': 'Pr+', 'Pr- [uC/cm2]': 'Pr-',
    'Prrel+ [uC/cm2]': 'Prrel+', 'Prrel- [uC/cm2]': 'Prrel-',
    'Psw [uC/cm2]': 'Psw', 'Pnsw [uC/cm2]': 'Pnsw', 'dPsw [uC/cm2]': 'dPsw',
    'Pvmax+ [uC/cm2]': 'Pvmax+', 'Pvmax- [uC/cm2]': 'Pvmax-',
    'Vc+ [V]': 'Vc+', 'Vc- [V]': 'Vc-', 'VcShift [V]': 'VcShift',
    'Vmax+ [V]': 'Vmax+', 'Vmax- [V]': 'Vmax-',
    'Ipk+ [A]': 'Ipk+', 'Ipk- [A]': 'Ipk-',
    'Wloss [uJ/cm2]': 'Wloss', 'Rav [Ohm]': 'Rav',
    'Cls [F]': 'Cls', 'Epsls [1]': 'Epsls',
}


def _aix_meta(raw_meta: dict) -> dict:
    """Normalise an aixACCT metadata block: numbers as floats, tidy key names."""
    out = {}
    for key, val in raw_meta.items():
        if key in _AIX_NUMERIC_KEYS:
            try:
                out[key] = float(val.split()[0]) if val else np.nan
            except (ValueError, IndexError):
                out[key] = val
        else:
            out[key] = val

    # Short, stable aliases used by the plotting layer and plot_string
    if 'SampleName' in out:
        out['sample'] = out['SampleName']
    if isinstance(out.get('Area [mm2]'), float):
        out['area_m2'] = out['Area [mm2]'] * 1e-6        # mm^2 -> m^2
    if isinstance(out.get('Thickness [nm]'), float):
        out['thickness_nm'] = out['Thickness [nm]']
    for src, dst in (('Pund Amplitude [V]', 'amplitude_V'),
                     ('Hysteresis Amplitude [V]', 'amplitude_V'),
                     ('Pund Frequency [Hz]', 'frequency_Hz'),
                     ('Hysteresis Frequency [Hz]', 'frequency_Hz'),
                     ('Read Pulse Delay [s]', 'delay_s'),
                     ('Total Cycles', 'cycles')):
        if isinstance(out.get(src), float):
            out.setdefault(dst, out[src])
    return out


def _aix_instrument_params(meta: dict) -> dict:
    """Pull the instrument's own extracted scalars out of a metadata block."""
    return {short: meta[long] for long, short in _AIX_PARAM_MAP.items()
            if isinstance(meta.get(long), float)}


def _aix_pund_frame(columns: list, data: np.ndarray) -> pd.DataFrame:
    """Reshape an aixACCT PUND table into a single time-ordered trace.

    The export lays the pulses out side by side: ``Time V I P`` repeated once per
    pulse.  They are stacked back into one trace, with a ``pulse`` column (1-based)
    recording which pulse each sample came from, so the P/U/N/D split is exact
    and does not depend on thresholding the voltage.
    """
    n_groups = len(columns) // 4
    frames = []
    for g in range(n_groups):
        block = data[:, 4 * g:4 * g + 4]
        frames.append(pd.DataFrame({
            'Time': block[:, 0],
            'V': block[:, 1],
            'I': block[:, 2],
            'P_instr': block[:, 3],
            'pulse': g + 1,
        }))
    df = pd.concat(frames, ignore_index=True)
    return df.dropna(subset=['Time', 'V', 'I']).reset_index(drop=True)


def _aix_dhm_frame(columns: list, data: np.ndarray) -> pd.DataFrame:
    """Reshape an aixACCT DHM / fatigue loop table.

    Columns are ``Time, V+, V-, I1, P1, I2, P2, I3, P3``.  ``V-`` is the negated
    drive (aixACCT records both polarities); loop 1 is the main hysteresis loop
    and loops 2-3 are the relaxed-remanence loops.
    """
    names = ['Time', 'V_pos', 'V_neg']
    n_loops = (data.shape[1] - 3) // 2
    for k in range(1, n_loops + 1):
        names += [f'I{k}', f'P{k}']
    names = names[:data.shape[1]]
    df = pd.DataFrame(data[:, :len(names)], columns=names)
    return df.dropna(subset=['Time', 'V_pos']).reset_index(drop=True)


# =============================================================================
#  aixACCT import classes
# =============================================================================

class _AixACCTBase(PMUContainer):
    """Shared machinery for the three aixACCT importers."""

    #: aixACCT module code this class consumes
    module: str = ''

    def _dat_files(self) -> list:
        """.dat files in the folder whose sniffed module matches this class."""
        if not self.folder_path.is_dir():
            print(f"Error: {self.folder_path} is not a directory.")
            return []
        files = sorted(f for f in self.folder_path.iterdir()
                       if f.suffix.lower() == '.dat' and not f.name.startswith('~'))
        matched = [f for f in files if _aix_sniff(f) == self.module]
        if not matched and files:
            found = {f.name: _aix_sniff(f) for f in files}
            print(f"Note: no {self.module} files in {self.folder_path}. "
                  f"Found: {found}")
        return matched


class AixACCT_PUND(_AixACCTBase):
    """aixACCT TFAnalyzer PUND export (``TfaModule: PM``).

    Each ``Table N`` in the file becomes one run.  The pulse train is five pulses,
    ``P U N D P`` (polarity + + - - +): the fifth pulse is a second positive
    switching pulse, and because it follows D the film is definitively
    down-poled beforehand.  The first P has no such guarantee - its prior state
    depends on measurement history - so ``p_pulse='last'`` is the default.
    """

    meas_type = 'PUND'
    module = 'PM'

    def __init__(self, root_folder, folder_name=None, electrode_area=None,
                 film_thickness_nm=None, sample=None, p_pulse: str = 'last',
                 centre_branches: bool = True, med_filt: Optional[int] = None,
                 sigma_filt: Optional[int] = None, transform: bool = True,
                 verbose: bool = True):
        super().__init__(root_folder, folder_name, electrode_area,
                         film_thickness_nm, sample, verbose)
        if p_pulse not in ('last', 'first', 'both'):
            raise ValueError("p_pulse must be 'last', 'first' or 'both'")
        self.p_pulse = p_pulse
        self.centre_branches = centre_branches
        self.med_filt = med_filt
        self.sigma_filt = sigma_filt
        self._load_data()
        if transform:
            self._transform_data()
        self._report()

    def _transform_data(self):
        from PMU_functions import extract_pund
        for meas in self:
            extract_pund(meas, p_pulse=self.p_pulse,
                         centre_branches=self.centre_branches,
                         med_filt=self.med_filt, sigma_filt=self.sigma_filt,
                         verbose=self.verbose)

    def _load_data(self):
        for path in self._dat_files():
            _, _, blocks = _aix_read_blocks(path)
            for block in blocks:
                if block['data'] is None or block['columns'] is None:
                    continue
                if not block['name'].lower().startswith('table'):
                    continue

                meta = _aix_meta(block['meta'])
                meta['source_file'] = path.name
                meta['table_no'] = self._extract_int(block['name'], r'(\d+)')
                meta['p_pulse'] = self.p_pulse

                df = _aix_pund_frame(block['columns'], block['data'])
                if df.empty:
                    continue

                area, thick = self._resolve_geometry(
                    meta.get('area_m2'), meta.get('thickness_nm'),
                    context=f"{path.name} {block['name']}")

                run = self._next_run()
                self.measurements[run] = PMUdata(
                    file_name=str(path), run_number=run,
                    folder_path=self.folder_path,
                    instrument='aixacct', meas_type='PUND',
                    raw_data=df, metadata=meta,
                    electrode_area=area, film_thickness_nm=thick,
                    instrument_params=_aix_instrument_params(meta),
                    sample=self.sample)


class AixACCT_DHM(_AixACCTBase):
    """aixACCT dynamic hysteresis export (``TfaModule: DHM``).

    Each ``Table N`` holds three measured loops.  Loop 1 is the main hysteresis
    loop - its first sample equals the reported ``Pr-`` exactly - while loops 2
    and 3 are the relaxed-remanence loops (``Prrel-`` and ``Prrel+``).
    """

    meas_type = 'DHM'
    module = 'DHM'

    def __init__(self, root_folder, folder_name=None, electrode_area=None,
                 film_thickness_nm=None, sample=None, loop: int = 1,
                 transform: bool = True, verbose: bool = True):
        super().__init__(root_folder, folder_name, electrode_area,
                         film_thickness_nm, sample, verbose)
        self.loop = loop
        self._load_data()
        if transform:
            self._transform_data()
        self._report()

    def _transform_data(self):
        from PMU_functions import extract_dhm
        for meas in self:
            extract_dhm(meas, loop=self.loop, verbose=self.verbose)

    def _load_data(self):
        for path in self._dat_files():
            _, _, blocks = _aix_read_blocks(path)
            for block in blocks:
                if block['data'] is None or block['columns'] is None:
                    continue
                if not block['name'].lower().startswith('table'):
                    continue

                meta = _aix_meta(block['meta'])
                meta['source_file'] = path.name
                meta['table_no'] = self._extract_int(block['name'], r'(\d+)')

                df = _aix_dhm_frame(block['columns'], block['data'])
                if df.empty:
                    continue

                area, thick = self._resolve_geometry(
                    meta.get('area_m2'), meta.get('thickness_nm'),
                    context=f"{path.name} {block['name']}")

                run = self._next_run()
                self.measurements[run] = PMUdata(
                    file_name=str(path), run_number=run,
                    folder_path=self.folder_path,
                    instrument='aixacct', meas_type='DHM',
                    raw_data=df, metadata=meta,
                    electrode_area=area, film_thickness_nm=thick,
                    instrument_params=_aix_instrument_params(meta),
                    sample=self.sample)


class AixACCT_Fatigue(_AixACCTBase):
    """aixACCT fatigue export (``TfaModule: FM``).

    Each ``Result Table`` becomes one :class:`FatigueRun`.  A result table holds
    a summary of the instrument's extracted values against cycle number, plus one
    ``Data Table [R,k]`` per cycle point containing the full loop measured there.
    Those loops are stored as ordinary :class:`PMUdata` objects in ``sub_runs``,
    so the DHM and PUND plotting functions work on them directly.

    The loop type is read from the summary column prefix (``1-DHM Pr+ [uC/cm2]``
    -> DHM), since a fatigue run may interrogate the film with either DHM or PUND.
    """

    meas_type = 'Fatigue'
    module = 'FM'

    def __init__(self, root_folder, folder_name=None, electrode_area=None,
                 film_thickness_nm=None, sample=None, Pr_mode: str = 'raw',
                 loop: int = 1, transform: bool = True, verbose: bool = True):
        super().__init__(root_folder, folder_name, electrode_area,
                         film_thickness_nm, sample, verbose)
        self.Pr_mode = Pr_mode
        self.loop = loop
        self._load_data()
        if transform:
            self._transform_data()
        self._report()

    def _transform_data(self):
        """Build each run's Pr-vs-cycles table and process every stored loop.

        The sub-runs are transformed too, so the individual loops behind a
        fatigue curve are ready to plot without a second pass - that is what
        makes ``plot_fatigue(..., sub_plots=True)`` and ad-hoc inspection of a
        single cycle point cheap.
        """
        from PMU_functions import extract_fatigue, extract_dhm, extract_pund
        for run in self:
            for sub in run.sub_runs:
                if run.loop_type == 'PUND':
                    extract_pund(sub, verbose=False)
                else:
                    extract_dhm(sub, loop=self.loop, verbose=False)
            extract_fatigue(run, Pr_mode=self.Pr_mode, verbose=self.verbose)

    @staticmethod
    def _parse_summary(columns: list, data: np.ndarray) -> tuple:
        """Build the cycles DataFrame and work out what was measured each step.

        Column names look like ``1-DHM Pr+ [uC/cm2]``.  The prefix gives the
        measurement index and type; it is stripped so the frame carries plain
        names (``Pr+``, ``Vc-``, ...).
        """
        loop_type, names = None, []
        for col in columns:
            match = re.match(r'^(\d+)-(\w+)\s+(.*)$', col)
            if match:
                loop_type = loop_type or match.group(2).upper()
                name = match.group(3)
            else:
                name = col
            names.append(re.sub(r'\s*\[.*?\]\s*$', '', name).strip())

        n = min(len(names), data.shape[1])
        df = pd.DataFrame(data[:, :n], columns=names[:n])
        if 'Cycles' in df.columns:
            df = df.rename(columns={'Cycles': 'cycles'})
        return df, (loop_type or 'DHM')

    def _load_data(self):
        for path in self._dat_files():
            _, _, blocks = _aix_read_blocks(path)

            current = None      # the FatigueRun being assembled
            for block in blocks:
                name = block['name']

                # --- a new fatigue experiment ---
                if name.lower().startswith('result table'):
                    meta = _aix_meta(block['meta'])
                    meta['source_file'] = path.name
                    meta['result_table'] = self._extract_int(name, r'(\d+)')
                    if isinstance(meta.get('Total Cycles'), float):
                        meta['total_cycles'] = meta['Total Cycles']

                    if block['data'] is None or block['columns'] is None:
                        print(f"Warning: '{name}' in {path.name} has no summary "
                              f"table - skipping.")
                        current = None
                        continue

                    summary, loop_type = self._parse_summary(
                        block['columns'], block['data'])
                    area, thick = self._resolve_geometry(
                        meta.get('area_m2'), meta.get('thickness_nm'),
                        context=f"{path.name} {name}")

                    run = self._next_run()
                    current = FatigueRun(
                        file_name=str(path), run_number=run,
                        folder_path=self.folder_path,
                        instrument='aixacct', loop_type=loop_type,
                        summary=summary, metadata=meta,
                        electrode_area=area, film_thickness_nm=thick,
                        sample=self.sample)
                    self.measurements[run] = current
                    continue

                # --- extra per-step parameters, folded into the run metadata ---
                if name.lower().startswith('data measurement parameters'):
                    if current is not None:
                        current.metadata['measurement_parameters'] = dict(block['meta'])
                    continue

                # --- one loop, measured at a given cycle count ---
                if name.lower().startswith('data table'):
                    if current is None or block['data'] is None:
                        continue
                    meta = _aix_meta(block['meta'])
                    meta['source_file'] = path.name
                    meta['data_table'] = name
                    meta['result_table'] = current.metadata.get('result_table')
                    meta['point_index'] = len(current.sub_runs) + 1
                    if isinstance(meta.get('Total Cycles'), float):
                        meta['cycles'] = meta['Total Cycles']

                    if current.loop_type == 'PUND':
                        df = _aix_pund_frame(block['columns'], block['data'])
                    else:
                        df = _aix_dhm_frame(block['columns'], block['data'])
                    if df.empty:
                        continue

                    sub = PMUdata(
                        file_name=str(path),
                        run_number=meta['point_index'],
                        folder_path=self.folder_path,
                        instrument='aixacct', meas_type=current.loop_type,
                        raw_data=df, metadata=meta,
                        electrode_area=current.electrode_area,
                        film_thickness_nm=current.film_thickness_nm,
                        instrument_params=_aix_instrument_params(meta),
                        sample=current.sample)
                    current.sub_runs.append(sub)

    def _report(self):
        if not self.verbose:
            return
        if not self.measurements:
            print(f"AixACCT_Fatigue: no fatigue data found in {self.folder_path}")
            return
        print(f"AixACCT_Fatigue: loaded {len(self.measurements)} fatigue run(s) "
              f"from {self.folder_path}")
        for run in self:
            cycles = run.summary['cycles'] if 'cycles' in run.summary else None
            span = (f"{cycles.min():.3g} to {cycles.max():.3g} cycles"
                    if cycles is not None and len(cycles) else 'unknown range')
            print(f"  run {run.run_number}: {run.loop_type}, "
                  f"{len(run.summary)} points ({span}), "
                  f"{len(run.sub_runs)} loops")


# =============================================================================
#  Keithley 4200A import classes
# =============================================================================

class _KeithleyBase(PMUContainer):
    """Shared machinery for the Keithley 4200A Clarius .xlsx exports."""

    #: substrings that identify a file as belonging to this class
    filename_keys: tuple = ()

    def _xlsx_files(self) -> list:
        if not self.folder_path.is_dir():
            print(f"Error: {self.folder_path} is not a directory.")
            return []
        return sorted(f for f in self.folder_path.iterdir()
                      if f.suffix.lower() == '.xlsx'
                      and not f.name.startswith('~$')
                      and any(k in f.name.lower() for k in self.filename_keys))

    @staticmethod
    def _run_from_name(name: str) -> Optional[int]:
        match = re.search(r'Run(\d+)', name)
        return int(match.group(1)) if match else None


class Keithley_PUND(_KeithleyBase):
    """Keithley 4200A-SCS Clarius PUND export (the ``nvm`` pundTest module).

    Files are ``.xlsx`` with 'pund' in the name; data sheets are those carrying a
    ``V`` column.  The waveform is four trapezoidal pulses (P U N D) with no
    preset, so the state of the film before the P pulse depends on measurement
    history - see the note in ``PMU_functions.extract_pund``.

    Clarius records no device geometry, so ``electrode_area`` (m^2) must be
    supplied; ``film_thickness_nm`` is optional and enables the field axes.
    """

    meas_type = 'PUND'
    filename_keys = ('pund',)

    def __init__(self, root_folder, folder_name=None, electrode_area=None,
                 film_thickness_nm=None, sample=None,
                 centre_branches: bool = True, med_filt: Optional[int] = None,
                 sigma_filt: Optional[int] = None, transform: bool = True,
                 verbose: bool = True):
        super().__init__(root_folder, folder_name, electrode_area,
                         film_thickness_nm, sample, verbose)
        self.centre_branches = centre_branches
        self.med_filt = med_filt
        self.sigma_filt = sigma_filt
        self._load_data()
        if transform:
            self._transform_data()
        self._report()

    def _transform_data(self):
        from PMU_functions import extract_pund
        for meas in self:
            extract_pund(meas, centre_branches=self.centre_branches,
                         med_filt=self.med_filt, sigma_filt=self.sigma_filt,
                         verbose=self.verbose)

    def _load_data(self):
        for path in self._xlsx_files():
            try:
                xls = pd.ExcelFile(path)
            except Exception as exc:
                print(f"Error opening {path.name}: {exc}")
                continue

            # Haoran-format files also match 'pund' - they belong to the other class
            if 'PUND_Diff' in xls.sheet_names:
                if self.verbose:
                    print(f"  skipping {path.name} (Haoran format - use "
                          f"Keithley_Haoran_PUND)")
                continue

            file_run = self._run_from_name(path.name)
            for sheet in xls.sheet_names:
                try:
                    peek = pd.read_excel(xls, sheet_name=sheet, header=0, nrows=1)
                except Exception:
                    continue
                if 'V' not in peek.columns:
                    continue

                df_raw = pd.read_excel(xls, sheet_name=sheet, header=0)
                cols = {}
                for src, dst in (('t', 'Time'), ('V', 'V'), ('I', 'I')):
                    if src in df_raw.columns:
                        cols[dst] = pd.to_numeric(df_raw[src], errors='coerce')
                if not {'Time', 'V', 'I'} <= set(cols):
                    print(f"Warning: {path.name}[{sheet}] lacks t/V/I - skipping.")
                    continue

                df = pd.DataFrame(cols).dropna().reset_index(drop=True)
                if df.empty:
                    continue

                # Clarius writes its own scalar results in row 0 of P/Psw/Qsw,
                # with NaN below - keep them as scalars, not columns.
                instr = {}
                for src, dst in (('P', 'P'), ('Psw', 'Psw'), ('Qsw', 'Qsw')):
                    if src in df_raw.columns:
                        vals = pd.to_numeric(df_raw[src], errors='coerce').dropna()
                        if len(vals):
                            instr[dst] = float(vals.iloc[0])

                meta = {'source_file': path.name, 'sheet': sheet}
                area, thick = self._resolve_geometry(None, None,
                                                     context=path.name)

                run = self._next_run()
                self.measurements[run] = PMUdata(
                    file_name=str(path), run_number=run,
                    folder_path=self.folder_path,
                    instrument='keithley', meas_type='PUND',
                    raw_data=df, metadata=meta,
                    electrode_area=area, film_thickness_nm=thick,
                    instrument_params=instr,
                    sample=self.sample)
                # Keep the Clarius run number visible in the label
                if file_run is not None:
                    meta['clarius_run'] = file_run
                    self.measurements[run].plot_string = (
                        f"{self.sample + ' ' if self.sample else ''}Run{file_run}")


class Keithley_Haoran_PUND(_KeithleyBase):
    """Keithley PUND exported by the Haoran acquisition script.

    These workbooks already contain the processed PUND: the ``PUND_Diff`` sheet
    holds the difference current, charge, polarisation, branch and segment for
    both current ranges, so nothing is recomputed here.  ``Parameters`` carries
    the waveform settings and the electrode area (``area_cm2``).

    Two current ranges are recorded.  ``channel='I2'`` (the finer range, the one
    used for the script's own ``loop_i2diff`` plot) is the default;
    ``channel='I1'`` selects the coarse range.
    """

    meas_type = 'PUND'
    filename_keys = ('pund',)

    def __init__(self, root_folder, folder_name=None, electrode_area=None,
                 film_thickness_nm=None, sample=None, channel: str = 'I2',
                 transform: bool = True, verbose: bool = True):
        super().__init__(root_folder, folder_name, electrode_area,
                         film_thickness_nm, sample, verbose)
        if channel not in ('I1', 'I2'):
            raise ValueError("channel must be 'I1' or 'I2'")
        self.channel = channel
        self._load_data()
        if transform:
            self._transform_data()
        self._report()

    def _transform_data(self):
        # The loop itself is already in the file; this only fills in Pr, Vc,
        # Ec and imprint from it.
        from PMU_functions import extract_pund
        for meas in self:
            extract_pund(meas, verbose=self.verbose)

    def _load_data(self):
        for path in self._xlsx_files():
            try:
                xls = pd.ExcelFile(path)
            except Exception as exc:
                print(f"Error opening {path.name}: {exc}")
                continue
            if 'PUND_Diff' not in xls.sheet_names:
                continue

            meta = {'source_file': path.name}

            # --- waveform settings and geometry ---
            file_area = None
            if 'Parameters' in xls.sheet_names:
                par = pd.read_excel(xls, sheet_name='Parameters', header=0)
                if {'name', 'value'} <= set(par.columns):
                    for _, r in par.iterrows():
                        try:
                            meta[str(r['name'])] = float(r['value'])
                        except (TypeError, ValueError):
                            meta[str(r['name'])] = r['value']
                if isinstance(meta.get('area_cm2'), float):
                    file_area = meta['area_cm2'] * 1e-4       # cm^2 -> m^2
                if isinstance(meta.get('Vp'), float):
                    meta['amplitude_V'] = meta['Vp']
                if isinstance(meta.get('delay_time'), float):
                    meta['delay_s'] = meta['delay_time']

            # --- raw trace ---
            if 'Total' in xls.sheet_names:
                tot = pd.read_excel(xls, sheet_name='Total', header=0)
                raw = pd.DataFrame({
                    'Time': pd.to_numeric(tot.get('Time'), errors='coerce'),
                    'V': pd.to_numeric(tot.get('Voltage'), errors='coerce'),
                    'I1': pd.to_numeric(tot.get('CurrentI1'), errors='coerce'),
                    'I2': pd.to_numeric(tot.get('CurrentI2'), errors='coerce'),
                }).dropna().reset_index(drop=True)
                raw['I'] = raw[self.channel]
            else:
                raw = pd.DataFrame(columns=['Time', 'V', 'I'])

            # --- already-processed PUND ---
            diff = pd.read_excel(xls, sheet_name='PUND_Diff', header=0)
            i_col = 'DiffCurrent' if self.channel == 'I2' else 'DiffCurrentI1'
            p_col = 'Polarization' if self.channel == 'I2' else 'PolarizationI1'
            q_col = 'Charge' if self.channel == 'I2' else 'ChargeI1'
            missing = [c for c in (i_col, p_col, q_col) if c not in diff.columns]
            if missing:
                print(f"Warning: {path.name} PUND_Diff lacks {missing} - skipping.")
                continue

            pund = pd.DataFrame({
                'Time': pd.to_numeric(diff['Time'], errors='coerce'),
                'V': pd.to_numeric(diff['Voltage'], errors='coerce'),
                'I_sw': pd.to_numeric(diff[i_col], errors='coerce'),
                'Q': pd.to_numeric(diff[q_col], errors='coerce'),
                'P': pd.to_numeric(diff[p_col], errors='coerce'),
                'branch': diff.get('Segment', pd.Series(dtype=object)),
                'direction': diff.get('Branch', pd.Series(dtype=object)),
            }).dropna(subset=['V', 'P']).reset_index(drop=True)
            meta['channel'] = self.channel
            # The acquisition script already did the P-U / N-D subtraction and
            # the integration; extract_pund must not redo it from the raw trace.
            meta['pund_precomputed'] = True

            area, thick = self._resolve_geometry(file_area, None,
                                                 context=path.name)
            run = self._next_run()
            self.measurements[run] = PMUdata(
                file_name=str(path), run_number=run,
                folder_path=self.folder_path,
                instrument='keithley_haoran', meas_type='PUND',
                raw_data=raw, metadata=meta, pund=pund,
                electrode_area=area, film_thickness_nm=thick,
                sample=self.sample)


class Keithley_Custom(_KeithleyBase):
    """Keithley 4200A custom-waveform export (pulsed FeFET measurement).

    Channel 2 drives the gate with a PUND-like pulse train while channel 1
    monitors the drain, so this is a 3-terminal measurement and belongs in
    ``PMU_Analysis_3T.ipynb``.
    """

    meas_type = 'CUSTOM'
    filename_keys = ('customwaveform', 'custom_waveform')

    def __init__(self, root_folder, folder_name=None, electrode_area=None,
                 film_thickness_nm=None, sample=None, verbose: bool = True):
        super().__init__(root_folder, folder_name, electrode_area,
                         film_thickness_nm, sample, verbose)
        self._load_data()
        self._report()

    def _load_data(self):
        for path in self._xlsx_files():
            try:
                xls = pd.ExcelFile(path)
            except Exception as exc:
                print(f"Error opening {path.name}: {exc}")
                continue

            file_run = self._run_from_name(path.name)
            for sheet in xls.sheet_names:
                try:
                    peek = pd.read_excel(xls, sheet_name=sheet, header=0, nrows=1)
                except Exception:
                    continue
                if 'VMeasCh1' not in peek.columns:
                    continue

                df_raw = pd.read_excel(xls, sheet_name=sheet, header=0)
                cols = {}
                for src, dst in (('TimeOutput', 'Time'),
                                 ('VMeasCh1', 'V_Ch1'), ('IMeasCh1', 'I_Ch1'),
                                 ('VMeasCh2', 'V_Ch2'), ('IMeasCh2', 'I_Ch2')):
                    if src in df_raw.columns:
                        cols[dst] = pd.to_numeric(df_raw[src], errors='coerce')

                df = pd.DataFrame(cols).dropna(subset=['V_Ch1']).reset_index(drop=True)
                if df.empty:
                    continue

                meta = {'source_file': path.name, 'sheet': sheet}
                if file_run is not None:
                    meta['clarius_run'] = file_run

                run = self._next_run()
                meas = PMUdata(
                    file_name=str(path), run_number=run,
                    folder_path=self.folder_path,
                    instrument='keithley', meas_type='CUSTOM',
                    raw_data=df, metadata=meta,
                    sample=self.sample)
                if file_run is not None:
                    meas.plot_string = (
                        f"{self.sample + ' ' if self.sample else ''}Run{file_run}")
                self.measurements[run] = meas


class Keithley_DHM(_KeithleyBase):
    """Placeholder for Keithley 4200A DHM (dynamic hysteresis) exports.

    Not implemented: the only DHM workbooks available at the time of writing are
    in a folder marked 'wrongly exported', so there is no known-good file to
    write the reader against.  Once a good export exists, mirror
    :class:`Keithley_PUND` - read the data sheet, map its columns onto
    ``Time``/``V``/``I``/``P`` and store ``meas_type='DHM'``; the DHM transform
    and plotting functions will then work unchanged.
    """

    meas_type = 'DHM'
    filename_keys = ('dhm',)

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "Keithley_DHM is not implemented yet - no known-good export exists. "
            "Use AixACCT_DHM, or supply a valid Clarius DHM workbook.")

    def _load_data(self):
        raise NotImplementedError


# =============================================================================
#  Combining objects
# =============================================================================

def merge_PMU(*objects, sample: Optional[str] = None):
    """Combine several PMU containers of the same measurement type into one.

    Run numbers are reassigned sequentially in the order the objects are given,
    so the result can be indexed and plotted like any single import.
    """
    objects = [o for o in objects if o is not None]
    if not objects:
        raise ValueError("merge_PMU needs at least one object")

    types = {o.meas_type for o in objects}
    if len(types) > 1:
        raise ValueError(f"cannot merge different measurement types: {types}")

    merged = copy.copy(objects[0])
    merged.measurements = {}
    merged._run_counter = 0
    merged.sample = sample or objects[0].sample

    for obj in objects:
        for meas in obj:
            run = merged._next_run()
            new = copy.copy(meas)
            new.run_number = run
            if not new._plot_string_override:
                new._plot_string_override = meas.plot_string
            merged.measurements[run] = new

    print(f"merge_PMU: combined {len(objects)} objects into "
          f"{len(merged.measurements)} runs")
    return merged
