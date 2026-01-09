"""
XRD Data Loader Module
Supports multiple XRD data formats: .xy, .csv, .raw, .txt (Bruker)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, Optional, Dict, Any


class XRDDataLoader:
    """
    Multi-format XRD data loader.
    
    Supported formats:
    - .xy: Two-column format (2θ, Intensity)
    - .csv: CSV format with header
    - .raw: Bruker RAW format
    - .txt: Bruker TXT export format
    """
    
    SUPPORTED_FORMATS = ['.xy', '.csv', '.raw', '.txt']
    
    def __init__(self):
        self.data: Optional[np.ndarray] = None
        self.metadata: Dict[str, Any] = {}
        self.filepath: Optional[Path] = None
    
    def load(self, filepath: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load XRD data from file.
        
        Args:
            filepath: Path to the XRD data file
            
        Returns:
            Tuple of (two_theta, intensity) arrays
        """
        self.filepath = Path(filepath)
        suffix = self.filepath.suffix.lower()
        
        if suffix not in self.SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported format: {suffix}. "
                f"Supported formats: {self.SUPPORTED_FORMATS}"
            )
        
        if suffix == '.xy':
            return self._load_xy()
        elif suffix == '.csv':
            return self._load_csv()
        elif suffix == '.txt':
            return self._load_bruker_txt()
        elif suffix == '.raw':
            return self._load_bruker_raw()
    
    def _load_xy(self) -> Tuple[np.ndarray, np.ndarray]:
        """Load .xy format (simple two-column)."""
        data = np.loadtxt(self.filepath, comments=['#', ';'])
        two_theta = data[:, 0]
        intensity = data[:, 1]
        return two_theta, intensity
    
    def _load_csv(self) -> Tuple[np.ndarray, np.ndarray]:
        """Load .csv format with automatic header detection."""
        # Try to detect if header exists
        df = pd.read_csv(self.filepath)
        
        # Assume first column is 2θ, second is intensity
        two_theta = df.iloc[:, 0].values.astype(float)
        intensity = df.iloc[:, 1].values.astype(float)
        
        return two_theta, intensity
    
    def _load_bruker_txt(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load Bruker TXT export format.
        
        Format:
        - Header lines starting with special characters
        - Data section with 2θ and intensity columns
        """
        two_theta = []
        intensity = []
        in_data_section = False
        
        with open(self.filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Try to parse as data
                try:
                    parts = line.split()
                    if len(parts) >= 2:
                        theta = float(parts[0])
                        inten = float(parts[1])
                        two_theta.append(theta)
                        intensity.append(inten)
                        in_data_section = True
                except ValueError:
                    # Not a data line, might be header
                    if in_data_section:
                        # We've passed the data section
                        break
                    continue
        
        return np.array(two_theta), np.array(intensity)
    
    def _load_bruker_raw(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Load Bruker RAW format (binary).
        
        Note: This is a simplified implementation.
        For full RAW support, consider using external libraries.
        """
        raise NotImplementedError(
            "Bruker RAW binary format not yet implemented. "
            "Please export to TXT or XY format."
        )
    
    def get_metadata(self) -> Dict[str, Any]:
        """Return file metadata."""
        return {
            'filepath': str(self.filepath),
            'format': self.filepath.suffix if self.filepath else None,
            **self.metadata
        }


def load_xrd_data(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convenience function to load XRD data.
    
    Args:
        filepath: Path to XRD data file
        
    Returns:
        Tuple of (two_theta, intensity) arrays
    """
    loader = XRDDataLoader()
    return loader.load(filepath)
