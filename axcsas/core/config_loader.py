"""
AXCSAS Configuration Loader
===========================

Centralized configuration loading from YAML files.
集中式 YAML 配置載入。
"""

from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

from axcsas.core.constants import (
    CU_KA1,
    SCHERRER_K,
    MIN_RELIABLE_SIZE,
    MAX_RELIABLE_SIZE,
)


def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    從 YAML 檔案載入配置。
    
    Falls back to defaults from constants module if no config file.
    若無配置檔案則使用 constants 模組的預設值。
    
    Args:
        config_path: Path to config file. 配置檔案路徑。
        
    Returns:
        Configuration dictionary. 配置字典。
    """
    default_config = {
        "physical_constants": {
            "wavelength": CU_KA1,
            "scherrer_k": {
                "spherical": SCHERRER_K.spherical,
                "cubic": SCHERRER_K.cubic,
                "default": SCHERRER_K.default,
            },
        },
        "validation": {
            "size_limits": {
                "min_reliable": MIN_RELIABLE_SIZE,
                "max_reliable": MAX_RELIABLE_SIZE,
            },
        },
    }
    
    if config_path is None:
        config_path = Path("config.yaml")
    
    if YAML_AVAILABLE and config_path.exists():
        with open(config_path, encoding="utf-8") as f:
            user_config = yaml.safe_load(f)
        if user_config:
            return _deep_merge(default_config, user_config)
    
    return default_config


def _deep_merge(base: dict, override: dict) -> dict:
    """
    Deep merge two dictionaries.
    深度合併兩個字典。
    """
    result = base.copy()
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result
