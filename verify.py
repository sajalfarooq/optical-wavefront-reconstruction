
"""Verify the retained release without training or modifying any files."""
from pathlib import Path
import hashlib
import json
ROOT=Path(__file__).resolve().parent

def verify():
    manifest=json.loads((ROOT/"published/manifest.json").read_text())
    for name,expected in manifest["files"].items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT): raise ValueError("Invalid manifest path")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError("Retained artifact changed: "+name)
    import sys
    sys.path.insert(0,str(ROOT/"src"))
    import numpy as np
    from optics_portfolio.imaging import aligned_phase_rmse
    with np.load(ROOT/"results/study_data.npz") as a:
        value,_=aligned_phase_rmse(a["estimate"],a["truth"],a["mask"])
    expected=json.loads((ROOT/"results/metrics.json").read_text())["noise_and_plane_study"][1]["phase_rmse_rad"]
    np.testing.assert_allclose(value,expected,atol=1e-12)
    return len(manifest["files"])

if __name__=="__main__": print("Verified",verify(),"retained artifact files; no retraining")
