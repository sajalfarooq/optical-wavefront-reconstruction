"""Descriptive numerical summaries."""
import numpy as np

def aggregate(values):
    values = np.asarray(values, float)
    return {"n": int(len(values)), "mean": float(values.mean()),
            "std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            "minimum": float(values.min()), "maximum": float(values.max())}
