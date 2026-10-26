"""Paraxial Fresnel propagation and multi-plane alternating projections.

Lengths in metres. FFT boundary is periodic; negative z is the numerical
adjoint. This is coherent scalar optics, not vector or high-NA imaging.
"""
import numpy as np


def fresnel(field, dx, wavelength, z):
    field = np.asarray(field, complex)
    if field.ndim != 2 or min(field.shape,default=0)<2 or dx <= 0 or wavelength <= 0:
        raise ValueError("A 2D field and positive sampling/wavelength are required")
    if not all(np.all(np.isfinite(v)) for v in (field,dx,wavelength,z)):
        raise ValueError("Finite propagation inputs are required")
    fy = np.fft.fftfreq(field.shape[0], dx)
    fx = np.fft.fftfreq(field.shape[1], dx)
    h = np.exp(-1j * np.pi * wavelength * z * (fy[:, None]**2 + fx[None, :]**2))
    return np.fft.ifft2(np.fft.fft2(field) * h)


def retrieve(amplitude, intensities, distances, dx, wavelength, iterations=250, seed=7):
    amplitude = np.asarray(amplitude, float)
    intensities = np.asarray(intensities, float)
    distances = np.asarray(distances,float)
    if amplitude.ndim != 2 or not all(np.all(np.isfinite(v)) for v in (amplitude,intensities,distances)):
        raise ValueError("Finite 2D amplitude, intensities and distances are required")
    if np.sum(amplitude**2)==0 or np.any(np.sum(intensities,axis=(-2,-1))==0):
        raise ValueError("Nonzero illumination and measurements are required")
    if intensities.shape != (len(distances), *amplitude.shape) or np.any(intensities < 0):
        raise ValueError("One nonnegative intensity image per propagation distance required")
    if np.any(amplitude < 0) or iterations < 1 or not len(distances):
        raise ValueError("Nonnegative known amplitude, distances and iterations required")
    rng = np.random.default_rng(seed)
    field = amplitude * np.exp(1j * rng.normal(0, 0.15, amplitude.shape))
    losses = []
    measured = np.sqrt(intensities)
    for _ in range(iterations):
        estimates, loss = [], 0.0
        for amp, z in zip(measured, distances):
            out = fresnel(field, dx, wavelength, z)
            loss += np.linalg.norm(np.abs(out) - amp)**2 / max(np.linalg.norm(amp)**2, 1e-30)
            constrained = amp * np.exp(1j * np.angle(out))
            estimates.append(fresnel(constrained, dx, wavelength, -z))
        average = np.mean(estimates, axis=0)
        field = amplitude * np.exp(1j * np.angle(average))
        losses.append(float(loss / len(distances)))
    return field, np.array(losses)


def aligned_phase_rmse(estimate, truth, mask):
    if np.shape(estimate)!=np.shape(truth) or np.shape(mask)!=np.shape(truth) or not np.any(mask):
        raise ValueError("Matching fields and a nonempty mask are required")
    offset = np.angle(np.sum(estimate[mask] * np.conj(truth[mask])))
    difference = np.angle(estimate * np.exp(-1j * offset) * np.conj(truth))
    return float(np.sqrt(np.mean(difference[mask]**2))), difference
