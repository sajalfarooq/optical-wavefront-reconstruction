from pathlib import Path
import numpy as np
from optics_portfolio.imaging import fresnel, retrieve, aligned_phase_rmse
from optics_portfolio.output import setup, write_metrics, plt


def main():
    out = setup(Path(__file__).parent)
    n, dx, wavelength = 96, 8e-6, 532e-9
    x = (np.arange(n)-n/2)*dx
    xx, yy = np.meshgrid(x, x)
    radius = 0.25e-3
    r2 = (xx*xx+yy*yy)/radius**2
    amplitude = np.exp(-1.2*r2)*(r2 < 1)
    phase = 1.2*r2 + 0.7*(xx*xx-yy*yy)/radius**2 + 0.5*xx/radius*r2
    truth = amplitude*np.exp(1j*phase)
    distances = [0.025, 0.055, 0.085]
    clean = np.array([abs(fresnel(truth, dx, wavelength, z))**2 for z in distances])
    mask = amplitude > 0.2
    rng = np.random.default_rng(7)
    rows = []
    for photons in [None, 10000, 1000]:
        # Photon count is expected photoelectrons at the brightest sampled pixel.
        observed = clean if photons is None else rng.poisson(clean/clean.max()*photons)*clean.max()/photons
        for planes in [1, 3]:
            estimate, losses = retrieve(amplitude, observed[:planes], distances[:planes], dx, wavelength)
            rmse, difference = aligned_phase_rmse(estimate, truth, mask)
            rows.append({"peak_expected_photons": photons, "planes": planes,
                         "phase_rmse_rad": rmse, "initial_amplitude_loss": float(losses[0]),
                         "final_amplitude_loss": float(losses[-1])})
            if photons is None and planes == 3:
                best, best_difference, best_losses = estimate, difference, losses
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), layout="constrained")
    offset = np.angle(np.sum(best[mask]*np.conj(truth[mask])))
    phase_estimate = np.angle(best*np.exp(-1j*offset))
    for ax, value, title in zip(axes, [np.angle(truth), phase_estimate, best_difference],
                               ["True phase", "Reconstructed phase", "Wrapped phase error"]):
        im = ax.imshow(np.ma.array(value, mask=~mask), extent=[x[0]*1e3, x[-1]*1e3, x[0]*1e3, x[-1]*1e3],
                       origin="lower", cmap="RdBu_r")
        ax.set(xlabel="x (mm)", ylabel="y (mm)", title=title)
        fig.colorbar(im, ax=ax, label="Radians")
    fig.savefig(out / "phase_reconstruction.png")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
    ax.semilogy(best_losses)
    ax.set(xlabel="Iteration", ylabel="Relative amplitude residual", title="Three-plane noiseless convergence")
    fig.savefig(out / "convergence.png")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for planes in [1, 3]:
        selected = [row for row in rows if row["planes"] == planes]
        ax.semilogy([0, 1, 2], [row["phase_rmse_rad"] for row in selected], "o-", label=f"{planes} plane(s)")
    ax.set(xticks=[0, 1, 2], xticklabels=["Noiseless", "10,000", "1,000"],
           xlabel="Peak expected photoelectrons", ylabel="Phase RMSE (rad)", title="Noise and measurement diversity")
    ax.legend()
    fig.savefig(out / "noise_comparison.png")
    plt.close(fig)
    np.savez_compressed(out / "study_data.npz", truth=truth, estimate=best, intensity=clean,
                        distances_m=distances, dx_m=dx, wavelength_m=wavelength, mask=mask)
    write_metrics(out, {"seed": 7, "grid": [n, n], "pixel_pitch_m": dx,
        "wavelength_m": wavelength, "distances_m": distances, "iterations": 250,
        "noise_and_plane_study": rows, "assumptions": "known object amplitude, scalar paraxial periodic propagation",
        "phase_metric": "wrapped RMS after removing global piston on amplitude > 0.2",
        "uncertainty": "one noise realization and initialization per setting; no global uniqueness claim"})
    print("Wavefront study saved:", out)


if __name__ == "__main__":
    main()
