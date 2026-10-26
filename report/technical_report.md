# Multi-plane wavefront reconstruction technical report

## Research question and findings

This study asks how observation-plane diversity, photon noise, propagation calibration and sampling affect phase retrieval when object amplitude is known. A small hard-aperture baseline is retained for teaching, while the extended experiment uses a larger field of view and a smooth pupil taper after a sampling audit exposed substantial window dependence.

Three-plane retrieval improves mean phase error at higher photon counts in the retained extended experiment. At lower counts it performs worse than one plane under the tested averaging estimator. Measurement diversity provides information, but an estimator's noise handling, constraints and convergence can limit how effectively that information is used. The study preserves both results rather than assuming a monotonic benefit from more planes.

## Propagation model and sampling

The scalar paraxial Fresnel operator is

$$P_z[u]=\mathcal F^{-1}\{\mathcal F[u]\exp[-i\pi\lambda z(f_x^2+f_y^2)]\}.$$

A constant carrier phase is omitted. The operator is discrete, periodic and unitary. Negative distance is its adjoint/inverse. Analytical infinite-domain Gaussian-beam propagation checks the transfer-function signs and length scales independently of a forward/backward round trip.

The original 96×96 hard aperture produced about 10.4% intensity discrepancy when the computational window was doubled. This exposed optimistic matched-grid results. The extended case uses 192×192 samples at 8 micrometre pitch, wavelength 532 nm, and a Gaussian pupil with a raised-cosine taper over the outer annulus. Observation distances are 25,55,85 mm. A fourth plane at 70 mm is never used in reconstruction and evaluates predicted intensity.

Increasing the field of view and smoothing the pupil changes the physical numerical experiment; the improved audit is not presented as a correction that leaves the original object unchanged. The original sampling audit is retained separately when available. The new pupil still imposes a known-amplitude prior.

## Inverse method

Each iteration propagates the object estimate to every measured plane, replaces propagated amplitude by the measured square-root intensity, keeps phase, and back-propagates. Back-propagated fields are averaged and projected onto the known object amplitude. The standard budget is 250 iterations. This educational alternating-projection method is nonconvex and does not certify a unique or globally optimal solution.

Photon counts are Poisson with expected peak levels 10,000 and 1,000 per image's shared scaling. Three seeds each change the noise draw and initial phase. Errors therefore reflect both sources jointly. Phase is undefined in dark regions; evaluation uses object amplitude above 0.2 and removes global piston. The RMS is computed on wrapped differences, not on an unwrapped optical-path map.

## Retained repeated-run results

| Peak expected photoelectrons | Planes | Mean phase RMS rad | Sample std rad | Mean held-out intensity relative error |
|---|---:|---:|---:|---:|
| 10000 | 1 | 0.010750 | 0.000534 | 0.005503 |
| 10000 | 3 | 0.009209 | 0.000110 | 0.004599 |
| 1000 | 1 | 0.031125 | 0.001649 | 0.016453 |
| 1000 | 3 | 0.039190 | 0.000745 | 0.021185 |

The enlarged-window intensity discrepancy is **0.1240%**. Halving pixel pitch at a shared field of view gives **0.0081%** discrepancy at common sample coordinates. These are sensitivity measurements, not rigorous continuum-error bounds.


## Model mismatch and interpretation

The wavelength sensitivity experiment reconstructs noiseless intensities with the assumed wavelength shifted by ±1%. A low residual under the wrong model can coexist with a biased wavefront, so matched-model reconstruction error is not sufficient for experimental accuracy claims. The held-out plane tests predictive consistency but does not remove all model ambiguities.

The near-machine-precision result in the original noiseless baseline uses the same discrete forward model for generation and reconstruction and exact object amplitude. It is an algorithmic consistency case, not a physical precision estimate. Independent calibration, finite detector pixels, read noise and amplitude uncertainty would all matter in the laboratory.

## Reproduction and limitations

Run `python reproduce.py --profile standard`. Three notebooks cover Fresnel sampling, phase retrieval, and noise/model-error interpretation. Saved metrics include individual seeds, summaries, the held-out plane, wavelength mismatch and sampling comparisons. Noisy reconstruction and robustness figures accompany the ideal baseline figures.

The implementation excludes high-NA vector diffraction, partial coherence, pixel integration, dark current, detector saturation, and unknown object amplitude. Multi-start selection is not performed, and a joint change of noise and initialization is not a factorial uncertainty decomposition. More observations need a noise-aware estimator to guarantee no performance loss; this simple projection scheme provides no such guarantee.

## References

1. J. W. Goodman, *Introduction to Fourier Optics*, 4th ed., W. H. Freeman / Macmillan Learning (2017). [Publisher](https://www.macmillanlearning.com/ed/uk/product/Introduction-to-Fourier-Optics-4th-edition/p/1319119166). Background for scalar diffraction, Fresnel propagation and Fourier-domain optical modelling.
2. J. R. Fienup, “Phase retrieval algorithms: a comparison,” *Applied Optics* **21**(15), 2758–2769 (1982). [DOI](https://doi.org/10.1364/AO.21.002758). Iterative phase-retrieval constraints, error reduction and algorithmic limitations.
3. R. W. Gerchberg and W. O. Saxton, “A practical algorithm for the determination of phase from image and diffraction plane pictures,” *Optik* **35**, 237–246 (1972). Foundational alternating amplitude/phase projection method. The multi-plane averaging estimator used here is an educational variant rather than a reproduction of the original experiment.
