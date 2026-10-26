"""Repeated reconstructions, held-out planes and model mismatch."""
from pathlib import Path
import numpy as np
from optics_portfolio.imaging import fresnel,retrieve,aligned_phase_rmse
from optics_portfolio.profiles import profile
from optics_portfolio.validation import aggregate
from optics_portfolio.output import result_dir,save_json,plt


def object_field(n=192,dx=8e-6,apodized=True):
    x=(np.arange(n)-n/2)*dx
    xx,yy=np.meshgrid(x,x); radius=.25e-3
    r2=(xx*xx+yy*yy)/radius**2
    # Raised-cosine apodization removes the abrupt pupil-amplitude jump.
    taper=np.where(r2<.5,1,np.where(r2<1,.5*(1+np.cos(np.pi*(r2-.5)/.5)),0))
    amp=np.exp(-1.2*r2)*(taper if apodized else r2<1)
    phase=1.2*r2+.7*(xx*xx-yy*yy)/radius**2+.5*xx/radius*r2
    return amp,amp*np.exp(1j*phase)


def main(name="standard"):
    cfg=profile(name)
    out=result_dir(Path(__file__).parent); out.mkdir(parents=True,exist_ok=True)
    dx,lam=8e-6,532e-9
    amp,truth=object_field(); mask=amp>.2
    zs=[.025,.055,.085]; held_z=.07
    clean=np.array([abs(fresnel(truth,dx,lam,z))**2 for z in zs])
    held_intensity=abs(fresnel(truth,dx,lam,held_z))**2
    rows=[]; shown=None
    for noise_seed in cfg["seeds"]:
        rng=np.random.default_rng(noise_seed)
        for photons in [10000,1000]:
            observed=rng.poisson(clean/clean.max()*photons)*clean.max()/photons
            for planes in [1,3]:
                estimate,loss=retrieve(amp,observed[:planes],zs[:planes],dx,lam,
                                      iterations=cfg["phase_iterations"],seed=noise_seed)
                phase_error,difference=aligned_phase_rmse(estimate,truth,mask)
                predicted=abs(fresnel(estimate,dx,lam,held_z))**2
                rows.append({"seed":noise_seed,"planes":planes,"peak_photons":photons,
                    "phase_rmse_rad":phase_error,
                    "held_out_intensity_relative_error":float(np.linalg.norm(predicted-held_intensity)/np.linalg.norm(held_intensity)),
                    "last_pre_update_residual":float(loss[-1])})
                if photons==1000 and planes==3 and noise_seed==cfg["seeds"][0]:
                    shown=(estimate,difference)
    summaries=[]
    for photons in [10000,1000]:
        for planes in [1,3]:
            selected=[r for r in rows if r["peak_photons"]==photons and r["planes"]==planes]
            summaries.append({"peak_photons":photons,"planes":planes,
                "phase":aggregate([r["phase_rmse_rad"] for r in selected]),
                "held_out_intensity":aggregate([r["held_out_intensity_relative_error"] for r in selected])})
    mismatch=[]
    for wavelength_scale in [.99,1.,1.01]:
        estimate,_=retrieve(amp,clean,zs,dx,lam*wavelength_scale,iterations=cfg["phase_iterations"],seed=7)
        err,_=aligned_phase_rmse(estimate,truth,mask)
        predicted=abs(fresnel(estimate,dx,lam*wavelength_scale,held_z))**2
        mismatch.append({"wavelength_scale":wavelength_scale,"phase_rmse_rad":err,
                         "held_out_intensity_error":float(np.linalg.norm(predicted-held_intensity)/np.linalg.norm(held_intensity))})
    # Same physical sampling and object, enlarged zero-padded computational domain.
    padding=truth.shape[0]//2
    padded=np.pad(truth,((padding,padding),(padding,padding)))
    pad_plane=abs(fresnel(padded,dx,lam,.085))[padding:-padding,padding:-padding]**2
    discrepancy=float(np.linalg.norm(pad_plane-clean[-1])/np.linalg.norm(pad_plane))
    # Same field of view at half pitch: compare amplitudes at common sample coordinates.
    _,fine_truth=object_field(384,dx/2)
    fine_plane=abs(fresnel(fine_truth,dx/2,lam,.085))[::2,::2]**2
    resolution_error=float(np.linalg.norm(fine_plane-clean[-1])/np.linalg.norm(fine_plane))
    fig,axes=plt.subplots(1,3,figsize=(13,3.9),layout="constrained")
    estimated,difference=shown
    piston=np.angle(np.sum(estimated[mask]*np.conj(truth[mask])))
    for ax,val,title in zip(axes,[np.angle(truth),np.angle(estimated*np.exp(-1j*piston)),difference],
                            ["True phase","Three-plane, 1,000 photons","Piston-aligned phase error"]):
        im=ax.imshow(np.ma.array(val,mask=~mask),origin="lower",cmap="RdBu_r")
        ax.set(xlabel="Pixel x",ylabel="Pixel y",title=title)
        fig.colorbar(im,ax=ax,label="Radians")
    fig.savefig(out/"noisy_reconstruction.png"); plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout="constrained")
    for planes in [1,3]:
        sel=[s for s in summaries if s["planes"]==planes]
        axes[0].errorbar([s["peak_photons"] for s in sel],[s["phase"]["mean"] for s in sel],
                         yerr=[s["phase"]["std"] for s in sel],fmt="o-",label=f"{planes} plane(s)")
    axes[0].set(xscale="log",xlabel="Peak expected photoelectrons",ylabel="Phase RMS (rad)",title="Across noise/initialization seeds")
    axes[0].legend()
    axes[1].plot([(r["wavelength_scale"]-1)*100 for r in mismatch],[r["phase_rmse_rad"] for r in mismatch],"o-")
    axes[1].set(xlabel="Assumed wavelength error (%)",ylabel="Phase RMS (rad)",title="Noiseless model mismatch")
    fig.savefig(out/"robustness.png"); plt.close(fig)
    save_json(out/"advanced_metrics.json",{"profile":name,"config":cfg,"runs":rows,"summary":summaries,
        "wavelength_mismatch":mismatch,"held_out_distance_m":held_z,
        "grid":[192,192],"pixel_pitch_m":dx,"pupil":"Gaussian with raised-cosine taper over 0.5<r_squared/radius_squared<1",
        "padding_intensity_discrepancy":discrepancy,"half_pitch_intensity_discrepancy":resolution_error,
        "sampling_note":"Apodized pupil and larger field of view reduce aliasing; discrepancies remain sensitivity measurements, not certified error bounds.",
        "noise_note":"Each seed changes both Poisson data and starting phase; error bars mix both sources, not independent estimates of each."})
    np.savez_compressed(out/"advanced_data.npz",truth=truth,noisy_estimate=estimated,mask=mask,
                        held_out_intensity=held_intensity,padded_plane=pad_plane,fine_plane=fine_plane)
    print("Extended wavefront study complete",flush=True)


if __name__=="__main__":
    import os
    main(os.environ.get("OPTICS_PROFILE","standard"))
