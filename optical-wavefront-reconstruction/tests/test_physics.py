import numpy as np
import pytest
def relative_error(a,b): return float(np.linalg.norm(a-b)/np.linalg.norm(b))
from optics_portfolio.imaging import fresnel,retrieve,aligned_phase_rmse

def test_propagation_round_trip_and_plane_wave():
    rng = np.random.default_rng(7)
    field = rng.normal(size=(32, 32))+1j*rng.normal(size=(32, 32))
    forward = fresnel(field, 8e-6, 532e-9, 0.03)
    np.testing.assert_allclose(np.sum(abs(field)**2), np.sum(abs(forward)**2), rtol=1e-14)
    np.testing.assert_allclose(fresnel(forward, 8e-6, 532e-9, -0.03), field, atol=3e-15)
    np.testing.assert_allclose(fresnel(np.ones((32,32)), 8e-6, 532e-9, 0.03), 1)

def test_phase_retrieval_reduces_residual():
    x = (np.arange(32)-16)/10
    xx, yy = np.meshgrid(x, x)
    amp = np.exp(-xx*xx-yy*yy)
    truth = amp*np.exp(0.3j*(xx*xx-yy*yy))
    zs = [0.02, 0.04, 0.06]
    intensities = np.array([abs(fresnel(truth, 8e-6, 532e-9, z))**2 for z in zs])
    estimate, losses = retrieve(amp, intensities, zs, 8e-6, 532e-9, iterations=150)
    rmse, _ = aligned_phase_rmse(estimate, truth, amp > 0.2)
    assert losses[-1] < losses[0]*0.05
    assert rmse < 0.08

def test_fresnel_infinite_gaussian_beam():
    n=128; dx=8e-6; wavelength=532e-9; z=.01; waist=80e-6
    x=(np.arange(n)-n/2)*dx
    xx,yy=np.meshgrid(x,x)
    initial=np.exp(-(xx*xx+yy*yy)/waist**2)
    q=1+1j*wavelength*z/(np.pi*waist**2)
    exact=np.exp(-(xx*xx+yy*yy)/(waist**2*q))/q
    assert relative_error(fresnel(initial,dx,wavelength,z),exact)<1e-10

def test_phase_error_removes_global_piston():
    truth=np.ones((16,16),complex)*np.exp(.4j)
    value,_=aligned_phase_rmse(truth*np.exp(1.3j),truth,np.ones((16,16),bool))
    assert value<1e-14

