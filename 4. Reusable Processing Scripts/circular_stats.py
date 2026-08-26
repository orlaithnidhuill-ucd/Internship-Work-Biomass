"""
circular_stats.py
Circular (directional) statistics for wrapped phase data. Linear mean/std
are WRONG for phase — e.g. mean of [-179, 179] degrees should be 180, not
0 — so every averaging operation on phi/CPD values must go through the
complex-exponential (phasor) route implemented here.

Consolidated from repeated `circ`/`circ_deg`/`circ_smooth` definitions
across the El Nino case-study notebook.

Usage:
    from circular_stats import circular_stats, circ_smooth

    stats = circular_stats(phi_array)   # dict: n, mean_deg, sigma_deg, coh
    smoothed = circ_smooth(phi_2d, sigma=6.0)
"""
import numpy as np
from scipy.ndimage import gaussian_filter


def circular_stats(vals):
    """
    Circular mean, angular spread, and mean resultant length (coherence)
    of a wrapped-phase array, in radians.

    Returns dict: n, mean_rad, mean_deg, sigma_deg, coh
      - mean_rad / mean_deg : circular mean
      - sigma_deg           : angular standard deviation, derived from the
                               mean resultant length R via sqrt(-2 ln R)
      - coh                 : R, the mean resultant length (0-1). Low R
                               means the phase values are scattered/noisy,
                               not that the "coherence" in a SAR-processing
                               sense is necessarily low — treat this as a
                               dispersion statistic on the angle samples.
    NaNs are dropped before computing. Returns all-NaN dict for empty input.
    """
    v = np.asarray(vals)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return dict(n=0, mean_rad=np.nan, mean_deg=np.nan,
                    sigma_deg=np.nan, coh=np.nan)
    z = np.mean(np.exp(1j * v))
    R = float(np.abs(z))
    mean_rad = float(np.angle(z))
    sigma_rad = float(np.sqrt(max(0.0, -2.0 * np.log(max(R, 1e-12)))))
    return dict(
        n=int(v.size),
        mean_rad=round(mean_rad, 6),
        mean_deg=round(np.degrees(mean_rad), 4),
        sigma_deg=round(np.degrees(sigma_rad), 3),
        coh=round(R, 4),
    )


def circ_smooth(phi, sigma, min_weight=0.25):
    """
    NaN-aware Gaussian smoothing of a 2D wrapped-phase array, done
    correctly by smoothing the real/imaginary parts of the unit phasor
    separately (never smooth angles directly — that produces wraparound
    artefacts at the +-180 degree boundary).

    Parameters
    ----------
    phi : 2D array of phase in radians, NaN where invalid
    sigma : Gaussian smoothing radius, in pixels
    min_weight : pixels where the smoothed valid-data weight falls below
                 this fraction are set back to NaN, since they're mostly
                 extrapolated from neighbours rather than real data

    Returns smoothed phase array, same shape, NaN where under-supported.
    """
    ok = np.isfinite(phi)
    zr = gaussian_filter(np.where(ok, np.cos(phi), 0.0), sigma, mode='nearest')
    zi = gaussian_filter(np.where(ok, np.sin(phi), 0.0), sigma, mode='nearest')
    wt = gaussian_filter(ok.astype(float), sigma, mode='nearest')
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.arctan2(zi / wt, zr / wt)
    return np.where(wt > min_weight, out, np.nan)


def circular_mean_deg(vals):
    """Convenience: circular mean only, in degrees. NaN if input is empty."""
    return circular_stats(vals)['mean_deg']
