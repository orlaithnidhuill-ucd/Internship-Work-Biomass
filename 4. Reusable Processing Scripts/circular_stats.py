# circular stats for wrapped phase; normal mean/std don't work here (mean of -179 and 179 should be 180, not 0), so everything goes through the complex exponential instead
import numpy as np
from scipy.ndimage import gaussian_filter

def circular_stats(vals):
    # circular mean, spread, and coherence of a wrapped phase array (radians)
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
    # gaussian smoothing on the unit phasor (real/imag parts separately),
    # not on the angle itself -- smoothing angles directly breaks at +-180
    ok = np.isfinite(phi)
    zr = gaussian_filter(np.where(ok, np.cos(phi), 0.0), sigma, mode='nearest')
    zi = gaussian_filter(np.where(ok, np.sin(phi), 0.0), sigma, mode='nearest')
    wt = gaussian_filter(ok.astype(float), sigma, mode='nearest')
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.arctan2(zi / wt, zr / wt)
    return np.where(wt > min_weight, out, np.nan)
    
def circular_mean_deg(vals):
    return circular_stats(vals)['mean_deg']
