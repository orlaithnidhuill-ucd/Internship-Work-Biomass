import re
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from pathlib import Path
from scipy.ndimage import uniform_filter, binary_dilation
from scipy.interpolate import griddata
from lxml import etree

try:
    import rasterio
    BACKEND = 'rasterio'
except Exception:
    import tifffile
    BACKEND = 'tifffile'

BASE     = Path(r'C:\Users\Orlaith.Doyle\STEP_BACK')
ARCHIVE  = BASE / 'PLUME_ARCHIVE' / 'amazon'
SMOS_DIR = ARCHIVE / 'smos'
WC_CACHE = BASE / 'worldcover_cache'
LABEL    = 'Amazon Delta'
SCENE_DATE = '2026-07-26'
MULTILOOK = (14, 1)
WC_WATER_CLASSES = {0, 80}
BAND_ORDER = ('VV', 'VH', 'HV', 'HH')

def is_complete(g):
    return (g.is_dir() and next(g.glob('*i_abs.tiff'), None)
            and next(g.glob('*i_phase.tiff'), None)
            and next(g.glob('*.kml'), None))

granules = sorted(g for g in ARCHIVE.iterdir() if is_complete(g))
print(f'{len(granules)} granule(s) on disk for {LABEL}')
for g in granules:
    print(' ', g.name)
if not granules:
    raise RuntimeError('nothing on disk -- put the granule folders in ARCHIVE first')

def read_bands(path):
    if BACKEND == 'rasterio':
        with rasterio.open(path) as src:
            return {BAND_ORDER[b]: src.read(b + 1) for b in range(4)}
    arr = tifffile.imread(path)
    return {BAND_ORDER[b]: arr[b] for b in range(4)}

def corner_grid(granule_dir, n_lines, n_samples):
    kml = next(granule_dir.glob('*.kml'))
    m = re.search(r'<coordinates>(.*?)</coordinates>', kml.read_text(), re.S)
    pts = np.array([tuple(map(float, p.split(',')))
                    for p in m.group(1).split()])[:4, :2]
    edge_len = lambda a, b: np.hypot(*(pts[b] - pts[a]))
    if (edge_len(0,1)+edge_len(2,3)) < (edge_len(1,2)+edge_len(3,0)):
        row_groups, col_groups = ({0,1},{2,3}), ({1,2},{0,3})
    else:
        row_groups, col_groups = ({1,2},{3,0}), ({0,1},{2,3})
    row_lats = [(g, np.mean([pts[i][1] for i in g])) for g in row_groups]
    south = min(row_lats, key=lambda x: x[1])[0]
    north = max(row_lats, key=lambda x: x[1])[0]
    ann = etree.parse(str(next(granule_dir.glob('*annot.xml'))))
    orbit_pass = next((e.text.strip() for e in ann.iter()
                       if etree.QName(e).localname == 'orbitPass' and e.text), 'Ascending')
    row0, rowM = (north, south) if orbit_pass == 'Descending' else (south, north)
    col_lons = [(g, np.mean([pts[i][0] for i in g])) for g in col_groups]
    col0 = max(col_lons, key=lambda x: x[1])[0]
    colM = min(col_lons, key=lambda x: x[1])[0]
    A, B = pts[list(row0 & col0)[0]], pts[list(row0 & colM)[0]]
    C, D = pts[list(rowM & col0)[0]], pts[list(rowM & colM)[0]]
    fr = np.linspace(0, 1, n_lines)[:, None, None]
    fc = np.linspace(0, 1, n_samples)[None, :, None]
    g = (A*(1-fc)+B*fc)*(1-fr) + (C*(1-fc)+D*fc)*fr
    return g[..., 0], g[..., 1]

def ocean_mask_for(lon, lat, shape):
    wc_file = next(WC_CACHE.glob(f'*{LABEL.split()[0]}*.tif'), None)
    if wc_file is None:
        print('  no WorldCover tile cached, skipping land mask')
        return np.ones(shape, dtype=bool)
    with rasterio.open(wc_file) as src:
        wc_arr = src.read(1)
    mask = np.isin(wc_arr, list(WC_WATER_CLASSES))
    return binary_dilation(mask, iterations=2)

def pauli_components(hh, hv, vv):
    return hh + vv, hh - vv, 2 * hv

def multilook_amplitude(slc, looksa, looksr):
    return np.sqrt(uniform_filter(np.abs(slc) ** 2, [looksa, looksr]))

def pauli_rgb(k1, k2, k3, looksa, looksr, clip_sigma=2.5):
    amp1 = multilook_amplitude(k1, looksa, looksr)
    amp2 = multilook_amplitude(k2, looksa, looksr)
    amp3 = multilook_amplitude(k3, looksa, looksr)
    rgb = np.stack([amp2, amp3, amp1], axis=-1).astype('float32')
    for b in range(3):
        rgb[:, :, b] = np.clip(rgb[:, :, b], 0, clip_sigma * np.mean(rgb[:, :, b]))
        rgb[:, :, b] /= np.max(rgb[:, :, b])
    return rgb

def alpha_angle(k1, k2, k3):
    plen = np.sqrt(np.abs(k1)**2 + np.abs(k2)**2 + np.abs(k3)**2)
    return np.arccos(np.abs(k1) / plen)

frames = []
for nom in granules:
    amp = read_bands(next(nom.glob('*i_abs.tiff')))
    pha = read_bands(next(nom.glob('*i_phase.tiff')))
    SLC = {p: amp[p] * np.exp(1j * pha[p]) for p in BAND_ORDER}
    del amp, pha

    n_lines, n_samples = SLC['HH'].shape
    lon, lat = corner_grid(nom, n_lines, n_samples)
    ocean = ocean_mask_for(lon, lat, SLC['HH'].shape)

    # Pauli, alpha, 4 pols; one figure per granule
    k1, k2, k3 = pauli_components(SLC['HH'], SLC['HV'], SLC['VV'])
    rgb = pauli_rgb(k1, k2, k3, *MULTILOOK)
    alpha = alpha_angle(k1, k2, k3)

    fig, ax = plt.subplots(1, 6, figsize=(24, 4))
    ax[0].imshow(rgb, aspect='auto'); ax[0].set_title('Pauli RGB')
    ax[1].imshow(np.degrees(alpha), cmap='viridis', vmin=0, vmax=90, aspect='auto')
    ax[1].set_title('alpha (deg)')
    for a, pol in zip(ax[2:], BAND_ORDER):
        a.imshow(20*np.log10(np.abs(SLC[pol]) + 1e-6), cmap='gray', aspect='auto')
        a.set_title(pol)
    fig.suptitle(nom.name)
    plt.tight_layout()
    plt.show()

    # HH-VV phase, ocean only, with histogram
    phi = np.angle(SLC['HH'] * np.conj(SLC['VV']))
    phi_ocean = np.where(ocean, phi, np.nan)
    v = phi_ocean[np.isfinite(phi_ocean)]
    z = np.mean(np.exp(1j*v)); R = float(np.abs(z))
    cmean = float(np.angle(z))
    sigma = float(np.sqrt(max(0, -2*np.log(max(R, 1e-12)))))

    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].imshow(phi_ocean, cmap='gist_rainbow', vmin=cmean-sigma, vmax=cmean+sigma, aspect='auto')
    ax[0].set_title('HH-VV phase, ocean only')
    ax[1].hist(np.degrees(v), bins=90, color='steelblue')
    ax[1].axvspan(np.degrees(cmean-sigma), np.degrees(cmean+sigma), color='gold', alpha=0.3)
    ax[1].axvline(np.degrees(cmean), color='k', lw=1)
    ax[1].set_title(f'phase histogram, mean {np.degrees(cmean):+.2f} deg')
    ax[1].set_xlabel('phase (deg)')
    plt.tight_layout()
    plt.show()

    frames.append(dict(phi=phi_ocean, lon=lon, lat=lat, name=nom.name))
    del SLC, k1, k2, k3, rgb, alpha

# mosaic into one strip
w = max(f['phi'].shape[1] for f in frames)
def pad(a): return np.pad(a, ((0,0),(0,w-a.shape[1])), constant_values=np.nan)
phi_strip = np.vstack([pad(f['phi']) for f in frames])
lon_strip = np.vstack([pad(f['lon']) for f in frames])
lat_strip = np.vstack([pad(f['lat']) for f in frames])

v = phi_strip[np.isfinite(phi_strip)]
z = np.mean(np.exp(1j*v)); R = float(np.abs(z))
cmean = float(np.angle(z))
sigma = float(np.sqrt(max(0, -2*np.log(max(R, 1e-12)))))

fig, ax = plt.subplots(figsize=(6, 10))
ax.imshow(phi_strip, cmap='gist_rainbow', vmin=cmean-sigma, vmax=cmean+sigma, aspect='auto')
ax.set_title(f'{LABEL}: mosaic strip, {len(frames)} frame(s)')
plt.tight_layout()
plt.show()

# compare to SMOS
smos_file = SMOS_DIR / f'smos_{SCENE_DATE}.nc'
if smos_file.exists():
    with nc4.Dataset(smos_file) as ds:
        key = ('Sea_Surface_Salinity_Rain_Corrected'
               if 'Sea_Surface_Salinity_Rain_Corrected' in ds.variables
               else 'Sea_Surface_Salinity')
        sss = np.squeeze(np.ma.filled(ds.variables[key][:], np.nan)).astype(float)
        la = np.ma.filled(ds.variables['latitude'][:], np.nan).astype(float)
        lo = np.ma.filled(ds.variables['longitude'][:], np.nan).astype(float)
    sss = np.where((sss > 0) & (sss < 42), sss, np.nan)
    LA, LO = np.meshgrid(la, lo, indexing='ij') if sss.ndim == 2 and la.ndim == 1 else (la, lo)
    ok = np.isfinite(sss)
    pts = np.column_stack([LO[ok].ravel(), LA[ok].ravel()])
    vals = sss[ok].ravel()
    sss_on_strip = griddata(pts, vals, (lon_strip, lat_strip), method='linear')

    fig, ax = plt.subplots(1, 2, figsize=(10, 10))
    ax[0].imshow(phi_strip, cmap='gist_rainbow', vmin=cmean-sigma, vmax=cmean+sigma, aspect='auto')
    ax[0].set_title('HH-VV phase strip')
    ax[1].imshow(sss_on_strip, cmap='viridis', aspect='auto')
    ax[1].set_title('SMOS SSS on same swath')
    plt.tight_layout()
    plt.show()
else:
    print(f'no SMOS file at {smos_file}, skipping comparison')

print('done.')
