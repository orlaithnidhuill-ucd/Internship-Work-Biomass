# Demerara / Guiana Shelf
import re
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import cdsapi
from pathlib import Path
from datetime import datetime
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
exec(open(BASE / 'secrets.txt').read())   # maap_client_id, maap_client_secret, maap_offline_token, cds_url, cds_key

ARCHIVE   = BASE / 'PLUME_ARCHIVE' / 'demerara'
SAVE_DIR  = ARCHIVE / 'figures'
SMOS_DIR  = ARCHIVE / 'smos'
WC_CACHE  = BASE / 'worldcover_cache'
ERA5_DIR  = ARCHIVE / 'era5'
for d in (SAVE_DIR, SMOS_DIR, ERA5_DIR):
    d.mkdir(parents=True, exist_ok=True)

LABEL = 'Demerara / Guiana Shelf'
BBOX = [-59.0, 6.0, -57.0, 8.5]# [min_lon, min_lat, max_lon, max_lat]
SCENE_DATE = '2026-07-06'  # YYYY-MM-DD
MULTILOOK  = (14, 1)# (azimuth, range) looks;  exactly 14:1
WC_WATER_CLASSES = {0, 80}     # WorldCover: 0 = nodata (open ocean), 80 = permanent water

BAND_ORDER = ('VV', 'VH', 'HV', 'HH')   # TIFF band order, opposite to XML label order

# 1: download or check files are already on system 
def get_token():
    import requests
    r = requests.post('https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token',
                       data={'client_id': maap_client_id, 'client_secret': maap_client_secret,
                             'grant_type': 'refresh_token', 'refresh_token': maap_offline_token,
                             'scope': 'offline_access openid'}, timeout=60)
    r.raise_for_status()
    return r.json()['access_token']

def find_local():
    return sorted(g for g in ARCHIVE.iterdir() if g.is_dir()
                  and next(g.glob('*i_abs.tiff'), None)
                  and next(g.glob('*i_phase.tiff'), None)
                  and next(g.glob('*.kml'), None))

local = find_local()
if local:
    print(f'{len(local)} granule(s) already on system, skipping download')
    nom = local[0]
else:
    print('nothing local... searching catalogue')
    from pystac_client import Client
    token = get_token()
    cat = Client.open('https://catalog.maap.eo.esa.int/catalogue/',
                      headers={'Authorization': f'Bearer {token}'})
    items = list(cat.search(collections=['BiomassLevel1a'], bbox=BBOX,
                            datetime=f'{SCENE_DATE}/{SCENE_DATE}').items())
    if not items:
        raise RuntimeError('no granules found for this site/date')
    it = items[0]
    nom = ARCHIVE / it.id
    nom.mkdir(parents=True, exist_ok=True)
    for key in ['enclosure_i_abs_tiff', 'enclosure_i_phase_tiff', 'enclosure_annot_xml',
                'enclosure_orb_xml', 'enclosure_kml']:
        if key in it.assets:
            import requests
            url = it.assets[key].href
            dest = nom / url.split('/')[-1]
            if not dest.exists():
                with requests.get(url, headers={'Authorization': f'Bearer {token}'},
                                  stream=True, timeout=1800) as resp:
                    with open(dest, 'wb') as f:
                        for chunk in resp.iter_content(chunk_size=8 << 20):
                            f.write(chunk)
    print(f'downloaded: {nom.name}')

print(f'using granule: {nom.name}')

# read the 4 SLC channels
def read_bands(path):
    if BACKEND == 'rasterio':
        with rasterio.open(path) as src:
            return {BAND_ORDER[b]: src.read(b + 1) for b in range(4)}
    arr = tifffile.imread(path)
    return {BAND_ORDER[b]: arr[b] for b in range(4)}

amp = read_bands(next(nom.glob('*i_abs.tiff')))
pha = read_bands(next(nom.glob('*i_phase.tiff')))
SLC = {p: amp[p] * np.exp(1j * pha[p]) for p in BAND_ORDER}
del amp, pha

# 2. Pauli, H-alpha, 4 polarisations
def pauli_components(hh, hv, vv):
    return hh + vv, hh - vv, 2 * hv  # k1 (surface), k2 (dihedral), k3 (volume)
def multilook_amplitude(slc, looksa, looksr):
    return np.sqrt(uniform_filter(np.abs(slc) ** 2, [looksa, looksr]))
def pauli_rgb(k1, k2, k3, looksa, looksr, clip_sigma=2.5):
    amp1 = multilook_amplitude(k1, looksa, looksr)
    amp2 = multilook_amplitude(k2, looksa, looksr)
    amp3 = multilook_amplitude(k3, looksa, looksr)
    rgb = np.stack([amp2, amp3, amp1], axis=-1).astype('float32')  # R=k2 G=k3 B=k1
    for b in range(3):
        rgb[:, :, b] = np.clip(rgb[:, :, b], 0, clip_sigma * np.mean(rgb[:, :, b]))
        rgb[:, :, b] /= np.max(rgb[:, :, b])
    return rgb

def alpha_angle(k1, k2, k3):
    # single-pixel quick estimate, no eigen-decomposition
    plen = np.sqrt(np.abs(k1)**2 + np.abs(k2)**2 + np.abs(k3)**2)
    return np.arccos(np.abs(k1) / plen)

k1, k2, k3 = pauli_components(SLC['HH'], SLC['HV'], SLC['VV'])
rgb = pauli_rgb(k1, k2, k3, *MULTILOOK)
alpha = alpha_angle(k1, k2, k3)

fig, ax = plt.subplots(1, 2, figsize=(14, 6))
ax[0].imshow(rgb, aspect='auto')
ax[0].set_title(f'{LABEL}: Pauli RGB (R=HH-VV, G=2HV, B=HH+VV)')
im = ax[1].imshow(np.degrees(alpha), cmap='viridis', aspect='auto', vmin=0, vmax=90)
ax[1].set_title('alpha angle (deg)')
plt.colorbar(im, ax=ax[1], fraction=0.03)
plt.tight_layout()
plt.savefig(SAVE_DIR / 'pauli_alpha.png', dpi=170, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(1, 4, figsize=(18, 5))
for a, pol in zip(ax, BAND_ORDER):
    im = a.imshow(20*np.log10(np.abs(SLC[pol]) + 1e-6), cmap='gray', aspect='auto')
    a.set_title(pol)
plt.tight_layout()
plt.savefig(SAVE_DIR / 'four_pols.png', dpi=170, bbox_inches='tight')
plt.show()

# 3. Multilooking @ 14:1
print(f'multilook window: {MULTILOOK} (ratio {MULTILOOK[0]/MULTILOOK[1]:.0f}:1)')
print(f'native shape {SLC["HH"].shape} -> multilooked '
      f'{SLC["HH"].shape[0]//MULTILOOK[0]} x {SLC["HH"].shape[1]//MULTILOOK[1]} (approx)')


# 4. land masking (geocorrected) + geolocation
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

n_lines, n_samples = SLC['HH'].shape
lon, lat = corner_grid(nom, n_lines, n_samples)

# WorldCover-based ocean mask (class 0 (nodata, mostly open ocean) counts as water)
wc_file = next(WC_CACHE.glob(f'*{LABEL.split()[0]}*.tif'), None)
if wc_file is not None:
    with rasterio.open(wc_file) as src:
        wc_arr = src.read(1)
    ocean_mask = np.isin(wc_arr, list(WC_WATER_CLASSES))
    ocean_mask = binary_dilation(ocean_mask, iterations=2)  # buffer the coastline
else:
    print('no cached WorldCover tile for this site -- masking skipped (all pixels treated as valid)')
    ocean_mask = np.ones_like(lon, dtype=bool)

print(f'ocean fraction: {ocean_mask.mean()*100:.1f}%')


# HH-VV phase
phi = np.angle(SLC['HH'] * np.conj(SLC['VV']))
phi_masked = np.where(ocean_mask, phi, np.nan)
v = phi_masked[np.isfinite(phi_masked)]
z = np.mean(np.exp(1j*v)); R = float(np.abs(z))
cmean = float(np.angle(z))
sigma = float(np.sqrt(max(0, -2*np.log(max(R, 1e-12)))))

fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(phi_masked, cmap='gist_rainbow', vmin=cmean-sigma, vmax=cmean+sigma, aspect='auto')
plt.colorbar(im, ax=ax, label='HH-VV phase (rad)')
ax.set_title(f'{LABEL}: HH-VV phase, ocean only  (mean {np.degrees(cmean):+.2f} deg)')
plt.tight_layout()
plt.savefig(SAVE_DIR / 'hhvv_phase.png', dpi=170, bbox_inches='tight')
plt.show()


# ERA5 wind correction
scene_date = datetime.strptime(SCENE_DATE, '%Y-%m-%d')
scene_lat, scene_lon = float(np.nanmean(lat)), float(np.nanmean(lon))
era5_file = ERA5_DIR / f'era5_wind_{SCENE_DATE.replace("-","")}.nc'

if era5_file.exists() and era5_file.stat().st_size < 1024:
    era5_file.unlink()

if not era5_file.exists():
    print('downloading ERA5 wind...')
    c = cdsapi.Client(url=cds_url, key=cds_key)
    try:
        c.retrieve('reanalysis-era5-single-levels', {
            'product_type': 'reanalysis',
            'variable': ['10m_u_component_of_wind', '10m_v_component_of_wind'],
            'year': str(scene_date.year), 'month': f'{scene_date.month:02d}',
            'day': f'{scene_date.day:02d}',
            'time': ['00:00','03:00','06:00','09:00','12:00','15:00','18:00','21:00'],
            'area': [scene_lat+2, scene_lon-2, scene_lat-2, scene_lon+2],
            'data_format': 'netcdf',
        }, str(era5_file))
    except Exception as e:
        print(f'ERA5 download failed: {e} continuing without wind correction')
        era5_file = None

U10 = None
if era5_file is not None and era5_file.exists():
    with nc4.Dataset(era5_file) as ds:
        times = ds.variables['valid_time'][:]
        u_all = ds.variables['u10'][:]
        v_all = ds.variables['v10'][:]
    overpass_epoch = float(np.datetime64(f'{SCENE_DATE}T00:40:00', 's').astype(float))
    nearest = int(np.argmin(np.abs(np.array(times, dtype=float) - overpass_epoch)))
    U10 = float(np.nanmean(np.sqrt(u_all[nearest]**2 + v_all[nearest]**2)))
    print(f'U10 at nearest step: {U10:.2f} m/s')
else:
    print('no wind data available, roughness correction skipped')

# 7. SSS validation, SMOS first, SMAP as a plan B
def load_smos(date):
    f = SMOS_DIR / f'smos_{date}.nc'
    if not f.exists():
        return None
    with nc4.Dataset(f) as ds:
        key = ('Sea_Surface_Salinity_Rain_Corrected'
               if 'Sea_Surface_Salinity_Rain_Corrected' in ds.variables
               else 'Sea_Surface_Salinity')
        sss = np.squeeze(np.ma.filled(ds.variables[key][:], np.nan)).astype(float)
        la = np.ma.filled(ds.variables['latitude'][:], np.nan).astype(float)
        lo = np.ma.filled(ds.variables['longitude'][:], np.nan).astype(float)
    return sss, la, lo

def load_smap(date):
    f = SMOS_DIR / f'smap_{date}.nc'
    if not f.exists():
        return None
    with nc4.Dataset(f) as ds:
        sss = np.squeeze(np.ma.filled(ds.variables['smap_sss'][:], np.nan)).astype(float)
        la = np.ma.filled(ds.variables['latitude'][:], np.nan).astype(float)
        lo = np.ma.filled(ds.variables['longitude'][:], np.nan).astype(float)
    return sss, la, lo

ref = load_smos(SCENE_DATE)
source = 'SMOS'
if ref is None:
    print('SMOS unavailable, falling back to SMAP')
    ref = load_smap(SCENE_DATE)
    source = 'SMAP'

if ref is not None:
    sss, la, lo = ref
    sss = np.where((sss > 0) & (sss < 42), sss, np.nan)
    LA, LO = np.meshgrid(la, lo, indexing='ij') if sss.ndim == 2 and la.ndim == 1 else (la, lo)
    ok = np.isfinite(sss)
    pts = np.column_stack([LO[ok].ravel(), LA[ok].ravel()])
    vals = sss[ok].ravel()
    sss_on_swath = griddata(pts, vals, (lon, lat), method='linear')
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    ax[0].imshow(phi_masked, cmap='gist_rainbow', vmin=cmean-sigma, vmax=cmean+sigma, aspect='auto')
    ax[0].set_title('HH-VV phase')
    im = ax[1].imshow(sss_on_swath, cmap='viridis', aspect='auto')
    ax[1].set_title(f'{source} SSS resampled onto swath')
    plt.colorbar(im, ax=ax[1], fraction=0.03)
    plt.tight_layout()
    plt.savefig(SAVE_DIR / 'sss_validation.png', dpi=170, bbox_inches='tight')
    plt.show()
else:
    print('no SSS reference available (neither SMOS nor SMAP found)')
    sss_on_swath = None


# summary, all together
fig, ax = plt.subplots(2, 2, figsize=(14, 10))
ax[0,0].imshow(rgb, aspect='auto'); ax[0,0].set_title('Pauli RGB')
ax[0,1].imshow(np.degrees(alpha), cmap='viridis', aspect='auto'); ax[0,1].set_title('alpha angle')
ax[1,0].imshow(phi_masked, cmap='gist_rainbow', vmin=cmean-sigma, vmax=cmean+sigma, aspect='auto')
ax[1,0].set_title('HH-VV phase')
if sss_on_swath is not None:
    ax[1,1].imshow(sss_on_swath, cmap='viridis', aspect='auto')
    ax[1,1].set_title(f'{source} SSS')
else:
    ax[1,1].axis('off')
fig.suptitle(f'{LABEL} -- {SCENE_DATE}'
             + (f'  |  U10 {U10:.1f} m/s' if U10 is not None else ''), fontsize=13)
plt.tight_layout()
plt.savefig(SAVE_DIR / 'summary.png', dpi=170, bbox_inches='tight')
plt.show()
print('done.')
