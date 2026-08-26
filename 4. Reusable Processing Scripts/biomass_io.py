# Code cell for reading Biomass granules, TIFF bands, granule ID parsing, KML georeferencing
import re
import numpy as np
from datetime import datetime
from lxml import etree

BAND_ORDER = ('VV', 'VH', 'HV', 'HH')  # TIFF band order (opposite to XML label order)
try:
    import rasterio
    BACKEND = 'rasterio'
except Exception:
    import tifffile
    BACKEND = 'tifffile'

def read_bands(path):
    if BACKEND == 'rasterio':
        with rasterio.open(path) as src:
            return {BAND_ORDER[b]: src.read(b + 1).astype(np.float32) for b in range(4)}
    arr = tifffile.imread(path)
    return {BAND_ORDER[b]: arr[b].astype(np.float32) for b in range(4)}

def parse_granule_id(gid):
    dates = re.findall(r'(\d{8}T\d{6})', gid)
    start = datetime.strptime(dates[0], '%Y%m%dT%H%M%S') if dates else None
    trk = re.search(r'_(T\d{3})_', gid)
    frm = re.search(r'_(F\d{3})_', gid)
    return (start,
            int(trk.group(1)[1:]) if trk else None,
            int(frm.group(1)[1:]) if frm else None)

def corner_grid(granule_dir, n_lines, n_samples):
    # bilinear interp between the 4 KML corners, oriented using orbitPass
    kml = next(granule_dir.glob('*.kml'), None)
    if kml is None:
        return None, None
    m = re.search(r'<coordinates>(.*?)</coordinates>', kml.read_text(), re.S)
    if not m:
        return None, None
    pts = np.array([tuple(map(float, p.split(',')))
                    for p in m.group(1).split()])[:4, :2]

    edge_len = lambda a, b: np.hypot(*(pts[b] - pts[a]))
    if (edge_len(0, 1) + edge_len(2, 3)) < (edge_len(1, 2) + edge_len(3, 0)):
        row_groups, col_groups = ({0, 1}, {2, 3}), ({1, 2}, {0, 3})
    else:
        row_groups, col_groups = ({1, 2}, {3, 0}), ({0, 1}, {2, 3})

    row_lats = [(g, np.mean([pts[i][1] for i in g])) for g in row_groups]
    south = min(row_lats, key=lambda x: x[1])[0]
    north = max(row_lats, key=lambda x: x[1])[0]

    ann_path = next(granule_dir.glob('*annot.xml'), None)
    orbit_pass = 'Ascending'
    if ann_path is not None:
        ann = etree.parse(str(ann_path))
        orbit_pass = next((e.text.strip() for e in ann.iter()
                            if etree.QName(e).localname == 'orbitPass' and e.text),
                           'Ascending')
    row0, rowM = (north, south) if orbit_pass == 'Descending' else (south, north)

    col_lons = [(g, np.mean([pts[i][0] for i in g])) for g in col_groups]
    col0 = max(col_lons, key=lambda x: x[1])[0]
    colM = min(col_lons, key=lambda x: x[1])[0]

    A = pts[list(row0 & col0)[0]]
    B = pts[list(row0 & colM)[0]]
    C = pts[list(rowM & col0)[0]]
    D = pts[list(rowM & colM)[0]]

    frac_r = np.linspace(0, 1, n_lines)[:, None, None]
    frac_c = np.linspace(0, 1, n_samples)[None, :, None]
    top = A * (1 - frac_c) + B * frac_c
    bot = C * (1 - frac_c) + D * frac_c
    grid = top * (1 - frac_r) + bot * frac_r
    return grid[..., 0], grid[..., 1]

def hhvv_phase(amp, pha):
    hh = amp['HH'] * np.exp(1j * pha['HH'])
    vv = amp['VV'] * np.exp(1j * pha['VV'])
    return np.angle(np.exp(1j * (np.angle(hh) - np.angle(vv))))


def is_complete_granule(granule_dir):
    return bool(
        next(granule_dir.glob('*i_abs.tiff'), None)
        and next(granule_dir.glob('*i_phase.tiff'), None)
        and next(granule_dir.glob('*.kml'), None)
    )
