"""
biomass_io.py
Low-level readers for a single Biomass L1a SCS granule: polarimetric TIFF
bands, granule-ID parsing, and KML-corner-based per-pixel georeferencing.

Consolidated from repeated cells across the El Nino case-study notebook
(granule-load cell, cache-load cell, and new-frame-processing cell) — all
three defined near-identical copies of these three functions.

IMPORTANT: TIFF band order is (VV, VH, HV, HH) — opposite to the XML label
order. This is baked into BAND_ORDER below; do not reorder without checking
against the annotation XML.

Usage:
    from biomass_io import read_bands, parse_granule_id, corner_grid, hhvv_phase

    amp = read_bands(granule_dir.glob('*i_abs.tiff').__next__())
    pha = read_bands(granule_dir.glob('*i_phase.tiff').__next__())
    phi = hhvv_phase(amp, pha)
    lon, lat = corner_grid(granule_dir, *phi.shape)
"""
import re
import numpy as np
from datetime import datetime
from lxml import etree

BAND_ORDER = ('VV', 'VH', 'HV', 'HH')  # TIFF band order — opposite to XML label order

try:
    import rasterio
    _BACKEND = 'rasterio'
except Exception:
    import tifffile
    _BACKEND = 'tifffile'


def read_bands(path):
    """
    Reads a 4-band Biomass polarimetric TIFF and returns a dict keyed by
    polarisation: {'HH': array, 'HV': array, 'VH': array, 'VV': array}.
    Uses rasterio if available, falls back to tifffile (needed when
    rasterio's GDAL DLL is broken, e.g. on this Windows/conda setup).
    """
    if _BACKEND == 'rasterio':
        with rasterio.open(path) as src:
            return {BAND_ORDER[b]: src.read(b + 1).astype(np.float32)
                    for b in range(4)}
    arr = tifffile.imread(path)
    return {BAND_ORDER[b]: arr[b].astype(np.float32) for b in range(4)}


def parse_granule_id(gid):
    """
    Extracts (start_datetime, track_int, frame_int) from a Biomass granule
    ID string. Track/frame are returned as plain ints (e.g. 'T033' -> 33),
    not the zero-padded string, so they sort and compare naturally.
    """
    dates = re.findall(r'(\d{8}T\d{6})', gid)
    start = datetime.strptime(dates[0], '%Y%m%dT%H%M%S') if dates else None
    trk = re.search(r'_(T\d{3})_', gid)
    frm = re.search(r'_(F\d{3})_', gid)
    return (start,
            int(trk.group(1)[1:]) if trk else None,
            int(frm.group(1)[1:]) if frm else None)


def corner_grid(granule_dir, n_lines, n_samples):
    """
    Builds per-pixel (lon, lat) arrays of shape (n_lines, n_samples) by
    bilinear interpolation between the four corner coordinates in the
    granule's KML file, oriented using the orbitPass field in the
    annotation XML (ascending vs descending changes which corner is
    "first row").

    Returns (None, None) if no KML is present or it can't be parsed.

    KML gx:LatLonQuad corner order is: bottom-left, bottom-right,
    top-right, top-left (lon, lat pairs) — this function does not assume
    that order and instead derives north/south/near/far edges from the
    points themselves, then uses orbitPass to decide row 0 vs row -1.
    """
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
    """
    Given amplitude and phase band dicts (as returned by read_bands, for
    the *_i_abs.tiff and *_i_phase.tiff products respectively), returns
    the wrapped HH-VV copolar phase difference array.
    """
    hh = amp['HH'] * np.exp(1j * pha['HH'])
    vv = amp['VV'] * np.exp(1j * pha['VV'])
    return np.angle(np.exp(1j * (np.angle(hh) - np.angle(vv))))


def is_complete_granule(granule_dir):
    """True if a granule folder has amplitude, phase, and KML files."""
    return bool(
        next(granule_dir.glob('*i_abs.tiff'), None)
        and next(granule_dir.glob('*i_phase.tiff'), None)
        and next(granule_dir.glob('*.kml'), None)
    )
