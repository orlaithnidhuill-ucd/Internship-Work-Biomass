"""
granule_cache.py
Batch-loads a directory of Biomass granules into an in-memory `layers` list
(one dict per frame: track, frame, phi, lon, lat) plus a summary DataFrame,
with pickle caching so a kernel restart doesn't force a full ~50-minute
reload, and incremental extension when new granules are added later.

Generalizes the load / cache-load / extend-with-new-frames cells that were
duplicated (with only path constants changed) across the El Nino, Amazon,
and Ganges notebooks.

Usage:
    from pathlib import Path
    from granule_cache import load_or_build_cache, extend_cache

    layers, df = load_or_build_cache(
        archive_dir=Path(r'C:\\Users\\Orlaith.Doyle\\Nino\\Nino-prods\\2026-06-30'),
        cache_path=Path(r'C:\\Users\\Orlaith.Doyle\\Nino\\layers_2026-06-30.pkl'),
        decimation=6, min_frames=3)

    # later, after downloading more granules into the same archive_dir:
    layers, df = extend_cache(layers, df, archive_dir=..., cache_path=..., decimation=6)
"""
import gc
import pickle
import numpy as np
import pandas as pd
from pathlib import Path

from biomass_io import read_bands, parse_granule_id, corner_grid, hhvv_phase, is_complete_granule


def _process_one(nom, decimation):
    """Loads one granule, returns (layer_dict, row_dict) or raises."""
    _, trk, frm = parse_granule_id(nom.name)
    amp = read_bands(next(nom.glob('*i_abs.tiff')))
    pha = read_bands(next(nom.glob('*i_phase.tiff')))
    nl, ns = amp['HH'].shape
    phi = hhvv_phase(amp, pha)
    del amp, pha

    lon, lat = corner_grid(nom, nl, ns)
    if lon is None:
        raise ValueError('no KML / georeferencing failed')

    pl = phi[::decimation, ::decimation]
    v = pl[np.isfinite(pl)]
    z = np.mean(np.exp(1j * v))

    layer = dict(track=trk, frame=frm, phi=pl,
                 lon=lon[::decimation, ::decimation].copy(),
                 lat=lat[::decimation, ::decimation].copy())
    row = dict(track=trk, frame=frm,
               mean_deg=round(float(np.degrees(np.angle(z))), 3),
               coh=round(float(np.abs(z)), 4),
               lon_c=round(float(np.nanmean(lon)), 2),
               lat_c=round(float(np.nanmean(lat)), 2))
    del lon, lat, phi, pl
    return layer, row


def load_granules(archive_dir, decimation=6, min_frames=3, progress_every=20):
    """
    Loads every complete granule in archive_dir into `layers` (list of
    dicts) and a summary DataFrame `df`. Tracks with fewer than
    min_frames are dropped. Prints progress every `progress_every`
    granules rather than per-granule (quiet mode).
    """
    archive_dir = Path(archive_dir)
    noms = sorted(g for g in archive_dir.iterdir()
                  if g.is_dir() and is_complete_granule(g))
    print(f'{len(noms)} complete granule(s) in {archive_dir}\n')

    layers, rows, failed = [], [], []
    for n, nom in enumerate(noms, 1):
        try:
            layer, row = _process_one(nom, decimation)
            layers.append(layer)
            rows.append(row)
        except Exception as e:
            _, trk, frm = parse_granule_id(nom.name)
            failed.append((trk, frm, f'{type(e).__name__}: {str(e)[:40]}'))
        finally:
            gc.collect()
        if n % progress_every == 0 or n == len(noms):
            print(f'  loaded {n}/{len(noms)}  ({len(layers)} ok, {len(failed)} failed)')

    if not layers:
        raise RuntimeError('nothing loaded')

    df = pd.DataFrame(rows)
    counts = df.groupby('track')['frame'].size()
    keep = sorted(counts[counts >= min_frames].index)
    layers = [l for l in layers if l['track'] in keep]
    df = df[df['track'].isin(keep)].reset_index(drop=True)

    if failed:
        print(f'\n{len(failed)} granule(s) failed:')
        for trk, frm, msg in failed[:10]:
            print(f'  T{trk}/F{frm}: {msg}')

    return layers, df


def load_or_build_cache(archive_dir, cache_path, decimation=6, min_frames=3):
    """
    Loads layers/df from cache_path if it exists, otherwise builds it from
    archive_dir with load_granules() and writes the cache. This is the
    pattern to run first after every kernel restart.
    """
    cache_path = Path(cache_path)
    if cache_path.exists():
        with open(cache_path, 'rb') as f:
            blob = pickle.load(f)
        layers, df = blob['layers'], blob['df']
        print(f'loaded cache: {len(layers)} frames from {cache_path.name}')
        return layers, df

    layers, df = load_granules(archive_dir, decimation=decimation, min_frames=min_frames)
    with open(cache_path, 'wb') as f:
        pickle.dump({'layers': layers, 'df': df}, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f'\ncached to {cache_path.name} ({cache_path.stat().st_size / 1e6:.0f} MB)')
    return layers, df


def extend_cache(layers, df, archive_dir, cache_path, decimation=6):
    """
    Scans archive_dir for granules not already present in `layers`
    (matched on track+frame), processes only those, appends to layers/df,
    and rewrites the cache. Never removes or reprocesses existing frames.
    """
    archive_dir = Path(archive_dir)
    have = {(l['track'], l['frame']) for l in layers}
    noms = sorted(g for g in archive_dir.iterdir()
                  if g.is_dir() and is_complete_granule(g))
    new_noms = [n for n in noms if parse_granule_id(n.name)[1:] not in have]
    print(f'{len(new_noms)} new granule(s) to process')
    if not new_noms:
        return layers, df

    new_rows = []
    for nom in new_noms:
        try:
            layer, row = _process_one(nom, decimation)
            layers.append(layer)
            new_rows.append(row)
            print(f"  T{row['track']:03d}/F{row['frame']:03d}: "
                  f"mean {row['mean_deg']:+.2f}°, coh {row['coh']:.3f}, "
                  f"lat_c {row['lat_c']:+.1f}")
        except Exception as e:
            _, trk, frm = parse_granule_id(nom.name)
            print(f'  T{trk}/F{frm}: FAILED ({e})')
        finally:
            gc.collect()

    if new_rows:
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
        with open(Path(cache_path), 'wb') as f:
            pickle.dump({'layers': layers, 'df': df}, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f'\n{len(new_rows)} frame(s) added, cache updated')

    return layers, df
