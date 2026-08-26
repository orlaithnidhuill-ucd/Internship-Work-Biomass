"""
maap_access.py
Authentication, catalogue search (with retry), and asset download for the
ESA MAAP STAC catalogue (explorer.maap.eo.esa.int).

Consolidated from repeated cells across the El Nino case-study notebook
(catalogue-scan cell, wide re-scan cell, and download cell).

Credentials are loaded from C:\\Users\\Orlaith.Doyle\\STEP_BACK\\secrets.txt
via exec(), which must define: maap_client_id, maap_client_secret,
maap_offline_token. This module never stores or prints those values.

Usage:
    from maap_access import get_token, open_catalog, search_retry, download

    token = get_token()
    cat = open_catalog(token)
    items = search_retry(cat, bbox=[-170, -5, -120, 5], date='2026-06-30')
    download(items[0].assets['enclosure_i_abs_tiff'].href, dest_path, token)
"""
import time
import requests
from pathlib import Path
from pystac_client import Client

TOKEN_URL = 'https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token'
CATALOG_URL = 'https://catalog.maap.eo.esa.int/catalogue/'
SECRETS_PATH = r'C:\Users\Orlaith.Doyle\STEP_BACK\secrets.txt'

# Asset keys typically wanted for a full Biomass L1a granule
DEFAULT_ASSETS = [
    'enclosure_i_abs_tiff', 'enclosure_i_phase_tiff',
    'enclosure_annot_xml', 'enclosure_orb_xml', 'enclosure_kml',
]


def load_credentials(secrets_path=SECRETS_PATH):
    """Executes the secrets file into a private namespace and returns the
    three MAAP credential variables. Never prints or logs their values."""
    ns = {}
    exec(open(secrets_path).read(), ns)
    return ns['maap_client_id'], ns['maap_client_secret'], ns['maap_offline_token']


def get_token(secrets_path=SECRETS_PATH):
    """Fetches a fresh MAAP access token via refresh-token grant."""
    client_id, client_secret, offline_token = load_credentials(secrets_path)
    r = requests.post(TOKEN_URL, data={
        'client_id': client_id,
        'client_secret': client_secret,  # required — omitting causes 401 even with a valid refresh token
        'grant_type': 'refresh_token',
        'refresh_token': offline_token,
        'scope': 'offline_access openid'}, timeout=60)
    r.raise_for_status()
    return r.json()['access_token']


def open_catalog(token, catalog_url=CATALOG_URL):
    """Opens a pystac_client Client authenticated with the given token."""
    return Client.open(catalog_url, headers={'Authorization': f'Bearer {token}'})


def search_retry(bbox, date=None, collection='BiomassLevel1a', tries=5,
                  secrets_path=SECRETS_PATH):
    """
    Searches the MAAP STAC catalogue with automatic token refresh and
    exponential backoff on failure. Re-authenticates on every retry since
    the most common failure mode is a stale/invalid token.

    Parameters
    ----------
    bbox : [min_lon, min_lat, max_lon, max_lat]
    date : 'YYYY-MM-DD' or None for no date filter
    tries : max attempts before raising

    Returns list of pystac Items.
    """
    token = get_token(secrets_path)
    cat = open_catalog(token)
    kwargs = {'collections': [collection], 'bbox': bbox}
    if date:
        kwargs['datetime'] = f'{date}/{date}'
    for k in range(tries):
        try:
            return list(cat.search(**kwargs).items())
        except Exception as e:
            if k == tries - 1:
                raise
            print(f'  retry {k + 1}: {str(e)[:50]}')
            time.sleep(4 * (k + 1))
            token = get_token(secrets_path)
            cat = open_catalog(token)


def download(url, dest, token, tries=4):
    """
    Streams a single asset to disk with retry, exponential backoff, token
    refresh on failure, and a size check against Content-Length before the
    temp file is renamed into place. Never overwrites a partial file with
    a corrupt one — writes to `dest.name + '.part'` first.

    Returns True on success, False if all retries were exhausted.
    """
    dest = Path(dest)
    for k in range(tries):
        try:
            with requests.get(url, headers={'Authorization': f'Bearer {token}'},
                               stream=True, timeout=1800) as r:
                r.raise_for_status()
                expected = int(r.headers.get('content-length', 0))
                tmp = dest.with_suffix(dest.suffix + '.part')
                with open(tmp, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8 << 20):
                        f.write(chunk)
            if expected and tmp.stat().st_size != expected:
                tmp.unlink()
                raise IOError('size mismatch')
            tmp.rename(dest)
            return True
        except Exception as e:
            if k == tries - 1:
                print(f'      GAVE UP: {str(e)[:55]}')
                return False
            time.sleep(15 * (2 ** k))
            token = get_token()  # refresh before next attempt
    return False


def download_granule(item, dest_dir, token, assets=DEFAULT_ASSETS):
    """
    Downloads all requested assets for one STAC item into dest_dir/<item.id>/.
    Skips files that already exist and are non-empty. Never deletes anything.
    """
    nom = Path(dest_dir) / item.id
    nom.mkdir(parents=True, exist_ok=True)
    for key in assets:
        if key not in item.assets:
            continue
        url = item.assets[key].href
        dest = nom / url.split('/')[-1]
        if dest.exists() and dest.stat().st_size > 0:
            continue
        download(url, dest, token)
    return nom
