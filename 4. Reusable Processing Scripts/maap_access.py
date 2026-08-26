# auth, catalogue search, and downloads for the ESA MAAP STAC catalogue - I have cleaned this up with Co-Pilot as it was messy and hard to follow
import time
import requests
from pathlib import Path
from pystac_client import Client

TOKEN_URL = 'https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token'
CATALOG_URL = 'https://catalog.maap.eo.esa.int/catalogue/'
SECRETS_PATH = r'C:\Users\Orlaith.Doyle\STEP_BACK\secrets.txt'

DEFAULT_ASSETS = [
    'enclosure_i_abs_tiff', 'enclosure_i_phase_tiff',
    'enclosure_annot_xml', 'enclosure_orb_xml', 'enclosure_kml',
]

def load_credentials(secrets_path=SECRETS_PATH):
    ns = {}
    exec(open(secrets_path).read(), ns)
    return ns['maap_client_id'], ns['maap_client_secret'], ns['maap_offline_token']

def get_token(secrets_path=SECRETS_PATH):
    client_id, client_secret, offline_token = load_credentials(secrets_path)
    r = requests.post(TOKEN_URL, data={
        'client_id': client_id,
        'client_secret': client_secret,  # needed or you get a 401 even with a valid refresh token
        'grant_type': 'refresh_token',
        'refresh_token': offline_token,
        'scope': 'offline_access openid'}, timeout=60)
    r.raise_for_status()
    return r.json()['access_token']

def open_catalog(token, catalog_url=CATALOG_URL):
    return Client.open(catalog_url, headers={'Authorization': f'Bearer {token}'})

def search_retry(bbox, date=None, collection='BiomassLevel1a', tries=5,
                  secrets_path=SECRETS_PATH):
    # re-auths on every retry since a stale token is the usual failure mode
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
    # writes to a .part file first so a failed download never leaves a
    # corrupt file sitting at the real destination
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
            token = get_token()
    return False

def download_granule(item, dest_dir, token, assets=DEFAULT_ASSETS):
    # skips files already downloaded, never deletes anything
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
