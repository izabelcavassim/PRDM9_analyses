import time
import logging

from io import StringIO

import skbio as sk
import pandas as pd
from sys import argv
import requests

#from .utilities import dump_fasta_from_pairs
__all__ = ('cdd')


logger = logging.getLogger(__name__)


BASE_URL = 'https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi'


class CDDError(Exception):
    ERROR_MESSAGES = {
        1: 'Invalid search ID.',
        2: 'No effective input (usually no query proteins or search ID specified).',
        4: 'Queue manager (qman) service error.',
        5: 'Data is corrupted or no longer available (cache cleaned, etc.).'
    }

    @staticmethod
    def from_status_code(error_code):
        return CDDError(CDDError.ERROR_MESSAGES[error_code])


def _parse_cdd(output):
    df = pd.read_csv(StringIO(output), sep='\t', skiprows=7)
    df['Query'] = df['Query'].str.extract('Q#([0-9]+).*').astype(int)
    df_sub = df.loc[df["Hit type"] == "superfamily"]
    return df_sub


def cdd(queries, smode='auto', maxhit=500, tdata='hits', dmode='all', qdefl=True, cddefl=True, useid1=True, evalue=0.1):
    #queries = dump_fasta_from_pairs(ids, sequences)


    params = {
        'queries': queries,
        'smode': smode,
        'maxhit': maxhit,
        'evalue': evalue,
        'tdata': tdata,
        'dmode': dmode,
        'qdefl': 'false' if not qdefl else 'true',
        'cddefl': 'false' if not cddefl else 'true'
    }

    resp = requests.post(BASE_URL, data=params)
    resp.raise_for_status()

    lines = resp.text.splitlines()
    cdsid = lines[1].split()[1]

    logger.debug('Submitted query to CDD with params. Received id %s.', cdsid)

    status = 100
    while status > 0:
        resp = requests.post(BASE_URL, params={'cdsid': cdsid})
        lines = resp.text.splitlines()
        print(lines)

        status = int(lines[3].split()[1])
        logger.debug('Query %s has status %d', cdsid, status)

        if status == 0:
            print(status)
            logger.debug('Query completed!')
            print(_parse_cdd(resp.text))
            results = _parse_cdd(resp.text)
            return(results)

        if status in (1, 2, 4, 5):
            raise CDDError.from_status_code(status)

    else:
        logger.debug('Query in progress. Waiting.')
        #time.sleep(10)
        print('Query in progress. Waiting.')