"""
"""


import os
import gzip
import requests
from urllib import request
import shutil
from contextlib import closing


def ftp_download(url, dir):
    """Downloads using ftp protocol"""
    filename = url.split('/')[-1]
    with closing(request.urlopen(url)) as r:
        with open(dir + filename, 'wb+') as f:
            shutil.copyfileobj(r, f)
    logging.info("Downloaded data from" + url + " stored at " + dir + " as " + filename)
    return dir + filename


def ungzip(path):
    """Decompresses and removes gzip file in same directory."""
    unzip_path = path.split(".")[0]
    with gzip.open(path, 'rb') as f_in:
        with open(unzip_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(path)
    return unzip_path


def get_genome(url, dir):
    """Downloads genome from URL and processes: unzip +/- index generation."""
    zip_path = ftp_download(url, dir)
    unzip_path = ungzip(zip_path)
    return unzip_path


def get_gff(url, dir):
    """Downloads gff from URL"""
    zip_path = ftp_download(url, dir)
    unzip_path = ungzip(zip_path)
    return unzip_path


def get_cpgs(url, dir) -> str:
    """Downloads platfrom cpg info from URL"""
    r = requests.get(url)
    if r.status_code == 200:
        filename = r.headers['Content-Disposition'].split("=")[-1]
        with open(dir + filename, 'w+') as f:
            f.write(r.text)
    logging.info("Downloaded Illumina 450k cpg metadata from " + url + " stored at " + dir + " as " + filename)
    return dir + filename