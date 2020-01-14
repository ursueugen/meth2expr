"""
"""


import os
import gzip
import requests
from urllib import request
import shutil
from contextlib import closing

from Bio import SeqIO


def ftp_download(url, dir):
    """Downloads using ftp protocol"""
    filename = url.split('/')[-1]
    with closing(request.urlopen(url)) as r:
        with open(dir + filename, 'wb+') as f:
            shutil.copyfileobj(r, f)
    return dir + filename


def ungzip(path):
    """Decompresses and removes gzip file in same directory."""
    unzip_path = "".join(path.split(".")[0]+".fna")
    with gzip.open(path, 'rb') as f_in:
        with open(unzip_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(path)
    return unzip_path


def fasta_header(path, new_path):
    """Edits header by removing other info than chromosome id."""
    with open(path, 'r') as f_in:
        with open(new_path, 'w+') as f_out:
            records = SeqIO.parse(f_in, 'fasta')
            for record in records:
                record.id = record.id.split(" ")[0]
                record.description = record.id.split(" ")[0]
                SeqIO.write(record, f_out, 'fasta')
    return new_path


def get_genome(url, dir):
    """Downloads genome from URL and processes: unzip +/- index generation."""
    zip_path = ftp_download(url, dir)
    unzip_path = ungzip(zip_path)

    #edit fasta headers to be compatible with gff
    edited_path = "".join(unzip_path.split(".")[:-1] + ['_edited'] + ['.fna'])
    edited_path = fasta_header(unzip_path, edited_path)
    return edited_path


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
    return dir + filename