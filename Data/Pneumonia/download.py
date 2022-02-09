#!/usr/bin/env python3

from re import search
from urllib.request import urlopen
from os.path import exists
from os import remove
import gzip
from contextlib import closing
from shutil import copyfileobj


def download(url, filename):
    with closing(urlopen(url)) as r:
        with open(filename, 'wb') as f:
            copyfileobj(r, f)

def extract_gz(in_file, out_file):
    with gzip.open(in_file, 'rb') as f_in:
        with open(out_file, 'wb') as f_out:
            copyfileobj(f_in, f_out)

def main():
    with open('GSE164948_series_matrix.txt') as geo_file:
        for line in geo_file.readlines():
            if "!Series_supplementary_file" in line:
                url = search("(?P<url>ftp?://[^\s]+)", line).group("url").replace("\"", "")
                print(url)
                x = url.split("/")
                gz_file = x[len(x)-1]
                filename = gz_file.replace(".gz", "")
                print(filename)
                if not exists(filename):
                    download(url, gz_file)
                    extract_gz(gz_file, filename)
                    remove(gz_file)

                print("-------------------")

if __name__ == '__main__':
    main()
