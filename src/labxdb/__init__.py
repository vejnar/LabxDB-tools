# -*- coding: utf-8 -*-

#
# Copyright (C) 2018-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import json
import os

import requests
import urllib3

def get_http_args(http_url=None, http_login=None, http_password=None, http_path=None, db=None):
    if http_url is None:
        if 'LABXDB_HTTP_URL' in os.environ:
            http_url = os.environ['LABXDB_HTTP_URL']
        else:
            http_url = 'http://127.0.0.1:8080'
    if http_login is None and 'LABXDB_HTTP_LOGIN' in os.environ:
        http_login = os.environ['LABXDB_HTTP_LOGIN']
    if http_password is None and 'LABXDB_HTTP_PASSWORD' in os.environ:
        http_password = os.environ['LABXDB_HTTP_PASSWORD']
    if http_path is None:
        if db is None:
            hpath = 'LABXDB_HTTP_PATH'
        else:
            hpath = 'LABXDB_HTTP_PATH_' + db.upper()
        if hpath in os.environ:
            http_path = os.environ[hpath]
    return http_url, http_login, http_password, http_path

def urljoin(*args):
    """
    Joins URLs
    """
    trailing_slash = '/' if args[-1].endswith('/') else ''
    return '/'.join(map(lambda x: str(x).strip('/'), args)) + trailing_slash

class DBLink(object):
    """
    Connecting class to HTTP API returning JSON data implemented with urllib3 and requests.

    Args:
        http_url (str): URL to connect
        http_login (str): Login
        http_password (str): Password
        http_path (str): Path joined to http_url
        db (str): Database name

    Attributes:
        http_url (str): URL to connect
        http_login (str): Login
        http_password (str): Password
        session (Session): Session
    """
    def __init__(self, http_url=None, http_login=None, http_password=None, http_path=None, db=None):
        self.http_url, self.http_login, self.http_password, http_path = get_http_args(http_url, http_login, http_password, http_path, db)
        if http_path is not None:
            self.http_url = urljoin(self.http_url, http_path)
        self.session = requests.Session()
        self.session.mount('http://', requests.adapters.HTTPAdapter(max_retries=urllib3.Retry(total=10, backoff_factor=0.3)))
        self.session.mount('https://', requests.adapters.HTTPAdapter(max_retries=urllib3.Retry(total=10, backoff_factor=0.3)))

    def get(self, url):
        r = self.session.get(urljoin(self.http_url, url), auth=(self.http_login, self.http_password))
        r.raise_for_status()
        if r.headers['Query-status'] != 'OK':
            raise ValueError(r.headers['Query-status'])
        if 'application/json' in r.headers['Content-Type']:
            return r.json()
        else:
            return r.text

    def post(self, url, data=None, json=None):
        r = self.session.post(urljoin(self.http_url, url), data=data, json=json, auth=(self.http_login, self.http_password))
        r.raise_for_status()
        if r.headers['Query-status'] != 'OK':
            raise ValueError(r.headers['Query-status'])
        if 'application/json' in r.headers['Content-Type']:
            return r.json()
        else:
            return r.text
