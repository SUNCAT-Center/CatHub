#!/usr/bin/env python
import os
import ase.db

from cathub.postgresql import CathubPostgreSQL

def main():
    server_name = CathubPostgreSQL().server_name
    db = ase.db.connect(server_name)
    return db
