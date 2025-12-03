#!/usr/bin/env python
import os
from .config import server_name

def main(user):
    os.system('psql --host={} --port=5432 --username={} --dbname=catalysishub --password'.format(server_name, user))


if __name__ == '__main__':
    from sys import argv
    user = argv[1]
    main(user)
