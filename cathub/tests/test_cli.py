
import os
import sys
import unittest
import tempfile
import pprint
import sqlite3
import json
import click
from click.testing import CliRunner

import cathub


class CommandLineTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_imports(self):
        import cathub.cathubsqlite
        import cathub.postgresql
        import cathub.folderreader
        import cathub.tools
        import cathub.ase_tools
        import cathub.folder2db
        import cathub.db2server
        import cathub.make_folders_template
        import cathub.psql_server_connect
        import cathub.organize

    def test_cli_query(self):
        runner = CliRunner()
        from cathub.cli import reactions, publications
        runner.invoke(reactions)
        runner.invoke(publications)

    def test_cli_make_folders(self):
        from cathub.cli import make_folders
        runner = CliRunner()
        runner.invoke(make_folders, ['--create-template', 'template'])
        runner.invoke(make_folders, ['template'])

    def test1_cli_read_folders(self):
        from cathub.cli import folder2db
        runner = CliRunner()
        runner.invoke(folder2db, ['aayush/'])

    def test2_cli_db2server(self):
        from cathub.postgresql import CathubPostgreSQL
        from cathub.cli import db2server
        db = CathubPostgreSQL(user='postgres')
        con = db._connect()
        db._initialize(con)
        db.truncate_schema()
        runner = CliRunner()
        runner.invoke(db2server, ['--dbuser=postgres',
                                  'aayush/MontoyaChallenge2015.db'])
    def test3_cli_asedb(self):
        from cathub.cli import ase
        runner = CliRunner()
        runner.invoke(ase, ['--dbuser=postgres', '--dbpassword=None'])

    def test4_show_reactions(self):
        from cathub.cli import show_reactions
        runner = CliRunner()
        runner.invoke(show_reactions, ['aayush/MontoyaChallenge2015.db'])

    def test_reactions(self):
        from cathub.cli import reactions
        runner = CliRunner()
        runner.invoke(reactions, ['-q chemicalComposition=~Co', '-q reactants=H'])

    def test_publications(self):
        from cathub.cli import publications
        runner = CliRunner()
        runner.invoke(publications)

if __name__ == '__main__':
    unittest.main()
