#!/usr/bin/env python

'''REPLACE.PY

Recursively searches in directory BASE_DIR for all files ending in SUFFIX.
Replaces all occurences of OLD_STR in these files with NEW_STR.

Copyright Adam Stewart, 2014'''


import fnmatch
import glob
import os
import subprocess


BASE_DIR = '/Users/bmonger/python_programs/class_problems_idl'      # base directory to start searching for files from
SUFFIX   = '*.pro'  # suffix of files you want to edit
OLD_STR  = '	'   # string or character you want to get rid of
NEW_STR  = ' '      # string or character to replace OLD_STR with


all_fnames = []
for root, dirnames, filenames in os.walk(BASE_DIR):
    for filename in fnmatch.filter(filenames, SUFFIX):
        all_fnames.append(os.path.join(root, filename))

cmd = ''
for fname in all_fnames:
    cmd += "sed -i '' 's/" + OLD_STR + "/" + NEW_STR + "/g' " + fname + "\n"

subprocess.call(cmd, shell=True)
print '\nAll "' + SUFFIX + '" files in directory "' + BASE_DIR + '" have been fixed\n'
