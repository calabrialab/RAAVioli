#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys 
# acquire ID from file as list, adding the prefix @
ids, lineno, flag, exclude = set('@' + x.rstrip() for x in open(sys.argv[1])), 0, False, True if len(sys.argv)==3 else False
for l in sys.stdin:
    lineno += 1
    if lineno%4 == 1: flag = (l.split(' ')[0].rstrip() in ids) ^ exclude
    if flag: print(l, end='')
