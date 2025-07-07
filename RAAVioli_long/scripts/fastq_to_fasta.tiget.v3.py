#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# da usare quando vuoi estrarre FASTQ avendo i dati dal BED.

import sys 
# acquire ID from file as list, adding the prefix @
lineno = 0
for l in sys.stdin:
    lineno += 1
    if lineno%4 == 1: 
        l = l.replace("@",">")
        print(l.rstrip())
    elif lineno%4 == 2:
        print(l.rstrip())
