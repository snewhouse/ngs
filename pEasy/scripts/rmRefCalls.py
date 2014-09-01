#!/usr/bin/env python

'''removes variant calls with no alt allele (refcalls)'''

import sys

for line in sys.stdin:
    if not line.startswith('#'):
        f = line.split()
        if f[4]=='.':
            continue
    sys.stdout.write(line)