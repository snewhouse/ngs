#!/usr/bin/env python

import sys
import re

for line in sys.stdin:
    m = re.match('##(\w+)=<(.*)>', line)
    if m:
        fields = m.group(2).split(',')
        try:
            data = { k:v for (k, v) in [ f.split('=') for f in fields ]}
        except ValueError:
            print line
            sys.exit()
        else:
            print data