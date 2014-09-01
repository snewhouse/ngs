#!/usr/bin/awk -f
( $5!="." && $5!="N" || /^#/ ) { print $0 }

