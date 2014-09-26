#!/usr/bin/env python

import sys
import os
import time
import logging
import functools
import subprocess
#import matplotlib.pyplot as plt

'''get GIT revision hash'''
def githash():
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).rstrip()

'''timer (decorator)'''
def timejob(stream):
    def actual_time_job(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Start the timer.
            start = time.time()
            # Run the decorated function.
            ret = func(*args, **kwargs)
            # Stop the timer.
            end = time.time()
            elapsed = end - start
            if isinstance(stream,logging.Logger) or isinstance(stream,logging.RootLogger):
                stream.debug("[TIME] {} took {} seconds".format(func.__name__, elapsed))
            else:
                stream.write("[TIME] {} took {} seconds".format(func.__name__, elapsed))
        return wrapper
    return actual_time_job


'''parses picard metrics files into an dictionary (cols) of lists (rows)'''
def parsePicard(fh):
    allcols = {}
    allnames = {}
    cols = {}
    indexToName = {}
    lineNum = 0
    for line in fh:
        if line.startswith('#'):
            pass  # comment line
        elif len(line)>1:
            if lineNum == 0:  # is header
                headings = line.split()
                for i, heading in enumerate(headings):
                    heading = heading.strip()
                    cols[heading] = []
                    indexToName[i] = heading
            else:  # is data
                cells = line.split()
                for i, cell in enumerate(cells):
                    cell = cell.strip()
                    cols[indexToName[i]] += [cell]
            lineNum += 1
        else: # separator (empty line)
            # append
            for i, n in indexToName.items():
                allnames[i+len(allcols)] = n
            allcols.update(cols)
            # reset
            cols = {}
            indexToName = {}
            lineNum = 0
    # update allcols
    for i, n in indexToName:
        allnames[i+len(allcols)] = n
    allcols.update(cols)


    # try finding row labels (column with non numeric fields)
    for i in sorted(allnames.keys()):
        try:
            float(allcols[allnames[i]][0])
        except ValueError:
            # is not a number => set as row label
            allrows = allcols[allnames[i]]
            break
    # set default if no labels (only print length of rows as in first column)
    try:
        allrows
    except:
        allrows = ['metric' for x in allcols[allnames[0]]]

    return allcols, allnames, allrows


def plot():
    pass

# returns the tail of a file (buffered traversal)
def tail( f, window=4 ):
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if (bytes - BUFSIZ > 0):
            # Seek back one whole BUFSIZ
            f.seek(block*BUFSIZ, 2)
            # read BUFFER
            data.append(f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.append(f.read(bytes))
        linesFound = data[-1].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return '\n'.join(''.join(data).splitlines()[-window:])

# flattens arbitrarily nested containers (eg list)
def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i

# truncate a file to zero bytes, and preserve its original modification time
def zeroFile(file):
    if os.path.exists(file):
        # save the current time of the file
        timeInfo = os.stat(file)
        try:
            f = open(file,'w')
        except IOError:
            pass
        else:
            f.truncate(0)
            f.close()
            # change the time of the file back to what it was
            os.utime(file,(timeInfo.st_atime, timeInfo.st_mtime))

if __name__=="__main__":
    picardValues, ordering = parsePicard(sys.stdin)
    print picardValues
    print ordering

'''
if __name__=="__main__":
    print >> sys.stderr, " -- This contains only helper functions"
    print >> sys.stderr, " == USE ngsEasy.py!"
'''