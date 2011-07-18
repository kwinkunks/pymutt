'''Run the mtft examples, etc.

usage: python runme.py [options]

    -h         print help info for runme.py
    -v         be more verbose
    -N case    add case to the list of tests to execute

runme with no options prints usage information for pymutt.mtft.

test cases are:
    1    simple randomly-generated fourier line extraction
    2    analysis of a seismograph record
    3    willamette river flow data
'''

import sys
import os
import traceback

import pymutt

import simpletest as t1
import freeosc as t2
import willamettedata as t3

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    import getopt

    if argv is None:
        argv = sys.argv

    options = "hvN:"
    verbose = 0
    cases = []

    try:
        try:
            opts, args = getopt.getopt(argv[1:], options, ["help"])
            for opt, val in opts:
                if opt == '-h':
                    print>> sys.stderr,  __doc__
                    return 0
                elif opt == "-v":
                    verbose += 1
                elif opt == "-N":
                    cases.append(int(val))

        except getopt.error, msg:
            print >> sys.stderr, __doc__
            raise Usage(msg)

        if len(cases) == 0:
            print >> sys.stderr, pymutt.mtft.__doc__
            return 0

        for c in cases:
            print "\n --- Running test %d ---\n" % c
            if c == 1:
                t1.doit()
            elif c == 2:
                t2.doit()
            elif c == 3:
                t3.doit()
            else:
                print >> sys.stderr, "unknown case %s" % c
                return 2
        
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

    except RuntimeError, err:
        print >>sys.stderr, err
        return 3

    except:
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
    
