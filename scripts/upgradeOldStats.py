#!/usr/bin/env python2.7
"""
upgradeOldStats.py: condense old position.stats.tsv files and discard read names for correctly-mapped reads.

"""

import argparse, sys, os, os.path, random, itertools, string, re
import doctest

import tsv

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("--input_tsv", type=argparse.FileType("r"),
        default=sys.stdin,
        help="input stats TSV")
    parser.add_argument("--output_tsv", type=argparse.FileType("w"),
        default=sys.stdout,
        help="output stats TSV")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Need to convert from:
    # correct, mapq, condition, read
    # To:
    # correct, mapq, condition, read, count
    
    # Start the output file
    writer = tsv.TsvWriter(options.output_tsv)
    
    # The super easy way is to accumulate reads and basically run-length compress
    accumulating = None
    
    is_header = True
    for parts in tsv.TsvReader(options.input_tsv):
        if is_header:
            # Improve the header
            writer.list_line(parts + ['count'])
            is_header = False
            continue
    
        if accumulating is not None and accumulating[0] == parts[0] and accumulating[1] == parts[1] and accumulating[2] == parts[2]:
            # Merge in
            accumulating[4] += 1
        else:
            if accumulating is not None:
                # Flush the buffer
                writer.list_line(accumulating)
                accumulating = None
                
            if parts[0] == '1':
                # Correct reads are accumulated
            
                # Say we will accumulate on this read
                accumulating = parts
                # Clear the name
                accumulating[3] = '.'
                # Tack on a count
                accumulating.append(1)
            else:
                # wrong reads are just dumped
                writer.list_line(parts + [1])
                
    
    if accumulating is not None:
        # Flush the buffer
        writer.list_line(accumulating)
    
    return 0
        
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

