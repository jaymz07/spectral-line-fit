#!/usr/bin/env python

"""
  plot_atlas
  Plot line position vs strength in 'atlas' style.

  usage:
  plot_atlas [options] files

  -h          help
  -c char     comment character(s) used in files (default '#', several characters allowed)
  -I          vertical impulses (bars) to indicate line strength (default: use '+')
  -g int      (max.) number of graphs per page with one molecule per graph
              (default: plot all lines (of all molecules) in one graph)
  -o          output (print) file (default: interactive xmgr plot)
  -v          verbose

  NOTE:
  the current implementation uses the xmgr plotting tool (a.k.a. ACE/gr, further developed into the GRACE tool).
  If you prefer the new xmgrace, use the -v option (verbose) to print the command to the screen and replace xmgr -> xmgrace.
  (Unfortunately the option syntax is slightly different, but the "all in one" plot should work as before.)
"""

import os
import sys
from types import *
from string import *

catsPath = os.path.split(sys.path[0])[0]
sys.path.insert (0, os.path.join(catsPath,'aux'))
from pairTypes import Interval

verbose=0

####################################################################################################################################
def aceGr_atlas (lineFiles, printFile='', nGraphs=0, symbol=9, title=''):
	# postscript file to be produced
        if printFile:
		aceGr =  "grbatch -hardcopy -printfile " + printFile
		if printFile.endswith('eps'): aceGr =  aceGr + ' -eps'
		print 'plotting in batch mode, postscript file: ', printFile
	else:
		aceGr = "xmgr"
	# title, logarithmic scale, x and y axis titles
	aceGr = aceGr + """ -pexec 'xaxis label "position [cm\S-1\N]"'"""
	if nGraphs==0:
        	if title:    aceGr = aceGr + """ -pexec 'title "%s"'""" % title
		# one plot/graph for all molecules
		aceGr = aceGr + """ -log y  -legend load  -pexec 'yaxis label "S [cm\S-1\N / (molec.cm\S-2\N)]"'"""
		for nFile,file in enumerate(lineFiles):
			set = ' -pexec "s%i linestyle 0" -pexec "s%i symbol %i" %s' % (nFile, nFile, symbol, file)
			aceGr = aceGr + set
		if verbose: print aceGr
		os.system(aceGr)
	else:
		numFiles = len(lineFiles)
		if nGraphs>numFiles or nGraphs<0: nGraphs = numFiles
		lineFiles = lineFiles[::-1]
		for i in range(len(lineFiles)/nGraphs+1):
			Files = lineFiles[nGraphs*i:nGraphs*(i+1)]
			print Files, [os.path.splitext(file)[0] for file in Files]
			nFiles = len(Files)
			graphs = ' -rows ' + str(nFiles)
			if nFiles>10: graphs = graphs + ' -maxgraph %i' % nFiles
			elif nFiles<1: break
			set   = ' -pexec "s0 linestyle 0" -pexec "s0 symbol %i"  -pexec "s0 color 2" ' % symbol
			for nFile,file in enumerate(Files):
				set = set + """ -pexec 'subtitle "%s"'""" % os.path.splitext(file)[0]
				if nFile>0:
					set = set + """ -pexec 'xaxis ticklabel off' """
				if nFile==nFiles-1:
					set = set + """ -pexec 'yaxis label "S [cm\S-1\N / (molec.cm\S-2\N)]"' """
					if title: set = set + """ -pexec 'title "%s"'""" % title
				graphs = graphs + " -graph %i -log y  %s %s" % (nFile, set, file)
			grCommand = aceGr+graphs
        		if printFile and nGraphs>1:
			   	he = os.path.splitext(printFile)
				grCommand = replace(grCommand, printFile, '%s.%i%s' % (he[0], i, he[1]))
			if verbose: print '\n', grCommand
			os.system(grCommand)
	
####################################################################################################################################

if (__name__ == "__main__"):
	from IO import open_outFile
	from command_parser import parse_command, standardOptions

        opts = standardOptions + [
	       {'ID': 'v'},
	       {'ID': 'p'},
	       {'ID': 'I'},
	       {'ID': 'g', 'name': 'nGraphs', 'type': IntType, 'default': 0},
	       {'ID': 't', 'name': 'title', 'type': StringType}
               ]

	Files, options, commentChar, outFile = parse_command (opts,(1,99))

	if options.has_key('h'):
		print __doc__%globals();  raise SystemExit, " end of plot_atlas help"

	verbose = options.has_key('v')
	if options.has_key('I'): symbolNumber=12
	else:                    symbolNumber= 9

	aceGr_atlas (Files, outFile, options.get('nGraphs',0), symbolNumber, options.get('title'))
