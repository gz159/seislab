#!/usr/bin/env python

from sys import exit, argv, stderr

if len(argv) == 2:
   mw = argv[1]
elif len(argv) == 1:
    try:
        mw = raw_input('Enter Mw : ')
    except NameError:
        stderr.write('Syntax: %s PZ_filename\n'%argv[0])
        exit(1)
else:
   stderr.write('Syntax: %s PZ_filename\n'%argv[0])
   exit(1)
   
hd=1.2*pow(10.,-8.)*pow(pow(10.,1.5*float(mw)+16.1),1./3)

print "The half during is:"+str(hd)
