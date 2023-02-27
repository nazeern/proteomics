"""
Read a csv file with peptide sequences in a column called "sequence" and append
the requested EM(p) fits for every row.

For the infile, we only need a "sequence" column. We don't care what the other
columns are -- we automatically copy all of them.

Run this file with "-h" to see options.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Homepage: http://caseywstark.com
License:
  Copyright (C) 2011, 2012 Casey W. Stark. All Rights Reserved.

  This file is part of the MIDA-Kinemed pipeline.

"""

### imports
# python stdlib
import argparse
import copy
import csv
import os
import sys
import time

start_time = time.time()

# other libs
import numpy as np
import scipy
import scipy.optimize
from scipy.stats import pearsonr

# mida
from mida.analysis import convert_p_to_abundances, renormalize
from mida import Peptide

# local
from run_data import chemical_data, enriched_aa_abundances, enriched_aa_fractions

###
# Fitting functions
###

# We assume the y-intercept is 0 for these functions
quad = lambda p, x: p[0] * x**2 + p[1] * x
quad_error = lambda p, x, y: quad(p, x) - y

cubic = lambda p, x: p[0] * x**3 + p[1] * x**2 + p[2] * x
cubic_error = lambda p, x, y: cubic(p, x) - y

###
# Command line options
###
parser = argparse.ArgumentParser(description="The MIDA EM(p) fitter. Appends the requested data to the input csv file.")

parser.add_argument("infile", type=argparse.FileType("r"))
parser.add_argument("outfile", type=argparse.FileType("wb"))
parser.add_argument("-e", default=4,
                    dest="em_max", type=int,
                    help="Expects integer. The EM(p) curve to fit up to. Ex: -e 4 outputs EM0 - EM4 fits.")
parser.add_argument("-f", default="cubic",
                    dest="fit", type=str,
                    help="Expects the fit name. 'quad' or 'cubic'")
parser.add_argument("-s", default=0,
                    dest="sequence_format", type=int,
                    help="Expects integer. '0' for normal, '1' for (a)sequence(b)")

args = parser.parse_args()

# before we do anything, handle bad input
if args.em_max > 4:
    raise Exception("Max EM fit error. We can only fit up to EM4, got EM%i." % args.em_max)
if args.fit != "quad" and args.fit != "cubic":
    raise Exception("Fit name error. Expected quad or cubic, got %s." % args.fit)
if args.sequence_format != 0 and args.sequence_format != 1:
    raise Exception("Sequence format error. The only formats are '0' and '1', got %i." % arg.sequence_format)

###
# Define the fields we will add.
###
quad_fields = ["EM_i_ quad coeff 2", "EM_i_ quad coeff 1", "EM_i_ pearson r"]
cubic_fields = ["EM_i_ cubic coeff 3", "EM_i_ cubic coeff 2", "EM_i_ cubic coeff 1", "EM_i_ pearson r"]

fields_to_add = ["composition", "mass", "n", "M0", "M1", "M2", "M3", "M4"]

if args.fit == "quad":
    for i in range(args.em_max+1):
        for field in quad_fields:
            fields_to_add.append( field.replace("_i_", str(i)) )
else:
    for i in range(args.em_max+1):
        for field in cubic_fields:
            fields_to_add.append( field.replace("_i_", str(i)) )

###
# Open files and setup reader and writer
###
if args.sequence_format == 0:
    reader = csv.DictReader(args.infile)
else:
    reader = csv.DictReader(args.infile, delimiter=';')

writer = csv.writer(args.outfile)

###
# Intro text
###
print ""
print "MIDA-Kinemed EM(p) Fitter"
print "========================="
print ""

###
# Tell the user what we will be doing
###
ems_string = "EM_0(p)"
if args.em_max >= 1:
    ems_string += ", EM_1(p)"
if args.em_max >= 2:
    ems_string += ", EM_2(p)"
if args.em_max >= 3:
    ems_string += ", EM_3(p)"
if args.em_max == 4:
    ems_string += ", EM_4(p)"

fields_string = ""
for field in fields_to_add:
    fields_string += "%s, " % field
fields_string = fields_string[:-2]

print "Using the sequences in '%s', in format %i." % (args.infile.name, args.sequence_format)
print "Performing %s fits to %s." % (args.fit, ems_string)
print ""
print "Adding the columns:"
print fields_string
print ""

sys.stdout.write("Opening %s output and writing headers ..." % args.outfile.name)

###
# Figure out the headers, write them to the outfile
###
# build outfile headers
in_fields = copy.copy(reader.fieldnames)
# get rid of empty columns...
in_fields = filter(None, in_fields)
out_fields = copy.copy(in_fields)
out_fields.extend(fields_to_add)
writer.writerow(out_fields)

sys.stdout.write(" done.\n")
print ""
print "Generating EM(p) and computing fits for:"

###
# Loop through the input file rows, generate the isotopomer distribution vs. p
# for the peptide, make the EM(p) data, fit it, and write the results to the
# output file.
###
# Some variables used inside the loop -- construct here so we only have to do
# it once.
# the enrichment resolution
num_ps = 50
p_array = np.linspace(0.0, 0.05, num=num_ps)
exp_abundances = np.zeros((len(p_array), 5))

h_abundances = convert_p_to_abundances(p_array, chemical_data.natural_abundances[0])

# to keep track of how many peptides we have processed
peptide_count = 0

for row in reader:
    ###
    # Grab the peptide sequence. Fail gracefully if there is a problem reading.
    ###
    if args.sequence_format == 0:
        try:
            seq = row["sequence"]
        except KeyError:
            print "Error: could not read the column 'sequence' in this file."
            print "Are you using the correct sequence format? Try format '1'."
            print ""
            sys.exit()
    else:
        try:
            # the '1' format is a bit more complicated. The headers have spaces
            # after the ';'. We would set the delimiter to '; ' to avoid this,
            # but just ';' is the delimiter in the rest of the file.
            seq = row[" sequence"]
            # Now we have to get rid of the (x) at the beginning and end of the
            # sequence string. This could go in the above line, but this is
            # clearer.
            seq = seq[3:-3]  # just slice the first and last 3 chars out
        except KeyError:
            print "Error: could not read the column 'sequence' in this file."
            print "Are you using the correct sequence format? Try format '0'."
            print ""
            sys.exit()

    # Construct the peptide given in this row
    peptide = Peptide(seq, chemical_data)

    # update the peptide count before we print which number we are on
    peptide_count += 1

    sys.stdout.write("(%i) %s ... " % (peptide_count, seq))
    sys.stdout.flush()

    ###
    # Compute isotopomer distribution vs. overenrichment.
    ###
    abundances = peptide.get_distribution(labile_abundances=(h_abundances,),
                     en_aa_abundances=enriched_aa_abundances,
                     en_aa_fraction=enriched_aa_fractions, mass_cutoff=4)

    # experimentally renormalize and add to the list
    for i in xrange(len(abundances)):
        exp_abundances[i] = renormalize(peptide, abundances[i])

    ###
    # Generate the EM(p) data
    ###
    # 0-th element is the distribution at p=0 (natural abundances)
    ems_array = np.zeros((5, num_ps))
    ems_array[0] = exp_abundances[:, 0] - exp_abundances[0, 0]
    ems_array[1] = exp_abundances[:, 1] - exp_abundances[0, 1]
    ems_array[2] = exp_abundances[:, 2] - exp_abundances[0, 2]
    ems_array[3] = exp_abundances[:, 3] - exp_abundances[0, 3]
    ems_array[4] = exp_abundances[:, 4] - exp_abundances[0, 4]

    ###
    # Create the list that we will write
    ###
    write_row = []
    for field in in_fields:
        write_row.append(row[field])

    # append new stuff
    write_row.append(peptide.formula)
    write_row.append(peptide.base_mass)
    write_row.append(peptide.labile_groups[0].num_atoms)  # number of labile H sites
    write_row.append(abundances[0][0])  # fractional abundance of M0 isotopomer at natural abundances
    write_row.append(abundances[0][1])  # M1 isotopomer at natural abundances
    write_row.append(abundances[0][2])  # M2 isotopomer at natural abundances
    write_row.append(abundances[0][3])
    write_row.append(abundances[0][4])

    ###
    # Fit the EM(p) data with whichever curve the user asked for.
    ###
    if args.fit == "quad":
        for i in range(args.em_max+1):
            # perform the fit
            q_coeffs, success = scipy.optimize.leastsq(quad_error, [1, 1],
                                    args=(p_array, ems_array[i]))
            # compute the correlation of the fit
            r = pearsonr(ems_array[i], quad(q_coeffs, p_array))
            # add the fit coeffs and correlation to the list
            write_row.extend([q_coeffs[0], q_coeffs[1], r])

    else:
        for i in range(args.em_max+1):
            # perform the fit
            c_coeffs, success = scipy.optimize.leastsq(cubic_error, [1, 1, 1],
                                    args=(p_array, ems_array[i]))
            # compute the correlation of the fit
            r = pearsonr(ems_array[i], cubic(c_coeffs, p_array))
            # add the fit coeffs and correlation to the list
            write_row.extend([c_coeffs[0], c_coeffs[1], c_coeffs[2], r])

    ###
    # Write the current row list
    ###
    writer.writerow(write_row)

    # let the user know this row is done before starting the next.
    print "done."

    #break  # for faster testing

end_time = time.time()  # cheap profiling

total_time = end_time - start_time
rate = peptide_count / total_time

print ""
print "Run time: %f" % total_time
print "Rate: %f per second" % rate
print ""
print "Done with all rows in `%s`." % args.infile.name
print "The output with fits is in `%s`" % args.outfile.name
print "Have a nice day."
print ""
