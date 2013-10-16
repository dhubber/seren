#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import glob
import time
import psutil

# Set parameters
batch_min = 10         # Minimum number of arguments to xargs (per batch)
batch_max = 500        # Maximum number of arguments to xargs (per batch)
safety_time = 120      # Do not remove files newer than this (in seconds)
memory_mult = 1.2      # Estimated memory usage is file size multiplied by this
default_overwrite = "y"   # Overwrite existing files by default?
default_remove_old = "y"  # Remove old files by default?
write_start_message = "y" # Display parameters at start
      # Do not set default_overwrite="n" and default_remove_old="y"
      # as this is silly and will delete files without using them
sizes = {'':1,'b':1,'kb':1024,'mb':1048576,'gb':1073741824}
chosen_unit = (sizes['mb'],'MB')

def available_ram():
# Read available physical RAM and available RAM using psutil
    vmem_info = psutil.virtual_memory()
    total_ram = getattr(vmem_info,'total')
    available_ram = getattr(vmem_info,'available')
    return (total_ram, available_ram)

def numtasks(needed_mem, warn=False):
    needed_mem = int(needed_mem)
    total_ram, remaining_ram = available_ram()
    if warn and needed_mem > remaining_ram:
        print ("Warning: required ram %i %s greater than available ram %i %s" %
            (needed_mem/chosen_unit[0],chosen_unit[1],
            remaining_ram/chosen_unit[0],chosen_unit[1]))
    n = int(remaining_ram / needed_mem)
    return n
    
if (default_overwrite=='n') and (default_remove_old=='y'):
    print "Settings are silly"
    sys.exit(-4)
    
# Current epoch time:
epoch_time = time.time()
#print "Epoch time is %i" % (epoch_time)

# Find out the number of cores in the system
numprocs = os.sysconf('SC_NPROCESSORS_ONLN')

print "%i cores available" % (numprocs)

ram_now = available_ram()
total_ram = ram_now[0]
ram_now = [x/chosen_unit[0] for x in ram_now]
print "Currently available RAM: %i %s of %i %s" % (ram_now[1], chosen_unit[1],
                                                   ram_now[0], chosen_unit[1])

#print "%i command line arguments" % (len(sys.argv))

# Check we have command line arguments (stub name, [overwrite])
if len(sys.argv) == 1:
    print "No command arguments!"
    sys.exit(-2)

stubname = sys.argv[1]

# Check for overwrite command line argument
if len(sys.argv) == 3:
    if sys.argv[2] == "y":
        overwrite = "y"
        remove_old = "n"
    elif sys.argv[2] == "yes":
        overwrite = "y"
        remove_old = "n"
    elif sys.argv[2] == "r":
        overwrite = "y"
        remove_old = "y"
    elif sys.argv[2] == "rm":
        overwrite = "y"
        remove_old = "y"
    elif sys.argv[2] == "n":
        overwrite = "n"
        remove_old = "n"
    elif sys.argv[2] == "no":
        overwrite = "n"
        remove_old = "n"
    else:
        print "Unknown command argument (argument 2)!"
        sys.exit(-2)
else:
    overwrite = default_overwrite
    remove_old = default_remove_old

if (write_start_message=="y"):
    print "Current parameters: "
    if (overwrite=="y"):
        print "Overwrite existing files"
    else:
        print "Do not overwrite existing files"
    if (remove_old=="y"):
        print "Remove old files after joining (unless modified within last %i seconds)" % (safety_time)
    else:
        print "Do not remove old files after joining"
    print "Batch large numbers of files in batches of no more than %i" % (batch_max)

#print "stubname is %s" % (stubname)

# Construct search path from stub name
testpath = stubname + ".*.MPI.0"
testpath_debug = stubname + ".debug.*.0"

#print "testpath is %s" % (testpath)
#print "testpath_debug is %s" % (testpath)

# Produce file list of MPI files to join
filelist = glob.glob(testpath)
filelist = filter(lambda x: (x.find("debug") == -1),filelist)
filelist_debug = glob.glob(testpath_debug)

if len(filelist) == 0:
    if len(filelist_debug) == 0:
        print "Filelist is empty!"
        sys.exit(-3)
#else:
#   print filelist

numfiles = len(filelist)
numfiles_debug = len(filelist_debug)

#print "numfiles = %i" % (numfiles)
#print "numfiles_debug = %i" % (numfiles_debug)

strippedname = []
strippedname_debug = []
filenames_debug = []
numthreads = []
flag = []
max_out_filesize = 0

# Loop over all files to create, find number of threads to join and check it is not already done
for a in range(numfiles):
    #print "a = %i %s %s" % (a, filelist[a], filelist[a].replace(".MPI.0",""))
    testpath = filelist[a].replace(".MPI.0","")
    #print "new test path = %s" % (testpath)
    #print os.path.isfile(testpath)
    if (overwrite=="y") or (not os.path.isfile(testpath)):
        strippedname.append(testpath)
        n_threads = len(glob.glob(testpath+".MPI.[0-9]*"))
        numthreads.append(n_threads)
        if (remove_old=="y"):
            curflag = "r"
        else:
            curflag = overwrite
        #print "n_threads = %i" % (n_threads)
        out_filesize = 0
        for t in range(n_threads):
            t_filename = testpath+".MPI."+str(t)
            #print "Thread number %i, file %s" % (t, t_filename)
            t_filesize = os.path.getsize(t_filename)
            #print "File size %i" % (t_filesize)
            out_filesize = out_filesize + t_filesize
            # File modification time check
            if (remove_old=="y"):
                #print "Testing file %s" % (t_filename)
                file_age = epoch_time - os.path.getmtime(t_filename)
                #print "File age = %i" % (file_age)
                if (file_age < safety_time):
                    curflag = "y"
        #print "Output file size %i" % (out_filesize)
        max_out_filesize = max(max_out_filesize,out_filesize)
        flag.append(curflag)

for a in range(numfiles_debug):
    testpath = filelist_debug[a].replace(".MPI.0","")
    if (overwrite=="y") or (not os.path.isfile(testpath)):
        strippedname_debug.append(testpath)
        globstring = ""
        globlist = glob.glob(testpath+".MPI.[0-9]*")
        for item in globlist:
            globstring = globstring + " " + item
        filenames_debug.append(globstring)

#print "Maximum output file size %i" % (max_out_filesize)

#print strippedname

numfiles = len(strippedname)
numfiles_debug = len(strippedname_debug)

if (numfiles > 0):
# Work out how many of the maximum file size would fit in memory
    maxjoins = int(total_ram / (memory_mult*max_out_filesize))
#  print "Maximum number of concurrent joins %i" % (maxjoins)
    if maxjoins == 0:
        print "Insufficient total memory; will try one task at a time"
        max_numtasks = 1
    elif maxjoins < numprocs:
        max_numtasks = maxjoins
        print ("Using at most %i concurrent joins " % (max_numtasks) +
           "because of limited physical memory")
    else:
        max_numtasks = numprocs
else:
    max_numtasks = numprocs

# Construct command line to fire off xargs
def command(n):
    command_string = "'| xargs -t -P " + str(n) + \
        " -L 1 joinsims" # took out an -x
    return command_string
#   print "Command is: %s" % (command)
#   print "Argumentlist is: %s" % (argumentlist)

# xargs for output files

# Guess a good batch size
batch_size = max_numtasks * 3
batch_size = max(batch_size,batch_min)
batch_size = min(batch_size,batch_max)

if (numfiles == 0):
    print "No files to create"
    if numfiles_debug == 0:
        print "No debug files to create"
        sys.exit(0)

elif numfiles <= batch_size:
    argumentlist = "echo '"

#  print "numfiles is %i" % (numfiles)

    print "%i files to create..." % (numfiles)

    print "%i debug files to create..." % (numfiles_debug)

# Loop over files to create and construct argumentlist to pass to xargs
    for a in range(numfiles):
        curname = strippedname[a]
        n = numthreads[a]
        curflag = flag[a]
        #print "Adding to argument list: %s %i %s" % (curname, n, curname)
        addtolist = curname + " " + str(n) + " " + curname + " " + curflag
        #print "joinsims %s" % (addtolist)
        argumentlist = argumentlist + addtolist + "\n"

# Run xargs to run joinsims on all necessary files
    ntasks = min(max_numtasks,numtasks(max_out_filesize*memory_mult,True))
    ntasks = max(ntasks,1)
    os.system(argumentlist + command(ntasks))
    #print argumentlist + command(ntasks)

else:
    argumentlist = "echo '"
    print "Running in batches of %i" % (batch_size)
    batch_start = 0
    batch_end = batch_size
    while batch_start < numfiles:
        argumentlist = "echo '"
        # Do xargs in batches
        if (batch_end > numfiles):
            batch_end = numfiles
        for a in range(batch_start,batch_end):
            curname = strippedname[a]
            n = numthreads[a]
            curflag = flag[a]
            addtolist = curname + " " + str(n) + " " + curname + " " + curflag
            argumentlist = argumentlist + addtolist + "\n"
        # Run xargs to run joinsims on all necessary files
        ntasks = min(max_numtasks,numtasks(max_out_filesize*memory_mult))
        ntasks = max(ntasks,1)
        os.system(argumentlist + command(ntasks))
        #print argumentlist + command(ntasks)
        # Increment counters
        batch_start = batch_start + batch_size
        batch_end = batch_start + batch_size - 1

# Construct debug command line to fire off xargs
if (numfiles_debug > 0):
    command_debug = "'| xargs -t -P " + str(numprocs) + " -L 1 joindebugs.bash" # took out an -x
#   print "Command is: %s" % (command)
#   print "Argumentlist is: %s" % (argumentlist_debug)

if (numfiles_debug == 0):
    print "No debug files to create"
elif numfiles_debug <= batch_max:
    argumentlist_debug = "echo '"
    for a in range(numfiles_debug):
        curname = strippedname_debug[a]
        addtolist = curname + " " + filenames_debug[a]
        argumentlist_debug = argumentlist_debug + addtolist + "\n"

# Run xargs to run joindebugs.bash on all necessary files
    os.system(argumentlist_debug + command_debug)
else:
    print "Running in batches of %i" % (batch_max)
    batch_start = 0
    batch_end = batch_max
    while batch_start < numfiles_debug:
        argumentlist_debug = "echo '"
        # Do xargs in batches
        if (batch_end > numfiles_debug):
            batch_end = numfiles_debug
        for a in range(batch_start,batch_end):
            curname = strippedname_debug[a]
            addtolist = curname + " " + filenames_debug[a]
            argumentlist_debug = argumentlist_debug + addtolist + "\n"
        # Run xargs to run joinsims on all necessary files
        os.system(argumentlist_debug + command_debug)
        # Increment counters
        batch_start = batch_start + batch_max
        batch_end = batch_start + batch_max

print "Job completed; %i files created" % (numfiles+numfiles_debug)

