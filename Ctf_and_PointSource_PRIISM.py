#!/usr/bin/python3
#########################################################################################
# Python script to perform the relevant steps for ER Deconvolution (Croxford 2021)    	#
# using an aligned tilt series as the main input.				      	#
#										      	#
# This script expects a ".xf" and ".tlt" file with the same root filename as the      	# 
# aligned and unaligned stacks in order to run CTFFIND. IMOD is also required for some	#
# of the functionality. It also depends on the python packages imported at the 		#
# beginning. Finally, this script requires an existing installation IVE/Priism, which 	#
# can be obtained by request from John Sedat.          					#							
#										      	#		
# Example usage: Ctf_and_PointSource_and_Deconvolution.py -st UNALIGNEDSTACK \        	#
#	-ali ALIGNEDSTACK -b BINFACTOR -z ZDIMENSIONUNBINNED 				#
# 											#
#########################################################################################
import mrcfile
from subprocess import call # Allows to call shell functions
import os 		## Allows for some OS commands, to check if files/paths/directories exist, and to remove files
import numpy as np 	## Matlab-like stuff
import sys 		## Allows command line inputs 
import glob
import argparse
from subprocess import Popen, PIPE, STDOUT

parser = argparse.ArgumentParser(description='Deconvolution Parameters')

parser.add_argument('-st', '--unalignedStack', dest='unalignedStack', type=str)
parser.add_argument('-ali', '--alignedStack', dest='alignedStackDW', type=str)
parser.add_argument('-b', '--bin', dest='Bin', type=int, help='Binning factor', default=4)
parser.add_argument('-z', '--zReconstruction', type=int, dest='ReconZ', help='Full unbinned reconstruction Z dimension', default='1500')

args = parser.parse_args()
Bin = args.Bin
ReconZ = args.ReconZ

unalignedStack = args.unalignedStack #sys.argv[1]
alignedStackDW = args.alignedStackDW #sys.argv[2]
etomoRootName = alignedStackDW.split(".")[0:-1][0]
print(unalignedStack, alignedStackDW, etomoRootName)


header = "/data/software/repo/imod/4.10.28-cuda8.0/bin/header" # to not confuse it with the priism version

if os.path.isfile(unalignedStack) == False:
    raise NameError("Input stack was not found. Please check your inputs.")
    
tltfile = etomoRootName + ".tlt"
print("Tltfile name: " + tltfile)

if os.path.isfile(tltfile) == False:
    raise NameError("Tilt angle file (.tlt) was not found. The file is required to be in the same directory as the the aligned stack (second input). Usually, this will be the etomo directory.)")
    
if os.path.isfile(etomoRootName+".xf") == False:
    raise NameError("Transform-file (.xf-file) not found. The file is required to be in the same directory as the the aligned stack (second input). Usually, this will be the etomo directory.)")
    
if len(unalignedStack.split("/")[-1].split(".")) > 2:
    raise NameError("Please avoid dots in the name of the input file.")
    
# Clean up previous files
if os.path.isfile("temp")==True:
    os.remove("temp")
    
if os.path.isfile("ctfparam.log")==True:
    os.remove("ctfparam.log")
    
if os.path.isfile("here.sh")==True:
    os.remove("here.sh")
    
#if os.path.isfile("ctffind_results_summary.log")==True:
#    os.remove("ctffind_results_summary.log")

if os.path.isfile("avg_def.txt") == False:
    # Re-align the stack using the unaligned stack and the xf-file
    
    alignedStackDW = unalignedStack.split("/")[-1].split(".")[0:-1][0] + "_temp.ali"
    
    call(["newstack -in {} -LinearInterpolation -x {}.xf -out {} -TaperAtFill 1,0 -AdjustOrigin -OffsetsInXandY 0.0,0.0".format(unalignedStack,etomoRootName,alignedStackDW)],shell=True)
        
    # Read parameters of the original stack via imod:
        
    etomoRootName = alignedStackDW.split("/")[-1].split(".")[0:-1][0]
    
    
    # Read size of frames, write them into a temporary text-file (only way to read it)
    
    call(["{} -size {} > temp".format(header,alignedStackDW)],shell=True)
    
    with open("temp","r") as tempFile:
        line = tempFile.read()
        
    dimX = int(line.split()[0])
    dimY = int(line.split()[1])
    dimZ = int(line.split()[2])
    
    os.remove("temp")
    
    # Read pixel size
    
    call(["{} -PixelSize {} > temp".format(header,alignedStackDW)],shell=True)
    
    with open("temp","r") as tempFile:
        line = tempFile.read()
        pixSize = line.split()[0]
        
    os.remove("temp")
    
    # Box out boxes along the tilt axis
    # Write input stack into single frames via newstack
    # Create an "list of outputs" for newstack
    
    with open("temp","a") as tempFile:
        tempFile.write(str(dimZ) + "\n")
        for i in range(0,dimZ):
            tempFile.write("{}_{}.mrc\n{}\n".format(etomoRootName,str(i),str(1)))
    
    call(["newstack -in {} -FileOfOutputs temp".format(alignedStackDW)],shell=True)
    os.remove("temp")
    
    
    # For each projection: Box out 8 boxes along the tilt axis
    
    boxSize = 512 
    nbox = 8
    midX = dimX/2
    midY = dimY/2
    
    # Calculate point coordinates in y 
    pointsY = np.linspace(boxSize/2,dimY-boxSize/2,nbox)
    
    # Write out the boxes via imod's "clip -resize"
    for i in range(0,dimZ):
        clipName = "{}_{}".format(etomoRootName,str(i))
        for j in range(0,nbox):
            boxName = "boxed_{}_{}.mrc".format(clipName,str(j))
            call(["clip resize -ix {} -iy {} -iz 1 -cx {} -cy {} -cz 1 -ox {} -oy {} -oz 1 {}.mrc {}".format(str(dimX),str(dimY),str(midX),str(round(pointsY[j])),boxSize,boxSize,clipName,boxName)],shell=True)
            
    # Output: tomo.ali -> boxed_tomo_1_1.mrc (i.e.fist frame, first subbox)
        call(["newstack boxed_{}*.mrc  boxed_{}.st".format(clipName,clipName)],shell=True)
        
        for k in range(0,nbox):
            os.remove("boxed_{}_{}.mrc".format(clipName,str(k)))
        os.remove(clipName + ".mrc")
    
    # Run ctffind4
    if os.path.exists("ctf")==False:   
        os.mkdir("ctf")
    
    
    # Write a sh here-file to allow running a scripted version of ctffind4 for every frame of the input stack
    
    for i in range(0,dimZ):
        if os.path.isfile("here.sh")==True:
            os.remove("here.sh")
        with open("here.sh","a") as herefile:
            herefile.write("#!/bin/sh\n") # shebang
            herefile.write("for input_fn in boxed_{}_{}.st\n".format(etomoRootName,str(i))) # filename (without image number) as input
            herefile.write("do\n")
            herefile.write("/data/Users/Jan/software/ctffind << here.sh\n") # executes ctffind with the input from the generated file
            herefile.write("${input_fn}\n") # specifies input
            herefile.write("yes\n") # is input file a stack?
            herefile.write("1\n")
            herefile.write("./ctf/${input_fn}_ctfcor.mrc\n")  # specifies output
            herefile.write(str(pixSize) + "\n") # pixel size [A]
            herefile.write("300\n") # acceleration voltage
            herefile.write("2.27\n") # aperture
            herefile.write("0.07\n") # amplitude contrast
            herefile.write("512\n") # size of image for fitting
            herefile.write("50\n") # resolution min
            herefile.write("8\n") # resolution max
            herefile.write("30000\n") # defocus max
            herefile.write("130000\n") # defocus min
            herefile.write("1000\n") # defocus step
            herefile.write("no\n") # do you know how much astigmatism is present?
            herefile.write("no\n") # slower, more exhaustive search?
            herefile.write("yes\n") # use restraints on astigmatism something?
            herefile.write("2000\n") # tolerated astigmatism
            herefile.write("no\n") # Find additional phase shifts?
            herefile.write("no\n") # Expert options?
            herefile.write("here.sh\n")
            herefile.write("done")
        call(["sh here.sh"],shell=True)
        os.rename("dbg_average_spectrum_before_conv.mrc",str(i)+"_dbg_average_spectrum_before_conv.mrc")
    
    # Read the defocus values from the files that were created for every frame of the input stack
        
    df1=[]
    df2=[]
    
    for i in range(0,dimZ):
        filepath="./ctf/boxed_{}_{}.st_ctfcor.txt".format(etomoRootName,str(i))
        # Read each log file, append both defocus values to a list
        with open(filepath,"r") as logfile:
            lines = logfile.readlines()
            with open("ctffind_results_summary.log","a") as summary:
                if i == 0:
                    for subline in lines:
                        summary.write(subline)
                else: 
                    summary.write(lines[5])
                    
            df1.append(lines[5].split()[1])
            df2.append(lines[5].split()[2])
    
    # EDIT BY MATT 2020 05 04
    # Read ctffind_results_summary.log , average columns 1 and 2, make sure it ends up in microns (x10000)
    smry  = np.loadtxt("ctffind_results_summary.log") 
    print(smry)
    avg_def = (smry[:, 2] + smry[:, 3]) / -10000
    np.savetxt("avg_def.txt", avg_def)
    
    # Read the tlt-file, and append the tilts to a list to feed the ctfphaseflip defocus file
    
    with open(tltfile,"r") as tlt:
        tilts = tlt.readlines()
            
    for i in range(0,len(tilts)): # Getting rid of whitespace
        tilts[i] = tilts[i].rstrip()
        
    for i in range(0,dimZ):
        with open("ctfparam.log","a") as ctffile:
            meanDefocus = str(0.1*np.mean([float(df1[i]),float(df2[i])]))
            # Write file for ctfparam: Each line contains frame number (twice, start and end, in our case equal), the tilt angle (twice again) and the mean defocus value
            var = "{} {} {} {} {}\n".format(str(i+1),str(i+1),tilts[i],tilts[i],meanDefocus) 
            ctffile.write(var)
    
    
    # Clean up all temporary files that were written out
    
    if os.path.isfile("here.sh") == True:
        os.remove("here.sh")
        
    for boxfile in glob.glob("./boxed_*"):
        os.remove(boxfile)
        
    os.remove(alignedStackDW)
        
    print("CTF determination concluded.")
    
else:
    print("avg_def.txt already exists, skipping CTFFind") 

    # Read size of frames, write them into a temporary text-file (only way to read it)
    
    call(["{} -size {} > temp".format(header,alignedStackDW)],shell=True)
    
    with open("temp","r") as tempFile:
        line = tempFile.read()
        
    dimX = int(line.split()[0])
    dimY = int(line.split()[1])
    dimZ = int(line.split()[2])
    
    os.remove("temp")
    
    # Read pixel size
    
    call(["{} -PixelSize {} > temp".format(header,alignedStackDW)],shell=True)
    
    with open("temp","r") as tempFile:
        line = tempFile.read()
        pixSize = line.split()[0]
        
    os.remove("temp")

#######################################################################################
## Now Make the Point Source
#######################################################################################
# Refresh the alignedStackDW variable
alignedStackDW = args.alignedStackDW
etomoRootName = alignedStackDW.split(".")[0:-1][0]

with open(tltfile,"r") as tlt:
    tilts = tlt.readlines()

for i in range(0,len(tilts)): # Getting rid of whitespace
    tilts[i] = tilts[i].rstrip()

#os.remove("temp")

call(["{} -size {} > temp".format(header,alignedStackDW)],shell=True)

with open("temp","r") as tempFile:
    line = tempFile.read()

dimX = int(line.split()[0])
dimY = int(line.split()[1])
dimZ = int(line.split()[2])

os.remove("temp")

binDimX = int(np.round(dimX / Bin))
binDimY = int(np.round(dimY / Bin))
binDimZ = int(ReconZ/Bin)

zpixSize = pixSize
    
    
tfName = etomoRootName + "_tf.mrc"
if os.path.isfile(tfName) == False:
    
    ## Enter tilt information
    call(["insert_tilts {} {} -tilt_file={}".format(alignedStackDW, alignedStackDW, tltfile)], shell=True)
    
    ## Append Resolution (for binning)
    call(["AppendRes {} {} {}".format(alignedStackDW, Bin, Bin-1)], shell=True)
    
    ## 2D FFT
    
    # The range selection is to match the chosen reconstruction size.
    FFTFile = etomoRootName + "_FFT.mrc"
    call(["FTransform2D {}  {} -center_zero -real_complex_full -same_units -xpad=0 -ypad=0 -x=0:{} -y=0:{}".format(alignedStackDW, FFTFile, binDimX - 1, binDimY - 1)], shell=True)
    
    ## Thresh
    pointSourceName = etomoRootName + '_delta.mrc'
    call(["Threshold {} {} -not_below=6 -result=mask -mode=short".format(FFTFile, pointSourceName)], shell=True)
    #rm tmp.mrc
    
    ## PFocusRamp through Python interface
    ctfName = etomoRootName + '_delta_ctf.mrc'
    ## avgDef = etomoRootName + '_avg_def.txt'
    avgDef = 'avg_def.txt'
    call(["pfocusramp {} {} -amp=0.07 -axis=0  -cs=2.27 -deftxt={} -kv=300 -mtf=1:5.2279:0:50 -op=apply -multires".format(pointSourceName, ctfName, avgDef)], shell=True)
    
    
    ## EWBP
    ewbpName = etomoRootName + '_psf_ctf.xzyt'
    call(["ewbp {} {} -reconxz={}:{} -sizexz={}:{} -iy=0:{} -filter=2 -hdfilt=1 -moderec=2".format(ctfName, ewbpName, binDimX, binDimZ, binDimX, binDimZ, binDimY -1)], shell=True)
    
    ## 3D fourier transform
    tfName = etomoRootName + '_tf.mrc'
    xshift = round(binDimX/2 - 1)
    yshift = round(binDimY/2 - 1)
    zshift = round(binDimZ/2 - 1)
    
    call(["FTransform3D {} {} -same_units -shift={}:{}:{} -xpad=0 -ypad=0 -zpad=0".format(ewbpName, tfName, str(xshift), str(zshift), str(yshift))], shell=True)
    
else:
    print("Transfer Function Exists, skipping to Deconvolution")
    
#######################################################################################
## Reconstruct the Tomogram to be deconvolved
#######################################################################################
## Back projection of the aligned tilt series to generate DC volume
# First generate the dummy bprmn file
hed = np.array(["SEC", "ROT", "GMAG", "TX", "TY", "SMEAN", "SFIT", "SCALE", "BASE", "TILT"])[np.newaxis]
seq = np.array(range(1,(dimZ+1)))[np.newaxis]
xforms = np.tile([0,1,0,0,1,1,1,0], (dimZ, 1))
ttilts = np.array(tilts)[np.newaxis]

oprmFilePrelim = np.concatenate((seq.T, xforms, ttilts.T), axis=1)
oprmFile = np.concatenate((hed, oprmFilePrelim), axis=0)

np.savetxt("{}.bprmMn".format(etomoRootName), oprmFile, delimiter=' ', fmt='%s')
call(["AppendRes {} 2 1".format(alignedStackDW)], shell=True)
# Apply parameters to from mass norm to the aligned tilt series to generate MnAln
prmStack = etomoRootName + "_sub_appl.mrc"
call(["appl_prm {} {} -iprmfile={} -dimxy={}:{} -fullsize={}:{} -iv=0:{}:1 -shxyz=0:0:0 -iref=-1 -imform=2 -pcbase=0.05 -tilt_offset=0 -statfile=none -res=1 -rscale=2 -multires".format(alignedStackDW, prmStack, "{}.bprmMn".format(etomoRootName), str(2*binDimX), str(2*binDimY), str(dimX), str(dimY), str(dimZ-1))], shell=True)

# Reconstruct
BP = etomoRootName + '_ewbp.xzyt'
call(["ewbp {} {} -reconxz={}:{} -sizexz={}:{} -iy=0:{} -filter=2 -hdfilt=1 -moderec=2".format(prmStack, BP, str(binDimX), str(binDimZ), str(binDimX), str(binDimZ), str(binDimY-1))], shell=True)
