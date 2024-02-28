#!/usr/bin/python

# History: Jun 12 2019 Initial coding. This was modified from extract.py.

# To do: Make rows= option work for VCF files.

from __future__ import division;
from __future__ import print_function;
import sys;
import os;
from operator import *;
import gzip;
import math;
import time;
import ntpath;

# Function to return a file id for reading
def getFIDr(file, gz=None):
    
  if not os.path.isfile(file):
    print(file, ' does not exist.')
    sys.exit(0) 
  
  if gz is None:
    slen = len(file)
    gz   = 0
    if file[(slen-3):] == '.gz': gz = 1;
  if (gz):
    fid = gzip.open(file, 'r')
  else:
    fid = open(file, 'r')

  return(fid)

# Function to return the file format ("impute" or "vcf")
def getFileFormat(file, gz=None):

  fid  = getFIDr(file, gz=gz);
  temp = fid.readline(1);
  temp = temp.decode('utf-8');
  fid.close();
  if ((temp == "#") or (temp == b'#')):
    form = "vcf"
  else:
    form = "impute"
  
  return(form)

# Function to return a file id for writing
def getFIDo(file, gz=None):
    
  if gz is None:
    slen = len(file)
    gz   = 0
    if file[(slen-3):] == '.gz': gz = 1;
  if (gz):
    base = ntpath.basename(file)
    fobj = open(file, 'wb')
    fid  = gzip.GzipFile(base, fileobj=fobj)
  else:
    fid = open(file, 'w')

  return(fid)

# Function to close tmp files
def closeFiles(FIDs):

  for i in range(len(FIDs)): FIDs[i].close(); 

  return(None);

# Return the header row for a VCF file and position the file pointer
def getHeader_vcf(file):
 
  fid  = getFIDr(file);
  flag = 0;

  # Start reading and search for #CHROM
  for i, line in enumerate(fid):  
    line = line.decode('UTF-8');
    if (len(line) > 10):
      temp = line[0:10].split()[0];
      if (temp.upper() == "#CHROM"):
        flag = 1;
        break
  closeFiles([fid]);

  if not flag:
    print('ERROR with VCF file: no header row found');
    sys.exit(0);

  return(line)

def checkForSep(ff):
  
  n = len(ff);
  if n == 0:
    ff = './';
  else:
    if ff[n-1] != '/': ff = ff+'/';

  return(ff);

def parseOptions(opList, required=None):
  
  n = len(opList);
  if n == 0: return {};

  ret = {};
  for i in range(n):
    vec = opList[i].split('=');
    if len(vec) != 2:
      print('ERROR in parseOptions: The option ', opList[i], ' is not correct');
      sys.exit(0);
    ret[vec[0]] = vec[1];

  # Check the required options
  if required != None:
    n    = len(required);
    keys = ret.keys();
    for i in range(n):
      if required[i] not in keys:
        print('ERROR in parseOptions: ', required[i], '= is required');
        sys.exit(0);

  return(ret);

# Function to get the row numbers. Each line gives the row numbers to read in for each genotype file
def getFileAndLocs(fline):
  
  fline = fline.strip("\n");
  vec   = fline.split(",");
  n     = len(vec);
  if (n < 2):
     print('ERROR: no file or locations specified');
     sys.exit(0);
 
  gfile = vec[0];
  locs  = [int(vec[i]) for i in range(1, n)];
 
  return(gfile, locs);

def getCols(colFile):

  cols  = getFIDr(colFile).readlines();
  ncols = len(cols);
  for i in range(ncols): cols[i] = int(cols[i].strip("\n")) - 1;

  return(cols);

def getLocs(locList):

  locs  = getFIDr(locList).readlines();
  nlocs = len(locs);
  for i in range(nlocs): locs[i] = int(locs[i].strip("\n"));
 
  return(locs);


def getLineGivenCols(line, cols):

  line = line.strip("\n");
  vec  = line.split();
  vec2 = [vec[k] for k in cols]
  line = " ".join(vec2);
  line = line + "\n";

  return(line);


def readFileNames(colfile):

  fid    = getFIDr(colfile);
  files  = fid.readlines();
  fid.close();
  nfiles = len(files);
  if nfiles:
    for i in range(nfiles): files[i] = files[i].strip();

  return(files);

# Function to set the options
def setOptions(opDict, snpFile):

    if opDict is None: opDict = {};
    colFile = opDict.get('cols', '');
    colFlag = len(colFile) > 0;
    print01 = int(opDict.get('print', 1));
    colVec  = None;
    locCol  = opDict.get('locCol', 3);
    locCol  = int(locCol) - 1;
    range   = int(opDict.get('range', 0));

    # For cols. Subtract 1
    if (colFlag): colVec = getCols(colFile)

    opDict.update({'colFlag':colFlag, 'locCol':locCol, 'print':print01,
                   'colVec':colVec, 'range':range});

    return(opDict);

# Function to get the snp (location) from a line
def getLoc(line, locCol, maxChar, sep=None):

  ret = -1;
  vec = line[0:maxChar].split(sep);
  if (locCol >= len(vec)-1): vec = line.split(sep);
  if (locCol <= len(vec)): ret = int(vec[locCol]);
  
  return(ret);

# Function to return a flag to use the line and return the snp id
def getUse(line, locs, locCol, minLoc, maxLoc, maxChar=100, sep=None, isRange=0):

  # currentRow is the current row, this index starts from 1
  use       = False;
  loc       = -1;
  breakFlag = False;
  linelen   = len(line);
  if (line[0] == "#"): return(use, loc, breakFlag);
  if (not linelen): return(use, loc, breakFlag);
  if (linelen < maxChar): maxChar = linelen;    

  loc = getLoc(line, locCol, maxChar, sep);
  if (isRange):
    if ((loc >= minLoc) and ((loc <= maxLoc) or (maxLoc == -1))):    
      use = True;
    elif ((loc > maxLoc) and (maxLoc > 0)):
      breakFlag = True;
  else:
    if (loc > maxLoc):
      breakFlag = True;
    elif (loc > -1):
      if (loc in locs): use = True;
    else:
      use = False;
    
  return(use, loc, breakFlag);

# Function to output the line ans snp to files
def outputLine(line, outfid, vcfFlag, headerInOutFile, header):

   ret = headerInOutFile;
  
   if ((vcfFlag) and (not headerInOutFile)):
      outfid.write(header);
      ret = 1;

   # Write out the line
   outfid.write(line);
   outfid.flush();

   return(ret);

# Function to extract locs from a single file
def ExtractFromFile(infile, locCol, locs, outfid, colVec, sep=None, isRange=0):

   if ((isRange) and (len(locs) != 2)): isRange = 0; 
   if (isRange):
     minLoc = locs[0];
     maxLoc = locs[1];
   else:
     minLoc = min(locs);
     maxLoc = max(locs);
   if (colVec != None):
      colFlag = 1;
   else:
      colFlag = 0; 

   # Open the genotype file
   fid  = getFIDr(infile);
   loc0 = -1;

   # Read in row by row as a character string
   for i, line in enumerate(fid):  

      # To prevent errors on Windows with binary characters
      line = line.decode('UTF-8');

      # Determine if we want this line
      use, loc, breakFlag = getUse(line, locs, locCol, minLoc, maxLoc, maxChar=100, 
                                   sep=sep, isRange=isRange);

      # If the current row is greater than the maximum row number to read in then stop reading
      if breakFlag: break;

      # Check that locations are increasing
      if ((loc != -1) and (loc0 > loc)):
         print("ERROR: locations are not increasing");
         print([loc0, loc]);
         sys.exit(0);
      else: loc0 = loc; 

      if use:
         # Re-construct the line to output with the columns that are needed
         if colFlag: line = getLineGivenCols(line, colVec);
         outfid.write(line);
         outfid.flush();

   fid.close();

def openAndInitOutFile(outfile, gfile):
   
   form = getFileFormat(gfile);
   if (form == "vcf"):
      vcfFlag = 1;
   else:
      vcfFlag = 0;
   outfid = getFIDo(outfile);
   if (vcfFlag):
      header = getHeader_vcf(gfile);   
      outfid.write(header);
      outfid.flush();

   return(outfid, vcfFlag);

# Main function
def ExtractByLoc(snpFile, outFile, opDict=None):

    # Check the options
    opDict  = setOptions(opDict, snpFile);

    print01 = opDict.get('print');
    colVec  = opDict.get('colVec', None);
    locCol  = opDict.get('locCol');
    sep     = opDict.get('sep', None);
    isRange = opDict.get('range', 0);

    del opDict;
    headerInOutFile = 0;
    
    # Open file of geno files and locations
    snpfid = getFIDr(snpFile);
    for j, fline in enumerate(snpfid):
       # To prevent errors on Windows with binary characters
       #fline = fline.decode('UTF-8');
      
       # Get the file to read and desired locations
       gfile, locs = getFileAndLocs(fline);

       # Open and initialize output file
       if (headerInOutFile == 0):
          outfid, vcfFlag = openAndInitOutFile(outFile, gfile);
          headerInOutFile = 1;

       # Extract locs from gfile
       if print01: print("Reading file: ", gfile);
       ExtractFromFile(gfile, locCol, locs, outfid, colVec, sep=sep, isRange=isRange);
    snpfid.close();
    outfid.close();

    return(None);
      

    


if __name__=="__main__":

    M = len(sys.argv);
    if M < 3:
        print("\nUsage: python extractByLoc.py files=<genotype files> out=<file> <options> \n");
        print("Required:");
        
        print("files=<file>  Comma seperated file with file names and locations to extract.");
        print("              The format is:");
        print("              file_1,loc11,loc12,...,loc1n:");
        print("              file_2,loc21,loc22,...,loc2n:");
        print("              Each file_i must have the same format (same delimiter and location column).");

        print("out=<file>    Output file. File will be compressed it has a .gz extension.\n");

        print("Optional:");
        print("locCol=<int>  Column number (starting from 1) for the location column in each file. Default is 3.");
        print("cols=<file>   File with a single column of column numbers (starting from 1) to extract.");
        print("print=<0,1>   To print information to the console. Default is 1.");
        print("sep=<character> File delimiter. Default is white space.");
        print("range=<0,1> For only two locations that will define the range loc1-loc2.");
        
        print("\n");
        sys.exit(0);

   
    myDict = parseOptions(sys.argv[1:M], required=['files', 'out']);
    ExtractByLoc(myDict['files'], myDict['out'], opDict=myDict);
    