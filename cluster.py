import re
import os
import argparse
import subprocess
import pandas as pd

'''
Parse arguments.
Input file is required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Cluster cell types using mcmicro marker expression data.')
    parser.add_argument('-i', '--input', help="Input CSV of mcmicro marker expression data for cells", type=str, required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved', type=str, required=False)
    parser.add_argument('-m', '--markers', help='A text file with a marker on each line to specify which markers to use for clustering', type=str, required=False)
    parser.add_argument('-v', '--verbose', help='Flag to print out progress of script', action="store_true", required=False)
    parser.add_argument('-c', '--method', help='Include a column with the method name in the output files.', action="store_true", required=False)
    parser.add_argument('-n', '--num-metaclusters', help='number of clusters for meta-clustering. Default is 25.', type=int, required=False, default=25)
    args = parser.parse_args()
    return args


'''
Get the path to the directory where this script is located and return it.
'''
def get_path():
    full_path = os.path.realpath(__file__)
    path_list = full_path.split('/')[:-1]
    s = '/'
    return s.join(path_list)


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1] # get filename from end of input path
    dataName = fileName[:fileName.rfind('.')] # get data name by removing extension from file name
    return dataName


'''
Get markers to use for clustering from a text file where each marker is on a line and corresponds exactly to the column name in the input data file.
Returns a list of markers to use for clustering.
'''
def get_markers(markers_file):
    markers = [] # list of markers in file

    # read markers from file
    f = open(markers_file, 'r')
    for line in f:
        markers.append(line.strip())

    return markers


'''
Clean data in input file.
NOTE: Currently we are doing this with pandas however, using csv might be faster, or numpy.
Exclude the following data from clustering:
    - X_centroid, …, Extent, Orientation - morphological features
    - Any of the following DNA stains
        - DNA0, DNA1, DNA2, …
        - Hoechst0, Hoechst1, ....
        - DAPI0, DAPI1, …
    - AF* (e.g., AF488, AF555, etc.) - autofluorescence
    - A* (e.g., A488, A555, etc.) - secondary antibody staining only
To include any of these markers in the clustering, provide their exact names in a file passed in with the '-m' flag
'''
def clean(input_file):

    # constants

    # a default list of features to exclude from clustering
    FEATURES_TO_REMOVE = ['X_centroid', 'Y_centroid', # morphological features
                        'column_centroid', 'row_centroid', 
                        'Area', 'MajorAxisLength', 
                        'MinorAxisLength', 'Eccentricity', 
                        'Solidity', 'Extent', 'Orientation', 
                        'DNA.*', 'Hoechst.*', 'DAP.*', # DNA stain
                        'AF.*', # autofluorescence
                        'A\d{3}.*'] # secondary antibody staining only (iy has to have 3 digist after)

    if args.verbose:
        print('Cleaning data...')

    # load csv
    data = pd.read_csv(input_file)

    # if markers provided, keep only those features and the Cell IDs. It is important that the CellID column is first.
    if args.markers:
        if CELL_ID not in markers: # add cell ID to list of columns to keep
            markers.insert(0, CELL_ID)
        elif markers.index(CELL_ID) != 0: # if cell ID column is included but not first, move it to the front
            markers.insert(0, markers.pop(markers.index(CELL_ID)))
        data = data[markers]
    else:
        # find any columns in the input csv that should be excluded from clustering be default
        # NOTE: may want to replace this with regex, it might be faster.
        col_to_remove = []
        cols = data.columns
        for feature in FEATURES_TO_REMOVE:
            r = re.compile(feature)
            col_to_remove.extend(list(filter(r.match, cols)))
        
        # drop all columns that should be excluded
        data = data.drop(columns=col_to_remove, axis=1)

    # save cleaned data to csv
    data.to_csv(f'{output}/{clean_data_file}', index=False)

    if args.verbose:
        print(f'Done. Cleaned data is in {output}/{clean_data_file}.')


'''
Run an R script that converts the clean data csv into fcs format
'''
def convertToFSC():
    if args.verbose:
        print('Converting CSV to FCS...')

    path = get_path() # get the path where the r script is located

    # Build subprocess command
    r_script = ['Rscript', f'{path}/CSVtoFCS.r'] # use CSVtoFCS.r script

    # pass args
    r_args = [output, clean_data_file]

    # Build subprocess command
    command = r_script + r_args

    # run it
    subprocess.run(command, universal_newlines=True)

    if args.verbose:
        print('Done.')


'''
Run an R script that runs flowSOM
'''
def runFlowSOM():
    if args.verbose:
        print('Running flowSOM...')

    path = get_path() # get the path where the r script is located

    r_script = ['Rscript', f'{path}/runFlowSOM.r'] # use FastPG.r script
    # pass input data file, k value, number of cpus to use for the k nearest neighbors part of clustering, output dir, cells file name, clusters file name
    r_args = [f'{output}/{clean_data_fcs_file}', str(args.num_metaclusters), str(args.method), output, cells_file, clusters_file]

    # Build subprocess command
    command = r_script + r_args

    # run it
    subprocess.run(command, universal_newlines=True)

    if args.verbose:
        print('Done.')


'''
Main.
'''
if __name__ == '__main__':
    args = parseArgs() # parse arguments

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    # get list of markers if provided
    if args.markers is not None:
        markers = get_markers(args.markers)

    # constants
    CELL_ID = 'CellID' # column name holding cell IDs
    
    # output file names
    data_prefix = getDataName(args.input) # get the name of the input data file to add as a prefix to the output file names
    clean_data_file = f'{data_prefix}-clean.csv' # name of output cleaned data CSV file
    clean_data_fcs_file = f'{data_prefix}-clean.fcs' # name of output cleaned data CSV file
    clusters_file = f'{data_prefix}-clusters.csv' # name of output CSV file that contains the mean expression of each feaute, for each cluster
    cells_file = f'{data_prefix}-cells.csv' # name of output CSV file that contains each cell ID and it's cluster assignation
    
    # clean input data file
    clean(args.input)

    # convert the clean input data to FCS format
    convertToFSC()

    # run flowSOM algorithm
    runFlowSOM()