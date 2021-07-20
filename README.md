# mcmicro-flowSOM
An MCMICRO module for clustering cell types using the flowSOM algorithm.

## Example Usage
```
docker run --rm -v "$PWD":/data labsyspharm/mc-flowsom:1.0.0 python3 /app/cluster.py -i /data/unmicst-2.csv -o /data/ -c
```

## Parameter Reference
```
usage: cluster.py [-h] -i INPUT [-o OUTPUT] [-m MARKERS] [-v] [-c] [-n NUM_METACLUSTERS]

Cluster cell types using mcmicro marker expression data.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV of mcmicro marker expression data for cells
  -o OUTPUT, --output OUTPUT
                        The directory to which output files will be saved
  -m MARKERS, --markers MARKERS
                        A text file with a marker on each line to specify which markers to use for clustering
  -v, --verbose         Flag to print out progress of script
  -c, --method          Include a column with the method name in the output files.
  -n NUM_METACLUSTERS, --num-metaclusters NUM_METACLUSTERS
                        number of clusters for meta-clustering. Default is 25.
```