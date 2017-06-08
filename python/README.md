# Python Utilities

## `append_prism.py`

```
usage: append_prism.py [-h] --lat LAT --lon LON [--sep SEP]
                       [--prism_loc PRISM_LOC] [--timescale {day,month,year}]
                       [--pandas_query PANDAS_QUERY]
                       pointfile prismvar

Add PRISM climate variables to geospatial point data.

positional arguments:
  pointfile             CSV Datafile containting point data.
  prismvar              PRISM Climate Variable to append.

optional arguments:
  -h, --help            show this help message and exit
  --lat LAT             Latitude Column Name
  --lon LON             Longitude Column Name
  --sep SEP             Datafile column separation character.
  --prism_loc PRISM_LOC
                        PRISM climate files location.
  --timescale {day,month,year}
                        Resolution of climate data to use.
  --pandas_query PANDAS_QUERY
                        provide a pandas query expression to apply to the data
                        prior to adding climate.
```

## `prism_utils.py`
Contains functions used by `append_prism.py` to extract PRISM data. 

