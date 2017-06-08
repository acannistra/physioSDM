import pandas as pd
import numpy
import geopandas as gpd
import argparse
import os
from itertools import chain
from prism_utils import get_prism_climate

parser = argparse.ArgumentParser(
    description="Add PRISM climate variables to geospatial point data.")

parser.add_argument("pointfile", help="CSV Datafile containting point data.")
parser.add_argument("prismvar", help="PRISM Climate Variable to append.")
parser.add_argument("--lat", help="Latitude Column Name",
                    nargs=1,
                    required=True)
parser.add_argument("--lon", help="Longitude Column Name",
                    nargs=1,
                    required=True)
parser.add_argument("--sep",
                    help='Datafile column separation character.',
                    default='\t')
parser.add_argument("--prism_loc", help="PRISM climate files location.",
                    default='~/prismtmp')
parser.add_argument("--timescale", help="Resolution of climate data to use.",
                    choices=['day', 'month', 'year'], default='month')
parser.add_argument("--pandas_query", help="provide a pandas query expression"
                    " to apply to the data prior to adding climate.", 
                    nargs=1)
parser.add_argument("--outsuffix", help="add suffix to outfile name.")
args = parser.parse_args()

#---------------------


pointdata = pd.read_csv(args.pointfile, sep=args.sep)
print("Points: %d"%len(pointdata))
print("Columns: %s"%pointdata.columns)

if(args.pandas_query):
    args.pandas_query = args.pandas_query[0]
    nrows_before = len(pointdata)
    pointdata = pointdata.query(args.pandas_query)
    nrows_after = nrows_before-len(pointdata)
    print("Removed %d rows with query %s" %(nrows_after, args.pandas_query))


points_with_climate = pd.DataFrame()

if(args.timescale == 'day'):
    raise NotImplementedError("Day is not implemented yet. ")
elif(args.timescale == 'month'):
    groups = pointdata.groupby(['year', 'month'])
    for name, group in groups:
        year, month = int(name[0]), int(name[1])
        group = pd.DataFrame(group)
        lons = list(chain(*group[args.lon].values))
        lats = list(chain(*group[args.lat].values))
        points = zip(lons, lats)
        try:
            climate = get_prism_climate(points, year, month, args.prismvar,
                                        args.prism_loc)
            group[args.prismvar] = list(chain(*climate))
        except FileNotFoundError as e:
            print(str(e) + " appending NaNs.")
            group[args.prismvar] = np.nan

        points_with_climate  = points_with_climate.append(group)


elif(args.timescale == 'year'):
    raise NotImplementedError("Year is not implemented yet.")
else:
    raise Exception("Should ne'er be here!")

## -----
if (not args.outsuffix):
    outfile = "%s.%s.csv"%(os.path.splitext(args.pointfile)[0],
                           args.prismvar)
else:
    outfile = "%s.%s_%s.csv"%(os.path.splitext(args.pointfile)[0],
                              args.prismvar, args.outsuffix)




points_with_climate.to_csv(outfile, sep=args.sep)
