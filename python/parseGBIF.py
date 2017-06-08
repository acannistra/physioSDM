import pandas as pd
import geopandas as gpd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("occurrence_file")
parser.add_argument("--date_range",
                    action='append',
                    nargs=1,
                    help="Supply date range in YYYY:YYYY format."
                         " Multiple ranges can be supplied, each producing"
                         " A different output file.",
                    dest="date_ranges")

print(parser.parse_args())
