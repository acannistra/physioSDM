## get_prism_climate function for extraction of prism climate data

def get_prism_climate(points, year, month, var, prismloc="~/prismtmp/"):
    import rasterio
    import glob
    from os.path import expanduser
    prismloc = expanduser(prismloc)
    search_str = prismloc + \
        "*/PRISM_%s_stable_*_%d%02d_bil.bil" % (var, year, month)
    climfile = glob.glob(search_str)

    if len(climfile) is not 1:
        import errno
        import os
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), search_str)

    climate_data = []
    with rasterio.open(climfile[0]) as climate:
        climate_data = climate.sample(points)
        return list(climate_data)