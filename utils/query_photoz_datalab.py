# Author: Igor Andreoni
# taken from https://github.com/igorandreoni/snippets
# modified by A Moller Feb 2022
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

from dl import queryClient as qc


def query_coords_ls(
    idx,
    ra,
    dec,
    radius_arcsec=5,
    radius_nuclear=1.0,
    catalog="dr8_north",
    datalab=True,
    check_quality=True,
):
    """Query the database to search for matches at the given coordinates"""

    # convert to int
    ra = float(ra)
    dec = float(dec)

    # Crossmatch with photoz database
    if datalab is True:
        radius_deg = radius_arcsec / 3600.0
        query = qc.query(
            sql=f"SELECT z_phot_median, z_phot_std, z_phot_l95, ra, dec, \
                             type, flux_z from ls_dr8.photo_z INNER JOIN ls_dr8.tractor \
                             ON ls_dr8.tractor.ls_id = ls_dr8.photo_z.ls_id \
                             where ra > ({ra-radius_deg}) and \
                             ra < ({ra+radius_deg}) and \
                             dec > ({dec-radius_deg}) and \
                             dec < ({dec+radius_deg})"
        )
        result0 = query.split("\n")
        result0 = [r.split(",") for r in result0][1:-1]

        ras = [float(r[3]) for r in result0]
        decs = [float(r[4]) for r in result0]
        result = []
        if len(ras) > 0:
            # Add separation
            gal_coords = SkyCoord(ra=ras * u.deg, dec=decs * u.deg)
            cand_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
            sep = cand_coord.separation(gal_coords)
            for i in np.arange(len(ras)):
                result0[i].append(float(sep[i].arcsec))
                # Check that the separation is less than required
                if float(sep[i].arcsec) < radius_arcsec:
                    result.append(result0[i])

        df = pd.DataFrame(
            columns=[
                "id",
                "ls_z_phot_median_1",
                "ls_z_phot_std_1",
                "ls_z_phot_l95_1",
                "ls_ra_1",
                "ls_dec_1",
                "ls_type_1",
                "ls_flux_z_1",
                "ls_sep_arcsec_1",
                "ls_z_phot_median_2",
                "ls_z_phot_std_2",
                "ls_z_phot_l95_2",
                "ls_ra_2",
                "ls_dec_2",
                "ls_type_2",
                "ls_flux_z_2",
                "ls_sep_arcsec_2",
            ],
            index=[0],
        )
        if len(result) > 1:
            df.loc[0] = [idx] + result[0] + result[1]
        if len(result) == 0:
            df.loc[0] = [idx] + [None for i in range(16)]
        else:
            df.loc[0] = [idx] + result[0] + [None for i in range(8)]
        return df
        # # Print the results
        # print(
        #     "z_phot_median, z_phot_std, z_phot_l95, ra, dec, type, flux_z, sep_arcsec"
        # )
        # for r in result:
        #     print(", ".join([str(x) for x in r]))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Legacy Survey phtotoz")
    parser.add_argument(
        "radec", metavar="RA, Dec", type=str, nargs="+", help="RA and Dec (degrees)"
    )
    parser.add_argument(
        "-r",
        dest="radius",
        type=float,
        required=False,
        help="Search radius (arcsc)",
        default=2,
    )
    args = parser.parse_args()

    # RA and Dec
    ra, dec = float(args.radec[0]), float(args.radec[1])

    # Radius
    search_rad = args.radius
    df = query_coords_ls(ra, dec, radius_arcsec=search_rad)
