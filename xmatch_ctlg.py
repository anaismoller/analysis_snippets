""" XMatching a catalogue
    Ingest candidates coordinates & cross-match to catalogues
    Code inspired from Fink broker utilities
    questions amoller@swin.edu.au
"""
import os
import argparse
import pandas as pd
from pathlib import Path
from utils import xmatch
from utils import query_photoz_datalab as photoz


def convert_to_degrees(row, ra_key="RA", dec_key="dec", key_toprocess=""):
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(row[ra_key], row[dec_key], unit=(u.hourangle, u.deg))
    if key_toprocess == ra_key:
        return c.ra.degree
    elif key_toprocess == dec_key:
        return c.dec.degree
    else:
        print("Wrong key")
        print(ra_key, dec_key, key_toprocess)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process candidates metadata")

    parser.add_argument(
        "--path_in",
        default="./test_files/transients_coo.txt",
        help="Candidate metadata (coords) from Mary",
    )
    parser.add_argument(
        "--path_out", default="./dump", help="Path to output enriched metadata",
    )
    parser.add_argument(
        "--outname",
        default="./transients_coo_enriched",
        help="Output name enriched metadata (no extension)",
    )

    parser.add_argument(
        "--test", action="store_true", help="one file processed only",
    )
    parser.add_argument(
        "--delimiter", default=" ", help="delimiter to be used",
    )
    parser.add_argument(
        "--ra_key", default="ra", help="ra key",
    )
    parser.add_argument(
        "--dec_key", default="dec", help="dec key",
    )
    parser.add_argument(
        "--id_key", default="id", help="id key",
    )
    args = parser.parse_args()

    ra_key = args.ra_key
    dec_key = args.dec_key
    id_key = args.id_key

    # df = du.read_mary_masterlist(args.path_in)
    df = pd.read_csv(args.path_in, delimiter=args.delimiter)

    df["ra"] = df.apply(
        lambda row: row[ra_key]
        if ":" not in str(row[ra_key])
        else convert_to_degrees(
            row, key_toprocess=ra_key, ra_key=ra_key, dec_key=dec_key
        ),
        axis=1,
    )

    df["dec"] = df.apply(
        lambda row: row[dec_key]
        if ":" not in str(row[dec_key])
        else convert_to_degrees(
            row, key_toprocess=dec_key, ra_key=ra_key, dec_key=dec_key
        ),
        axis=1,
    )

    df = df.rename(columns={id_key: "indx"})

    if args.test:
        df = df[:10]
    print(f"Processing {args.path_in}: {len(df)} candidates")

    print("SIMBAD xmatch")
    z, sptype, typ, ctlg = xmatch.cross_match_simbad(
        df["indx"].to_list(), df["ra"].to_list(), df["dec"].to_list()
    )
    # save in df
    df["simbad_type"] = typ
    df["simbad_ctlg"] = ctlg
    df["simbad_sptype"] = sptype
    df["simbad_redshift"] = z

    print("GAIA2 + e3 xmatch")
    source, ragaia, decgaia, plx, plxerr, gmag, angdist = xmatch.cross_match_gaia(
        df["indx"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/345/gaia2",
    )
    (
        source_edr3,
        ragaia_edr3,
        decgaia_edr3,
        plx_edr3,
        plxerr_edr3,
        gmag_edr3,
        angdist_edr3,
    ) = xmatch.cross_match_gaia(
        df["indx"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/350/gaiaedr3",
    )
    # save in df
    df["gaia_DR2_source"] = source
    df["gaia_DR2_ra"] = ragaia
    df["gaia_DR2_dec"] = decgaia
    df["gaia_DR2_parallax"] = plx
    df["gaia_DR2_parallaxerr"] = plxerr
    df["gaia_DR2_gmag"] = gmag
    df["gaia_DR2_angdist"] = angdist
    df["gaia_eDR3_source"] = source_edr3
    df["gaia_eDR3_ra"] = ragaia_edr3
    df["gaia_eDR3_dec"] = decgaia_edr3
    df["gaia_eDR3_parallax"] = plx_edr3
    df["gaia_eDR3_parallaxerr"] = plxerr_edr3
    df["gaia_eDR3_gmag"] = gmag_edr3
    df["gaia_eDR3_angdist"] = angdist_edr3

    print("USNO-A.20 xmatch")
    (source_usno, angdist_usno,) = xmatch.cross_match_usno(
        df["indx"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:I/252/out",
    )
    df["USNO_source"] = source_usno
    df["USNO_angdist"] = angdist_usno

    print("Legacy Survey xmatch")
    list_ls_df = []
    for (idx, ra, dec) in df[["indx", "ra", "dec"]].values:
        list_ls_df.append(photoz.query_coords_ls(idx, ra, dec, radius_arcsec=10))
    df_ls = pd.concat(list_ls_df)
    df_ls = df_ls.rename(columns={"id": "indx"})
    print("Finished Legacy Survey xmatch")
    df = pd.merge(df, df_ls, on="indx")

    print("WISE xmatch")
    df_wise = xmatch.cross_match_alerts_raw_generic(
        df["indx"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:II/311/wise",
        distmaxarcsec=5,
    )
    df_wise = df_wise.rename(
        columns={
            col: col + "_wise"
            for col in df_wise.columns
            if col not in ["indx", "ra", "dec"]
        }
    )
    df_wise["indx"] = df_wise["indx"].astype(int)
    df_wise["ra"] = df_wise["ra"].astype(float)
    df_wise["dec"] = df_wise["dec"].astype(float)

    df["ra"] = df["ra"].astype(float)
    df["dec"] = df["dec"].astype(float)
    df = pd.merge(df, df_wise, on=["indx", "ra", "dec"], how="left")

    print("SkyMapper xmatch")
    df_smss = xmatch.cross_match_alerts_raw_generic(
        df["indx"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:II/358/smss",
        distmaxarcsec=5,
    )
    df_smss = df_smss.rename(
        columns={
            col: col + "_SMSS"
            for col in df_smss.columns
            if col not in ["indx", "ra", "dec"]
        }
    )
    df_smss["indx"] = df_smss["indx"]
    df_smss["ra"] = df_smss["ra"].astype(float)
    df_smss["dec"] = df_smss["dec"].astype(float)
    df = pd.merge(df, df_smss, on=["indx", "ra", "dec"], how="left")

    print("TESS xmatch")
    df_tess = xmatch.cross_match_alerts_raw_generic(
        df["indx"].to_list(),
        df["ra"].to_list(),
        df["dec"].to_list(),
        ctlg="vizier:IV/39/tic82",
        distmaxarcsec=5,
    )
    df_tess = df_tess.rename(
        columns={
            col: col + "_TESS"
            for col in df_tess.columns
            if col not in ["indx", "ra", "dec"]
        }
    )
    df_tess["indx"] = df_tess["indx"]
    df_tess["ra"] = df_tess["ra"].astype(float)
    df_tess["dec"] = df_tess["dec"].astype(float)
    df = pd.merge(df, df_tess, on=["indx", "ra", "dec"], how="left")

    # output
    os.makedirs(args.path_out, exist_ok=True)
    outname = f"{args.path_out}/{Path(args.path_in).stem}_xmatched.csv"
    df.to_csv(outname)
    print(f"Xmatched catalgoue saved {outname}")

