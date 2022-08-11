# Year 2022
# Based on fink-broker.org code
# https://github.com/astrolabsoftware/fink-science/tree/master/fink_science/xmatch
# Adapted by Anais Möller

import io
import csv
import logging
import requests
import numpy as np
import pandas as pd


def generate_csv(s: str, lists: list) -> str:
    """ Make a string (CSV formatted) given lists of data and header.
    Parameters
    ----------
    s: str
        String which will contain the data.
        Should initially contain the CSV header.
    lists: list of lists
        List containing data.
        Length of `lists` must correspond to the header.
    Returns
    ----------
    s: str
        Updated string with one row per line.
    Examples
    ----------
    >>> header = "toto,tata\\n"
    >>> lists = [[1, 2], ["cat", "dog"]]
    >>> table = generate_csv(header, lists)
    >>> print(table)
    toto,tata
    1,"cat"
    2,"dog"
    <BLANKLINE>
    """
    output = io.StringIO()
    writer = csv.writer(output, quoting=csv.QUOTE_NONNUMERIC)
    _ = [writer.writerow(row) for row in zip(*lists)]
    return s + output.getvalue().replace("\r", "")


def refine_search(
    ra: list,
    dec: list,
    oid: list,
    id_out: list,
    names: list,
    types: list,
    sptypes: list,
    redshift: list,
) -> list:
    """ Create a final table by merging coordinates of objects found on the
    bibliographical database, with those objects which were not found.
    Parameters
    ----------
    ra: list of float
        List of RA
    dec: list of float
        List of Dec of the same size as ra.
    oid: list of str
        List of object ID (custom)
    id_out: list of str
        List of object ID returned by the xmatch with CDS
    names: list of str
        For matches, names of the celestial objects found
    types: list of str
        For matches, astronomical types of the celestial objects found
    sptypes: list of str
        For matches, spectral types of the celestial objects found
    redshift: list of str
        For matches, astronomical redshifts of the celestial objects found
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (index, ra, dec, name, type,sptype,redshift).
        If the object is not found in Simbad, name & type
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    """
    out = []
    for ra_in, dec_in, id_in in zip(ra, dec, oid):
        # cast for picky Spark
        ra_in, dec_in = float(ra_in), float(dec_in)
        id_in = str(id_in)

        # Discriminate with the index
        if id_in in id_out:
            # Return the closest object in case of many
            # (smallest angular distance)
            index = id_out.index(id_in)
            sp_type_tmp = sptypes[index] if sptypes[index] != "" else "Unknown"
            redshift_tmp = redshift[index] if redshift[index] != "" else "Unknown"
            out.append(
                (
                    id_in,
                    ra_in,
                    dec_in,
                    str(names[index]),
                    str(types[index]),
                    str(sp_type_tmp),
                    str(redshift_tmp),
                )
            )

        else:
            # Mark as unknown if no match
            out.append(
                (id_in, ra_in, dec_in, "Unknown", "Unknown", "Unknown", "Unknown")
            )

    return out


def xmatch(
    ra: list, dec: list, id: list, extcatalog: str = "simbad", distmaxarcsec: int = 1
) -> (list, list):
    """ 
    
    Build a catalog of (ra, dec, id) in a CSV-like string,
    cross-match with `extcatalog`, and decode the output.
    See http://cdsxmatch.u-strasbg.fr/ for more information.
    Parameters
    ----------
    ra: list of float
        List of RA
    dec: list of float
        List of Dec of the same size as ra.
    id: list of str
        List of object ID (custom)
    extcatalog: str
        Name of the catalog to use for the xMatch.
        See http://cdsxmatch.u-strasbg.fr/ for more information.
    distmaxarcsec: int
        Radius used for searching match. extcatalog sources lying within
        radius of the center (ra, dec) will be considered as matches.
    Returns
    ----------
    data: list of string
        Unformatted decoded data returned by the xMatch
    header: list of string
        Unformatted decoded header returned by the xmatch
    """
    # Build a catalog of alert in a CSV-like string
    table_header = """ra_in,dec_in,index\n"""
    table = generate_csv(table_header, [ra, dec, id])

    # Send the request!
    r = requests.post(
        "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync",
        data={
            "request": "xmatch",
            "distMaxArcsec": distmaxarcsec,
            "selection": "all",
            "RESPONSEFORMAT": "csv",
            "cat2": extcatalog,
            "colRA1": "ra_in",
            "colDec1": "dec_in",
        },
        files={"cat1": table},
    )

    # Decode the message, and split line by line
    # First line is header - last is empty
    data = r.content.decode().split("\n")[1:-1]
    header = r.content.decode().split("\n")[0].split(",")

    return data, header


def cross_match_alerts_raw_simbad(oid: list, ra: list, dec: list) -> list:
    """ Query the CDSXmatch service to find identified objects
    in alerts. The catalog queried is the SIMBAD bibliographical database.
    We can also use the 10,000+ VizieR tables if needed :-)
    Parameters
    ----------
    oid: list of str
        List containing object ids (custom)
    ra: list of float
        List containing object ra coordinates
    dec: list of float
        List containing object dec coordinates
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (index, ra, dec, name, type).
        If the object is not found in Simbad, name & type
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    Examples
    ----------
    >>> ra = [26.8566983, 26.24497]
    >>> dec = [-26.9677112, -26.7569436]
    >>> id = ["1", "2"]
    >>> objects = cross_match_alerts_raw(id, ra, dec)
    >>> print(objects) # doctest: +NORMALIZE_WHITESPACE
    [('1', 26.8566983, -26.9677112, 'TYC 6431-115-1', 'Star'),
     ('2', 26.24497, -26.7569436, 'Unknown', 'Unknown')]
    """
    if len(ra) == 0:
        return []

    # Catch TimeoutError and ConnectionError
    try:
        data, header = xmatch(ra, dec, oid, extcatalog="simbad", distmaxarcsec=5)
    except (ConnectionError, TimeoutError, ValueError) as ce:
        logging.warning("XMATCH failed " + repr(ce))
        return []

    # Sometimes the service is down, but without TimeoutError or ConnectionError
    # In that case, we grab the error message from the data.
    if len(data) > 0 and "504 Gateway Time-out" in data[0]:
        msg_head = "CDS xmatch service probably down"
        msg_foot = "Check at http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
        logging.warning(msg_head)
        logging.warning(data[0])
        logging.warning(msg_foot)
        return []

    # Fields of interest (their indices in the output)
    if "main_id" not in header:
        return []

    # Fields of interest (their indices in the output)
    main_id = header.index("main_id")
    main_type = header.index("main_type")
    oid_ind = header.index("index")
    redshift_ind = header.index("redshift")
    sp_type_ind = header.index("sp_type")

    # Get the index of matches
    id_out = [np.array(i.split(","))[oid_ind] for i in data]

    # Get the names of matches
    names = [np.array(i.split(","))[main_id] for i in data]

    # Get the types of matches
    types = [np.array(i.split(","))[main_type] for i in data]

    # Get the types of matches
    sp_types = [np.array(i.split(","))[sp_type_ind] for i in data]

    # Get the z of matches
    redshifts = [np.array(i.split(","))[redshift_ind] for i in data]

    # Assign names and types to inputs
    out = refine_search(ra, dec, oid, id_out, names, types, sp_types, redshifts)

    return out


def cross_match_simbad(list_idx, list_ra, list_dec):
    """ Cross-match list with SIMBAD
    """
    # xmatch better done in list (da,dec in deg)
    matches = cross_match_alerts_raw_simbad(list_idx, list_ra, list_dec)
    xmatch_simbad_redshift = np.transpose(matches)[-1]
    xmatch_simbad_sptype = np.transpose(matches)[-2]
    xmatch_simbad_type = np.transpose(matches)[-3]
    xmatch_simbad_ctlg = np.transpose(matches)[-4]

    return (
        xmatch_simbad_redshift,
        xmatch_simbad_sptype,
        xmatch_simbad_type,
        xmatch_simbad_ctlg,
    )


def refine_search_gaia(
    ra: list,
    dec: list,
    oid: list,
    id_out: list,
    source: list,
    ragaia: list,
    decgaia: list,
    plx: list,
    plxerr: list,
    gmag: list,
    angDist: list,
) -> list:
    """ Create a final table by merging coordinates of objects found on the
    bibliographical database, with those objects which were not found.
    Parameters
    ----------
    ra: list of float
        List of RA
    dec: list of float
        List of Dec of the same size as ra.
    oid: list of str
        List of object ID (custom)
    id_out: list of str
        List of object ID returned by the xmatch with CDS
    source: list of str
        List of source ID returned by the xmatch with Gaia
    ragaia: list of float
        List of source ra returned by the xmatch with Gaia
    decgaia: list of float
        List of source ra returned by the xmatch with Gaia
    plx: list of float
        List of source parallax returned by the xmatch with Gaia
    plxerr: list of float
        List of source parallax error returned by the xmatch with Gaia
    gmag: list of float
        List of source g magnitude error returned by the xmatch with Gaia
    angDist: list of float
        List of source angular distance returned by the xmatch with Gaia
    
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (index, ra, dec, source, ragaia, decgaia, plx, plxerr, gmag, angdist).
        If the object is not found in Gaia, source, ragaia, decgaia, plx, plxerr, gmag, angdist
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    """
    out = []
    for ra_in, dec_in, id_in in zip(ra, dec, oid):
        # cast for picky Spark
        ra_in, dec_in = float(ra_in), float(dec_in)
        id_in = str(id_in)

        # Discriminate with the index
        if id_in in id_out:
            # Return the closest object in case of many
            # (smallest angular distance)
            index = id_out.index(id_in)
            source_tmp = source[index] if source[index] != "" else "Unknown"
            ragaia_tmp = ragaia[index] if ragaia[index] != "" else "Unknown"
            decgaia_tmp = decgaia[index] if decgaia[index] != "" else "Unknown"
            plx_tmp = plx[index] if plx[index] != "" else "Unknown"
            plxerr_tmp = plxerr[index] if plxerr[index] != "" else "Unknown"
            gmag_tmp = gmag[index] if gmag[index] != "" else "Unknown"
            angdist_tmp = angDist[index] if angDist[index] != "" else "Unknown"

            out.append(
                (
                    id_in,
                    ra_in,
                    dec_in,
                    source_tmp,
                    ragaia_tmp,
                    decgaia_tmp,
                    plx_tmp,
                    plxerr_tmp,
                    gmag_tmp,
                    angdist_tmp,
                )
            )

        else:
            # Mark as unknown if no match
            out.append(
                (
                    id_in,
                    ra_in,
                    dec_in,
                    "Unknown",
                    "Unknown",
                    "Unknown",
                    "Unknown",
                    "Unknown",
                    "Unknown",
                    "Unknown",
                )
            )

    return out


def cross_match_alerts_raw_gaia(oid: list, ra: list, dec: list, ctlg: str) -> list:
    """ Query the CDSXmatch service to find identified objects
    in alerts. The catalog queried is the Gaia DR database.
    Parameters
    ----------
    oid: list of str
        List containing object ids (custom)
    ra: list of float
        List containing object ra coordinates
    dec: list of float
        List containing object dec coordinates
    dec: string
        string with catalogue name
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (index, ra, dec, sourcename, gmag, gmagerr, parallax,parallaxerr,separation).
        If the object is not found in Gaia, sourcename, gmag, gmagerr, parallax,separation
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    Examples
    ----------
    >>> ra = [26.8566983, 26.24497]
    >>> dec = [-26.9677112, -26.7569436]
    >>> id = ["1", "2"]
    >>> objects = cross_match_alerts_raw(id, ra, dec)
    """
    if len(ra) == 0:
        return []

    # Catch TimeoutError and ConnectionError
    try:
        data, header = xmatch(ra, dec, oid, extcatalog=ctlg, distmaxarcsec=2)
    except (ConnectionError, TimeoutError, ValueError) as ce:
        logging.warning("XMATCH GAIA failed " + repr(ce))
        return []

    # Sometimes the service is down, but without TimeoutError or ConnectionError
    # In that case, we grab the error message from the data.
    if len(data) > 0 and "504 Gateway Time-out" in data[0]:
        msg_head = "CDS xmatch service probably down"
        msg_foot = "Check at http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
        logging.warning(msg_head)
        logging.warning(data[0])
        logging.warning(msg_foot)
        return []

    # Fields of interest (their indices in the output)
    if "source_id" not in header:
        return []

    # Fields of interest (their indices in the output)
    # sourcename, gmag, gmagerr, parallax,separation
    oid_ind = header.index("index")
    source_ind = header.index("source_id")
    ragaia_ind = header.index("ra")
    decgaia_ind = header.index("dec")
    parallax_ind = header.index("parallax")
    parallaxerr_ind = header.index("parallax_error")
    gmag_ind = header.index("phot_g_mean_mag")
    angDist_ind = header.index("angDist")

    # Get the index of matches
    id_out = [np.array(i.split(","))[oid_ind] for i in data]

    # Get the names of matches
    source = [np.array(i.split(","))[source_ind] for i in data]
    ragaia = [np.array(i.split(","))[ragaia_ind] for i in data]
    decgaia = [np.array(i.split(","))[decgaia_ind] for i in data]
    plx = [np.array(i.split(","))[parallax_ind] for i in data]
    plxerr = [np.array(i.split(","))[parallaxerr_ind] for i in data]
    gmag = [np.array(i.split(","))[gmag_ind] for i in data]
    angDist = [np.array(i.split(","))[angDist_ind] for i in data]

    out = refine_search_gaia(
        ra, dec, oid, id_out, source, ragaia, decgaia, plx, plxerr, gmag, angDist
    )

    return out


def cross_match_gaia(list_idx, list_ra, list_dec, ctlg="vizier:I/345/gaia2"):
    """ Cross-match list with SIMBAD
    """
    # xmatch better done in list (da,dec in deg)
    matches = cross_match_alerts_raw_gaia(list_idx, list_ra, list_dec, ctlg)

    xmatch_gaia_source = np.transpose(matches)[3]
    xmatch_gaia_ragaia = np.transpose(matches)[4]
    xmatch_gaia_decgaia = np.transpose(matches)[5]
    xmatch_gaia_plx = np.transpose(matches)[6]
    xmatch_gaia_plxerr = np.transpose(matches)[7]
    xmatch_gaia_gmag = np.transpose(matches)[8]
    xmatch_gaia_angdist = np.transpose(matches)[9]

    return (
        xmatch_gaia_source,
        xmatch_gaia_ragaia,
        xmatch_gaia_decgaia,
        xmatch_gaia_plx,
        xmatch_gaia_plxerr,
        xmatch_gaia_gmag,
        xmatch_gaia_angdist,
    )


def refine_search_usno(
    ra: list, dec: list, oid: list, id_out: list, source: list, angDist: list,
) -> list:
    """ Create a final table by merging coordinates of objects found on the
    bibliographical database, with those objects which were not found.
    Parameters
    ----------
    ra: list of float
        List of RA
    dec: list of float
        List of Dec of the same size as ra.
    oid: list of str
        List of object ID (custom)
    id_out: list of str
        List of object ID returned by the xmatch with CDS
    source: list of str
        List of source ID returned by the xmatch with USNO
    angDist: list of float
        List of source angular distance returned by the xmatch with USNO
    
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (index, ra, dec, source, angdist).
        If the object is not found in Gaia, source, angdist
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    """
    out = []
    for ra_in, dec_in, id_in in zip(ra, dec, oid):
        # cast for picky Spark
        ra_in, dec_in = float(ra_in), float(dec_in)
        id_in = str(id_in)

        # Discriminate with the index
        if id_in in id_out:
            # Return the closest object in case of many
            # (smallest angular distance)
            index = id_out.index(id_in)
            source_tmp = source[index] if source[index] != "" else "Unknown"
            angdist_tmp = angDist[index] if angDist[index] != "" else "Unknown"

            out.append((id_in, ra_in, dec_in, source_tmp, angdist_tmp,))

        else:
            # Mark as unknown if no match
            out.append((id_in, ra_in, dec_in, "Unknown", "Unknown",))

    return out


def cross_match_alerts_raw_usno(oid: list, ra: list, dec: list, ctlg: str) -> list:
    """ Query the CDSXmatch service to find identified objects
    in alerts. The catalog queried is the USNO A2.0 database.
    Parameters
    ----------
    oid: list of str
        List containing object ids (custom)
    ra: list of float
        List containing object ra coordinates
    dec: list of float
        List containing object dec coordinates
    dec: string
        string with catalogue name
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (index, ra, dec, sourcename ,angulardistance).
        If the object is not found in usno, index, ra, dec, sourcename ,angulardistance
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    Examples
    ----------
    >>> ra = [26.8566983, 26.24497]
    >>> dec = [-26.9677112, -26.7569436]
    >>> id = ["1", "2"]
    >>> objects = cross_match_alerts_raw_usno(id, ra, dec)
    """
    if len(ra) == 0:
        return []
    # Catch TimeoutError and ConnectionError
    try:
        data, header = xmatch(ra, dec, oid, extcatalog=ctlg, distmaxarcsec=2)
    except (ConnectionError, TimeoutError, ValueError) as ce:
        logging.warning("XMATCH USNO failed " + repr(ce))
        return []
    # Sometimes the service is down, but without TimeoutError or ConnectionError
    # In that case, we grab the error message from the data.
    if len(data) > 0 and "504 Gateway Time-out" in data[0]:
        msg_head = "CDS xmatch service probably down"
        msg_foot = "Check at http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
        logging.warning(msg_head)
        logging.warning(data[0])
        logging.warning(msg_foot)
        return []

    # Fields of interest (their indices in the output)
    if "USNO-A2.0" not in header:
        return []

    # Fields of interest (their indices in the output)
    # sourcename, gmag, gmagerr, parallax,separation
    oid_ind = header.index("index")
    source_ind = header.index("USNO-A2.0")
    angDist_ind = header.index("angDist")

    # Get the index of matches
    id_out = [np.array(i.split(","))[oid_ind] for i in data]

    # Get the names of matches
    source = [np.array(i.split(","))[source_ind] for i in data]
    angDist = [np.array(i.split(","))[angDist_ind] for i in data]

    out = refine_search_usno(ra, dec, oid, id_out, source, angDist)

    return out


def cross_match_usno(list_idx, list_ra, list_dec, ctlg="vizier:I/345/usno2"):
    """ Cross-match list with USNO A2.0 star catalogue
    """
    # xmatch better done in list (da,dec in deg)
    matches = cross_match_alerts_raw_usno(list_idx, list_ra, list_dec, ctlg)

    xmatch_usno_source = np.transpose(matches)[3]
    xmatch_usno_angdist = np.transpose(matches)[4]

    return (
        xmatch_usno_source,
        xmatch_usno_angdist,
    )


def cross_match_alerts_raw_generic(
    oid: list, ra: list, dec: list, ctlg: str, distmaxarcsec: float
) -> list:
    """ Query the CDSXmatch service to find identified objects
    in alerts. 
    Parameters
    ----------
    oid: list of str
        List containing object ids (custom)
    ra: list of float
        List containing object ra coordinates
    dec: list of float
        List containing object dec coordinates
    dec: string
        string with catalogue name
    Returns
    ----------
    out: List of Tuple
        Each tuple contains (indx, ra, dec, sourcename ,angulardistance).
        If the object is not found in usno, indx, ra, dec, sourcename ,angulardistance
        are marked as Unknown. In the case several objects match
        the centroid of the alert, only the closest is returned.
    Examples
    ----------
    >>> ra = [26.8566983, 26.24497]
    >>> dec = [-26.9677112, -26.7569436]
    >>> id = ["1", "2"]
    >>> objects = cross_match_alerts_raw_usno(id, ra, dec)
    """

    if len(ra) == 0:
        return []
    # Catch TimeoutError and ConnectionError
    try:
        data, header = xmatch(
            ra, dec, oid, extcatalog=ctlg, distmaxarcsec=distmaxarcsec
        )
    except (ConnectionError, TimeoutError, ValueError) as ce:
        logging.warning("XMATCH failed " + repr(ce))
        return []

    # Sometimes the service is down, but without TimeoutError or ConnectionError
    # In that case, we grab the error message from the data.
    if len(data) > 0 and "504 Gateway Time-out" in data[0]:
        msg_head = "CDS xmatch service probably down"
        msg_foot = "Check at http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
        logging.warning(msg_head)
        logging.warning(data[0])
        logging.warning(msg_foot)
        return []

    df_out_tmp = pd.DataFrame()
    df_out_tmp["objectId"] = oid
    df_out_tmp["objectId"] = df_out_tmp["objectId"].astype(str)
    df_out_tmp["ra"] = ra
    df_out_tmp["dec"] = dec

    if len(data) == 0:
        print("No match found")
        return df_out_tmp

    data = [x.split(",") for x in data]
    df_search_out = pd.DataFrame(data=np.array(data), columns=header)
    if "angDist" not in df_search_out.keys():
        print("Xmatch failure")
        raise Exception
    else:
        df_search_out["angDist"] = df_search_out["angDist"].astype(float)
        df_search_out = df_search_out.rename(columns={"objectId": "idx_to_match"})
        df_search_out_tmp = df_search_out.sort_values("angDist", ascending=True)
        df_search_out_tmp = df_search_out_tmp.groupby("idx_to_match").first()
        df_search_out_tmp = df_search_out_tmp.rename(
            columns={"ra": "ra_out", "dec": "dec_out"}
        )
        df_search_out_tmp["objectId"] = df_search_out_tmp.index

        df_out = pd.merge(df_out_tmp, df_search_out_tmp, on="objectId", how="left")
        df_out = df_out.fillna("Unknown")
        df_out = df_out.drop(["ra_in", "dec_in"], axis=1)

        return df_out
