import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord


def convert_coords(ra, dec):
    try:
        if "h" in str(ra):
            coord = SkyCoord(ra, dec, frame="icrs")
        elif ":" in str(ra) or " " in str(ra):
            coord = SkyCoord(ra, dec, frame="icrs", unit=(u.hourangle, u.deg))
        else:
            coord = SkyCoord(ra, dec, frame="icrs", unit="deg")
    except ValueError as e:
        print("Format not supported")
    return coord


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert coordinates")

    parser.add_argument("--ra", default="ra", help="ra in any format")
    parser.add_argument("--dec", default="dec", help="dec in any format")
    parser.add_argument(
        "--format", choices=["deg", "hms"], default="deg",
        help="output format: deg (decimal degrees) or hms (hh:mm:ss.s dd:mm:ss.s)",
    )
    args = parser.parse_args()

    coord = convert_coords(args.ra, args.dec)

    if args.format == "hms":
        ra_str = coord.ra.to_string(unit=u.hourangle, sep=":", precision=1, pad=True)
        dec_str = coord.dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True)
        print(ra_str, dec_str)
    else:
        print(coord.ra.deg, coord.dec.deg)
