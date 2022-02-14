import argparse

from . import idr_tools, misc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="My tools", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='a command. Add -h for command help', dest="command")
    
    parser_peak_to_narrowpeak = subparsers.add_parser('peaks_to_narrowpeak', help='Homer .peaks -> ENCODE narrowPeak',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_peak_to_narrowpeak.add_argument("peaks", metavar="PEAKS", type=str, help="peaks file")
    parser_peak_to_narrowpeak.add_argument("narrowPeak", metavar="NARROWPEAK", type=str, nargs="?", help="narrowPeak file")

    args = parser.parse_args()
    if args.command == "peaks_to_narrowpeak":
        idr.peaks_to_narrowpeak(args.peaks, args.narrowPeak)
    else:
        parser.print_help()
