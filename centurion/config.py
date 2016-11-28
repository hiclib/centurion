import ConfigParser
import os


def get_default_options():
    options = {"output_name": "centromeres_calls.txt",
               "resolution": 10000,
               "counts": "data/counts.npy",
               "lengths": "data/lengths.npy",
               "seed": 0,
               }
    return options


def parse(filename=None):
    """
    Parses a configuration file.

    Parameters
    ----------
    filename : str, optional, default: None
        If a filename is provided, reads the configuration file, and returns
        the options. If None, returns the default options.

    Returns
    -------
    options : dict
    """
    options = get_default_options()
    if filename is None:
        return options

    if not os.path.exists(filename):
        raise IOError("File %s doesn't existe" % filename)

    config = ConfigParser.ConfigParser()
    config.readfp(open(filename))
    for key in options.iterkeys():
        try:
            options[key] = type(options[key])(
                config.get("all", key))
        except ConfigParser.NoOptionError:
            pass
    return options
