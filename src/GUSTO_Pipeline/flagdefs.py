from enum import IntFlag, auto

class RowFlags(IntFlag):
    BAD_DATA = auto()
    BAD_DATA_LEN = auto()
    NO_HK = auto()
    BAD_PHASE = auto()
    MISSING_INT = auto()
    MIXER_MISPUMPED = auto()
    MIXER_UNPUMPED = auto()
    UNSTABLE_POINTING = auto()
    UNSTABLE_TIMEBASE = auto()
    PROCESSOR_ERROR = auto()
    MIXERCURRENT_MISMATCH = auto()
    ALL_NAN = auto()
    BAD_BASELINE = auto()
    UNABLE_TO_PROCESS = auto()
    LO_SYNTH_UNLOCKED = auto()
    RINGING_BIT0 = auto()
    RINGING_BIT1 = auto()
    QUALITY_BIT0 = auto()
    QUALITY_BIT1 = auto()
    QUALITY_BIT2 = auto()
    DAC_CAL_FIXED = 1 << 29
    DAC_CAL_FAKED = 1 << 30

class ChanFlags(IntFlag):
    BADCHAN = auto()
    INTERPOLATED = auto()
    BADCAL = auto()
    SPUR_CANDIDATE = auto()
    SPUR = auto()
    FRINGE = auto()
    LINE = auto()
    VARIABLE_SPUR = auto()
    OOB = auto()
    WINDOW = auto()
    IGNORE = auto()

class SeqFlags(IntFlag):
    MISSING_HK = auto()
    MAYBE_HUNG_LO = auto()
    MISSING_HOT = auto()
    MISSING_LEADING_REF = auto()
    MISSING_TRAILING_REF = auto()
    NOREFS = auto()
    NODATA = auto()

B2Unlocked = [(13756, 14114)]
B1Unlocked = [(13756, 14114)]


def string_to_enum_combination(istring):
    """Converts a string of color names to an enum combination.
    
    """
    if type(istring)==type('m'):
        inames = istring.split()
    else:
        inames = str(istring).split()
    print(inames)
    cflags = RowFlags(0)  # Initialize with no value
    for name in inames:
        if '|' in name:
            continue
        else:
            try:
                flag = RowFlags[name.replace("'","").replace('"','').split('.')[1].upper()]
                cflags |= flag
            except KeyError:
                raise ValueError(f"Invalid flag: {name}")
    return cflags


def anaFlagString(aa):
    r"""Function .

    Parameters
    ----------
    aa : int or string
            rowflagfilter integer or string representation of allowed flags

    Examples
    --------
    aa = 'RowFlags.NO_HK | RowFlags.MISSING_INT | RowFlags.MIXER_UNPUMPED'
    enum_combination = anaFlagString(aa)
    
    print(repr(enum_combination), type(enum_combination))
    #<RowFlags.NO_HK|MISSING_INT|MIXER_UNPUMPED: 84> <flag 'RowFlags'>
    
    """
    if aa.isnumeric():
        return RowFlags(int(aa))
    else:
        aa = aa.replace("|","")
        return string_to_enum_combination(aa)

