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
    DAC_CAL_FIXED = 1 << 25
    RINGING_BIT0 = 1 << 26
    RINGING_BIT1 = 1 << 27
    QUALITY_BIT0 = 1 << 28
    QUALITY_BIT1 = 1 << 29
    QUALITY_BIT2 = 1 << 30
    DAC_CAL_FAKED = 1 << 31

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
    WINDOW = 1 << 14
    IGNORE = 1 << 15

class SeqFlags(IntFlag):
    MISSING_HK = auto()
    MAYBE_HUNG_LO = auto()
    MISSING_HOT = auto()
    MISSING_REF = auto()
    NOREFS = auto()
    NODATA = auto()
    
B2Unlocked = [(13756, 14114)]



def getFlagNames(flag_value, enum_type):
    """Returns a list of flag names that compose the given value."""
    if not isinstance(flag_value, enum_type):
        flag_value = enum_type(flag_value)

    names = []
    for member in enum_type:
        if member.value & flag_value.value:
            names.append(member.name)
    return names

def getRowFlagNames(flag_value):
    return getFlagNames(flag_value, RowFlags)


