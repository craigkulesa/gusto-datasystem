from enum import IntFlag, auto

class RowFlags(IntFlag):
    BAD_DATA = auto()
    BAD_DATA_LEN = auto()
    MISSING_ACS = auto()
    NO_HK = auto()
    HK_MISALIGNED = auto()
    BAD_PHASE = auto()
    MISSING_INT = auto()
    MIXER_MISPUMPED = auto()
    UNSTABLE_POINTING = auto()
    UNSTABLE_TIMEBASE = auto()
    PROCESSOR_ERROR = auto()
    MIXERCURRENT_MISMATCH = auto()
    ALL_NAN = auto()
    BAD_BASELINE = auto()
    UNABLE_TO_PROCESS = auto()
    RINGING_BIT0 = 1 << 26
    RINGING_BIT1 = 1 << 27
    QUALITY_BIT0 = 1 << 28
    QUALITY_BIT1 = 1 << 29
    QUALITY_BIT2 = 1 << 30
    IGNORE_OBSLOG = 1 << 31

    
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
    MISSING_HOT = auto()
    MISSING_REF = auto()
    NOREFS = auto()
    NODATA = auto()
    
