from enum import IntFlag, auto

class RowFlags(IntFlag):
    BAD_DATA = auto()                # 1
    BAD_DATA_LEN = auto()
    NO_HK = auto()
    BAD_PHASE = auto()
    MISSING_INT = auto()             # 5
    MIXER_MISPUMPED = auto()
    MIXER_UNPUMPED = auto()
    UNSTABLE_POINTING = auto()
    UNSTABLE_TIMEBASE = auto()
    PROCESSOR_ERROR = auto()         # 10
    MIXERCURRENT_MISMATCH = auto()
    ALL_NAN = auto()
    BAD_BASELINE = auto()
    UNABLE_TO_PROCESS = auto()
    LO_SYNTH_UNLOCKED = auto()       # 15
    RINGING_BIT0 = auto()
    RINGING_BIT1 = auto()
    QUALITY_BIT0 = auto()
    QUALITY_BIT1 = auto()
    QUALITY_BIT2 = auto()            # 20
    DAC_CAL_FIXED = 1 << 29
    DAC_CAL_FAKED = 1 << 30

class ChanFlags(IntFlag):
    BADCHAN = auto()                 # 1
    INTERPOLATED = auto()
    BADCAL = auto()
    SPUR_CANDIDATE = auto()
    SPUR = auto()                    # 5
    FRINGE = auto()
    LINE = auto()
    VARIABLE_SPUR = auto()
    OOB = auto()
    WINDOW = auto()                  # 10
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