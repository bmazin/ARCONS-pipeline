import tables

timeMaskReasonList = []
timeMaskReasonList.append("unknown");
timeMaskReasonList.append("Flash in r0");
timeMaskReasonList.append("Flash in r1");
timeMaskReasonList.append("Flash in r2");
timeMaskReasonList.append("Flash in r3");
timeMaskReasonList.append("Flash in r4");
timeMaskReasonList.append("Flash in r5");
timeMaskReasonList.append("Flash in r6");
timeMaskReasonList.append("Flash in r7");
timeMaskReasonList.append("Merged Flash");
timeMaskReasonList.append("cosmic");
timeMaskReasonList.append("poofing");

timeMaskReason = tables.Enum(timeMaskReasonList)


class TimeMask(tables.IsDescription):
    """The pytables derived class that stores time intervals to be masked"""
    tBegin = tables.UInt32Col() # beginning time of this mask (clock ticks)
    tEnd   = tables.UInt32Col() # ending time of this mask (clock ticks)
    reason = tables.EnumCol(timeMaskReason, "unknown", base='uint8')