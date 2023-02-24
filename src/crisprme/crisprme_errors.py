"""
"""


class CRISRPRmeError(Exception):
    pass


class CompleteSearchError(CRISRPRmeError):
    pass


class CRISPRitzError(CRISRPRmeError):
    pass


class VerbosityHandlerError(CRISRPRmeError):
    pass


class TSTError(CRISRPRmeError):
    pass
