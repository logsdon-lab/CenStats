from enum import StrEnum


class Orientation(StrEnum):
    """
    Orientation of a contig.

    `| <- Forward | Reverse -> |`
    """

    Forward = "fwd"
    Reverse = "rev"
