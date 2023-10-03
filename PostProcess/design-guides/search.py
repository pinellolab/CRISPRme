"""
"""

from encoder import encode_pam


def search_pam(pam: str, genome: str, mm: int) -> None:
    pam_bits = encode_pam(pam)  # encode PAM in bits for efficient search
