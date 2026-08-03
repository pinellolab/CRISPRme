"""Lightweight, non-fatal system memory check for CRISPRme.

CRISPRme's complete-search / complete-test pipelines are memory hungry.
When run inside Docker Desktop the container often inherits a small default
memory limit (e.g. 2-8 GB), which leads to confusing OOM failures midway
through a search.  This module prints a *non-fatal* warning to stderr when
the detected total RAM is below a recommended threshold, pointing Docker
Desktop users to raise their memory allocation.

The check is intentionally defensive: it never raises, and if total memory
cannot be determined it silently returns.

Environment overrides:
    CRISPRME_SKIP_MEM_CHECK=1   skip the check entirely
    CRISPRME_MIN_GB=<number>    override the warning threshold (in GB)
"""

import os
import sys

# Recommended minimum total RAM (in GB) for CRISPRme's heavy pipelines.
DEFAULT_MIN_GB = 32


def _total_memory_bytes():
    """Return total physical RAM in bytes, or None if it cannot be detected.

    Tries /proc/meminfo first (Linux, incl. Docker containers) and falls back
    to sysconf (macOS / other Unix). Never raises.
    """
    # Linux (and Linux Docker containers): parse /proc/meminfo.
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    # Format: "MemTotal:       16318360 kB"
                    parts = line.split()
                    return int(parts[1]) * 1024  # kB -> bytes
    except (OSError, ValueError, IndexError):
        pass
    # Fallback (macOS / other Unix): page size * number of physical pages.
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        phys_pages = os.sysconf("SC_PHYS_PAGES")
        if page_size > 0 and phys_pages > 0:
            return page_size * phys_pages
    except (ValueError, OSError, AttributeError):
        pass
    return None


def warn_low_memory(min_gb=DEFAULT_MIN_GB):
    """Print a non-fatal stderr WARNING when total RAM is below ``min_gb``.

    Honors CRISPRME_SKIP_MEM_CHECK (skip) and CRISPRME_MIN_GB (override
    threshold). Silently returns when memory cannot be detected. Never raises.
    """
    try:
        if os.environ.get("CRISPRME_SKIP_MEM_CHECK") == "1":
            return
        override = os.environ.get("CRISPRME_MIN_GB")
        if override:
            try:
                min_gb = float(override)
            except ValueError:
                pass
        total_bytes = _total_memory_bytes()
        if total_bytes is None:  # undetectable -> stay silent
            return
        total_gb = total_bytes / (1024 ** 3)
        if total_gb < min_gb:
            sys.stderr.write(
                "WARNING: CRISPRme detected only {:.1f} GB of total system "
                "memory, which is below the recommended {:.0f} GB.\n"
                "         Heavy searches may run out of memory and fail.\n"
                "         If you are running via Docker Desktop, raise the "
                "memory limit in\n"
                "         Settings -> Resources -> Memory to at least "
                "{:.0f} GB and restart the container.\n".format(
                    total_gb, float(min_gb), float(min_gb)
                )
            )
    except Exception:
        # Absolutely never let a memory check break the pipeline.
        return
