from __future__ import annotations

import pooch

__all__ = [
    "fetch_nisar_frame_to_bounds_file",
]

NISAR_FRAME_DB_VERSION = "0.1.0"
NISAR_FRAME_TO_BOUNDS_FILENAME = (
    f"opera-nisar-disp-{NISAR_FRAME_DB_VERSION}-frame-to-bounds.json"
)

NISAR_POOCH = pooch.create(
    path=pooch.os_cache("opera_utils"),
    # TODO: update to a release URL once the file is published.
    # Must be a raw-content URL: github.com/.../tree/... serves the HTML file
    # viewer, so pooch downloads a web page and rejects it against the hash
    # below ("does not match the known hash"), with no hint that the URL is at
    # fault. The hash here is correct for the JSON itself.
    base_url=(
        "https://raw.githubusercontent.com/opera-adt/disp-nisar/main/"
        "configs/static_ancillary_files/"
    ),
    version=NISAR_FRAME_DB_VERSION,
    version_dev="main",
    env="OPERA_UTILS_DATA_DIR",
    # $ shasum -a 256 opera-nisar-disp-0.1.0-frame-to-bounds.json
    # f9f2e64f34cedadb9a35d7a792990f684f153eb5463a80d5b26982a942ea1a03
    registry={
        NISAR_FRAME_TO_BOUNDS_FILENAME: (
            "f9f2e64f34cedadb9a35d7a792990f684f153eb5463a80d5b26982a942ea1a03"
        ),
    },
)


def fetch_nisar_frame_to_bounds_file() -> str:
    """Get the NISAR frame-to-bounds mapping file."""
    return NISAR_POOCH.fetch(NISAR_FRAME_TO_BOUNDS_FILENAME)
