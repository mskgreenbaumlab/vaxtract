#!/usr/bin/env python3
"""Back-compat shim — the renderer lives in ``vaxtract.report``.

    python cancervac_packet/make_report.py EXTRACTED.json [OUT.html]
    vaxtract-report EXTRACTED.json [OUT.html]
"""
from vaxtract.report import (  # noqa: F401
    DASH,
    SOFT_OVERRIDES,
    TCR_CLONE_CAP,
    build_html,
    chain_cdr3,
    count_needs_review,
    esc,
    load_record,
    mag_str,
    main,
    responders,
    tcr_target,
)

if __name__ == "__main__":
    raise SystemExit(main())
