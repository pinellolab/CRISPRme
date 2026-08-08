#!/usr/bin/env python3
"""Playwright end-to-end smoke test for the CRISPRme Dash web interface.

Guards the Dash 1.x -> 2.x migration (2.2.0): the migration's regressions all
manifested as *blank pages* (the removed ``dbc.FormGroup`` blanked the homepage;
a dangling callback output blanked every page). This test therefore drives a
real browser against a running web app and asserts that each page renders
non-blank, shows the shared navbar, and produces no uncaught JavaScript errors
or Dash error overlay.

It is intentionally dependency-light: a standalone script (no pytest), driven by
``playwright.sync_api``, that exits non-zero on the first batch of failures so it
can be run directly in CI:

    CRISPRME_WEB_URL=http://localhost:8080 python test/web/test_web_e2e.py

The web app must already be serving at ``CRISPRME_WEB_URL`` (default
http://localhost:8080).
"""

import os
import sys

from playwright.sync_api import sync_playwright

BASE = os.environ.get("CRISPRME_WEB_URL", "http://localhost:8080").rstrip("/")

# Pages served by the Dash router (see index.py + navbar_creation.py). Each must
# render non-blank with the shared navbar present.
PAGES = ["/index", "/history", "/user-guide", "/contacts", "/load"]

# Selectors that indicate a Dash server/callback error was rendered into the DOM.
DASH_ERROR_SELECTORS = "._dash-error, .dash-fe-error__title, .dash-debug-error"


def check_page(page, path, failures):
    url = f"{BASE}{path}"
    try:
        page.goto(url, wait_until="networkidle", timeout=45000)
    except Exception as e:  # noqa: BLE001 - report, don't crash the whole run
        failures.append(f"{path}: navigation failed ({e})")
        return
    body = (page.inner_text("body") or "").strip()
    # not blank — the core regression the Dash 2.x migration had to fix
    if len(body) < 20:
        failures.append(f"{path}: page is blank (body len {len(body)})")
    # shared navbar rendered (present in the base layout on every page)
    if page.locator("#navbar-toggler").count() == 0:
        failures.append(f"{path}: navbar not rendered")
    # no Dash error overlay
    if page.locator(DASH_ERROR_SELECTORS).count() > 0:
        failures.append(f"{path}: Dash error overlay present")


def main():
    failures = []
    js_errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("pageerror", lambda e: js_errors.append(str(e)))

        # 1. homepage renders with its key controls (Submit + brand)
        try:
            page.goto(BASE + "/index", wait_until="networkidle", timeout=45000)
        except Exception as e:  # noqa: BLE001
            failures.append(f"/index: initial load failed ({e})")
        home = (page.inner_text("body") or "")
        if "CRISPRme" not in home:
            failures.append("/index: 'CRISPRme' text not found (homepage blank?)")
        if page.get_by_text("Submit", exact=False).count() == 0:
            failures.append("/index: no Submit control (main form did not render)")

        # 2. every routed page renders non-blank with the navbar
        for path in PAGES:
            check_page(page, path, failures)

        # 3. no uncaught JS errors anywhere in the walk
        if js_errors:
            failures.append(f"uncaught JS error(s): {js_errors}")

        browser.close()

    if failures:
        print("WEB E2E FAILED:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print(f"WEB E2E PASSED — {len(PAGES)} pages rendered, no JS errors")


if __name__ == "__main__":
    main()
