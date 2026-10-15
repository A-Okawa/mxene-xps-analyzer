"""Browser E2E run + screenshots of each workflow page (synthetic data only).

Requires a running app:  .venv/bin/streamlit run app.py --server.headless true --server.port 8599
Then:                    .venv/bin/python tools/e2e_screenshots.py [out_dir]
Uses the system Google Chrome via Playwright (dev-only dependency).
"""
from __future__ import annotations

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "docs" / "screenshots"
URL = "http://localhost:8599"


def wait_idle(page, ms=1500):
    page.wait_for_timeout(400)
    page.wait_for_function("() => !document.querySelector('[data-testid=\"stStatusWidget\"]')", timeout=60000)
    page.wait_for_timeout(ms)


def shot(page, name):  # noqa: D401
    h = page.evaluate("() => Math.max(...[...document.querySelectorAll('[data-testid=stMain], section, div')]"
                      ".map(e => e.scrollHeight))")
    page.set_viewport_size({"width": 1500, "height": min(int(h) + 50, 12000)})
    page.wait_for_timeout(1200)
    page.screenshot(path=str(OUT / name), full_page=True)
    page.set_viewport_size({"width": 1500, "height": 1100})
    page.wait_for_timeout(300)


def tab(page, label):
    page.get_by_role("tab", name=label).click()
    wait_idle(page)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = []
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", headless=True)
        page = b.new_page(viewport={"width": 1500, "height": 1100})
        page.goto(URL)
        page.wait_for_selector("text=MAX/MXene XPS Analyzer", timeout=60000)
        wait_idle(page, 2500)
        shot(page, "1_request_generator.png")
        log.append("request tab rendered")

        tab(page, "2. Peak table & quantification")
        page.locator("[data-testid='stFileUploader']").first.locator("input[type=file]").set_input_files(
            str(FIX / "synthetic_NbTi_corrected_area_peak_table.csv"))
        wait_idle(page, 2500)
        page.get_by_role("combobox", name="Area input mode").click()
        page.keyboard.press("ArrowDown")
        page.wait_for_timeout(500)
        page.locator("[role=option]", has_text="corrected_area").first.click()
        wait_idle(page)
        page.get_by_role("button", name="Calculate", exact=True).click()
        wait_idle(page, 2500)
        body = page.inner_text("body")
        assert "apparent_M_deficiency_delta" in body, "quantification result missing"
        assert "unresolved mixed-metal" in body
        log.append("quantification ran (corrected_area); unresolved-O warning shown")
        shot(page, "2_quantification.png")

        tab(page, "3. Provider ↔ Python validation")
        ups = page.locator("[data-testid='stFileUploader'] input[type=file]")
        # uploaders on this tab: raw spectrum, provider table (first uploader of tab 2 is hidden in another tab panel)
        visible = [i for i in range(ups.count()) if ups.nth(i).locator("xpath=ancestor::div[@role='tabpanel']").first.is_visible()]
        ups.nth(visible[0]).set_input_files(str(FIX / "synthetic_Ti2p_raw_desc_noheader.txt"))
        wait_idle(page)
        ups.nth(visible[1]).set_input_files(str(FIX / "synthetic_Ti2p_provider_fit.csv"))
        wait_idle(page)
        page.get_by_role("button", name="Run Python fit").click()
        wait_idle(page, 4000)
        body = page.inner_text("body")
        assert "Role-grouped QC" in body, "validation tables missing"
        log.append("validation ran on descending/no-header raw file with provider table")
        shot(page, "3_validation.png")

        tab(page, "4. Knowledge base")
        shot(page, "4_knowledge_base.png")
        tab(page, "5. Report & export")
        shot(page, "5_report_export.png")
        log.append("KB and report tabs rendered")

        # Legacy Ti3C2Tx model path (synthetic, made-up areas – not workbook data)
        page.get_by_role("combobox", name="Nominal structure").click()
        page.keyboard.press("ArrowDown")
        page.wait_for_timeout(500)
        page.locator("[role=option]", has_text="Ti3C2Tx").first.click()
        wait_idle(page, 2500)
        tab(page, "2. Peak table & quantification")
        page.get_by_text("Legacy Ti3C2Tx workbook model").click()
        wait_idle(page, 2000)
        ups = page.locator("[data-testid='stFileUploader'] input[type=file]")
        visible = [i for i in range(ups.count()) if ups.nth(i).locator("xpath=ancestor::div[@role='tabpanel']").first.is_visible()]
        ups.nth(visible[0]).set_input_files(str(FIX / "synthetic_legacy_Ti3C2_areas.csv"))
        wait_idle(page, 2000)
        page.get_by_role("button", name="Calculate (legacy model)").click()
        wait_idle(page, 2500)
        body = page.inner_text("body")
        if not ("Download quantification results" in body and "Legacy/validated" in body):
            shot(page, "_debug_legacy.png")
            raise AssertionError("legacy results missing")
        log.append("legacy Ti3C2Tx model ran")
        shot(page, "6_legacy_model.png")
        b.close()
    print("\n".join(log))
    print("screenshots in", OUT)


if __name__ == "__main__":
    main()
