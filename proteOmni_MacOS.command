#!/bin/bash
# ── proteOmni Launcher (macOS / Linux) ──────────────────────────────────────
# Double-click this file to start proteOmni in your default browser.
# Requirements: R must be installed and available in PATH.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
echo ""
echo "  ╔══════════════════════════════════════════════════╗"
echo "  ║     proteOmni — Proteomics QC Dashboard          ║"
echo "  ╚══════════════════════════════════════════════════╝"
echo ""
echo "  🚀  Starting proteOmni from: $SCRIPT_DIR"
echo "  🌐  Opening in your default browser..."
echo ""

cd "$SCRIPT_DIR"
Rscript -e "shiny::runApp('proteOmni.r', launch.browser = TRUE)"
