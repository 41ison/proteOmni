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
# run.R installs any missing dependencies (including shiny itself) before
# starting the app. Do NOT call shiny::runApp() here: on a fresh R installation
# shiny does not exist yet.
Rscript run.R
