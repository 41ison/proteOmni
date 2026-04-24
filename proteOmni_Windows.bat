@echo off
REM ── proteOmni Launcher (Windows) ────────────────────────────────────────
REM Double-click this file to start proteOmni in your default browser.
REM Requirements: R must be installed and Rscript must be in PATH.

REM Enable UTF-8 so emoji render in Windows Terminal / PowerShell
chcp 65001 >nul

echo.
echo   ╔══════════════════════════════════════════════════╗
echo          proteOmni - Proteomics QC Dashboard
echo   ╚══════════════════════════════════════════════════╝
echo.
echo   🚀  Starting proteOmni...
echo   🌐  Opening in your default browser...
echo.

cd /d "%~dp0"
Rscript -e "shiny::runApp('proteOmni.r', launch.browser = TRUE)"
pause
