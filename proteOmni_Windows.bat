@echo off
REM ── proteOmni Launcher (Windows) ────────────────────────────────────────
REM Double-click this file to start proteOmni in your default browser.
REM Requirements: R must be installed and Rscript must be in PATH.

REM Enable UTF-8 so emoji render in Windows Terminal / PowerShell
chcp 65001 >nul

echo.
echo   +--------------------------------------------------+
echo          proteOmni - Proteomics QC Dashboard
echo   +--------------------------------------------------+
echo.
echo   ^>^>  Starting proteOmni...
echo   ^>^>  Opening in your default browser...
echo.

cd /d "%~dp0"
REM run.R installs any missing dependencies (including shiny itself) before
REM starting the app. Do NOT call shiny::runApp() here: on a fresh R
REM installation shiny does not exist yet.
Rscript run.R
pause
