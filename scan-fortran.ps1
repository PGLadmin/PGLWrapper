<#
  scan-fortran.ps1
  Scans Fortran sources under the script folder (or a provided -Root).
  The script reports all uses of modules and common blocks.
  The script will default to scanning all subfolders.
  To limit to a single subfolder, move to the desired subfolder and run:
      PS> .\scan-fortran.ps1 -Root '.'
#>

param(
  [string]$Root = $PSScriptRoot   # default to the folder where this script lives
)

# Resolve the root to an absolute path
$Root = (Resolve-Path -Path $Root).Path

# Create timestamped output file
$timestamp  = Get-Date -Format "yyyyMMdd_HHmmss"
$OutputFile = Join-Path $PSScriptRoot "scan-output_$timestamp.txt"

# Run the scan in a script block and tee output to file + console
& {
    Push-Location $Root

    Write-Host "Scanning Fortran sources under $Root" -ForegroundColor Cyan
	Write-Output "Scanning Fortran sources under $Root"
    ""

    # File search pattern
    $patterns = @("*.f90","*.for","*.f")
    $notpatterns = @("build","*genmod.f90")

    # Helper to get files (handles spaces and long paths)
    $files = Get-ChildItem -Path $Root -Recurse -File -Include $patterns -ErrorAction SilentlyContinue | Where-Object FullName -notmatch build

    if (-not $files) {
        Write-Host "No Fortran files found under $Root" -ForegroundColor Red
		Write-Output "No Fortran files found under $Root"
        Pop-Location
        return
    }

    # Module definitions (case-insensitive)
    Write-Host "`nModule definitions (file : line):" -ForegroundColor Yellow
	Write-Output "`nModule definitions (file : line):"
    $files |
        Select-String -Pattern '(?i)^\s*module\b' |
        Format-Table Path, LineNumber, Line -AutoSize

    # Module uses
    Write-Host "`nModule uses (file : line):" -ForegroundColor Yellow
	Write-Output "`nModule uses (file : line):"
    $files |
        Select-String -Pattern '(?i)^\s*use\b' |
        Format-Table Path, LineNumber, Line -AutoSize

    # COMMON blocks (case-insensitive)
    Write-Host "`nCOMMON blocks (file : line):" -ForegroundColor Yellow
	Write-Output "`nCOMMON blocks (file : line):"
    $files |
        Select-String -Pattern '(?i)common\s*/' |
        Format-Table Path, LineNumber, Line -AutoSize

    # Unique file lists (leaf names)
    Write-Host "`nUnique files that define modules:" -ForegroundColor Yellow
	Write-Output "`nUnique files that define modules:"
    $files |
        Select-String -Pattern '(?i)^\s*module\b' |
        Select-Object -ExpandProperty Path |
        Sort-Object -Unique |
        ForEach-Object { Split-Path $_ -Leaf }

    Write-Host "`nUnique files that contain COMMON blocks:" -ForegroundColor Yellow
	Write-Output "`nUnique files that contain COMMON blocks:"
    $files |
        Select-String -Pattern '(?i)common\s*/' |
        Select-Object -ExpandProperty Path |
        Sort-Object -Unique |
        ForEach-Object { Split-Path $_ -Leaf }

    Pop-Location

} | Tee-Object -FilePath $OutputFile
