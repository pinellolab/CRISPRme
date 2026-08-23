# CRISPRme+ one-line installer (Windows).
#   irm https://pinellolab.github.io/CRISPRme/install.ps1 | iex
#
# Installs a clickable "CRISPRme" app (a small window with buttons) that manages the
# Docker container, downloads reference data, and opens the web interface — so after
# this one command the user never needs the terminal again. Docker Desktop required.
$ErrorActionPreference = 'Stop'
$Image   = 'pinellolab/crisprme:v2.4.0'
$DataDir = Join-Path $env:USERPROFILE 'CRISPRme'
$AppDir  = Join-Path $env:LOCALAPPDATA 'CRISPRme'

Write-Host "=== CRISPRme+ installer ===" -ForegroundColor Cyan

# ---- 1. Docker check -------------------------------------------------------
if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
  Write-Host "Docker is not installed. CRISPRme runs inside Docker." -ForegroundColor Yellow
  Write-Host "Please install Docker Desktop, then re-run this installer:"
  Write-Host "  https://www.docker.com/products/docker-desktop/" -ForegroundColor Blue
  Start-Process "https://www.docker.com/products/docker-desktop/"
  return
}
Write-Host "Docker found." -ForegroundColor Green

# ---- 2. Working dir + compose file ----------------------------------------
New-Item -ItemType Directory -Force -Path (Join-Path $DataDir 'crisprme-data') | Out-Null
@"
# Managed by the CRISPRme app. Data persists in ./crisprme-data.
services:
  crisprme:
    image: $Image
    working_dir: /DATA
    volumes:
      - ./crisprme-data:/DATA
    ports:
      - "8080:8080"
    command: crisprme.py web-interface
    stdin_open: true
    tty: true
    restart: "no"
"@ | Set-Content -Encoding ascii (Join-Path $DataDir 'docker-compose.yml')
Write-Host "Data folder ready: $DataDir" -ForegroundColor Green

# ---- 3. Write the clickable app (a WinForms window) ------------------------
New-Item -ItemType Directory -Force -Path $AppDir | Out-Null
$AppScript = Join-Path $AppDir 'CRISPRme.ps1'
@'
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
$DataDir = Join-Path $env:USERPROFILE 'CRISPRme'
$Compose = Join-Path $DataDir 'docker-compose.yml'
$Image   = 'pinellolab/crisprme:v2.4.0'

function Has-Data { Test-Path (Join-Path $DataDir 'crisprme-data\Genomes') }
function Docker-Up { try { docker info *> $null; return $true } catch { return $false } }
function Run-InConsole($title, $cmd) {
  Start-Process powershell -ArgumentList @('-NoExit','-Command',
    "Write-Host '== $title =='; Write-Host 'You can close this window when it finishes.'; $cmd; Write-Host ''; Write-Host 'DONE - you can close this window.'")
}

$form = New-Object System.Windows.Forms.Form
$form.Text = 'CRISPRme+'
$form.Size = New-Object System.Drawing.Size(430,340)
$form.StartPosition = 'CenterScreen'
$form.FormBorderStyle = 'FixedSingle'; $form.MaximizeBox = $false

$status = New-Object System.Windows.Forms.Label
$status.SetBounds(20,15,390,40); $status.Text = 'CRISPRme+ off-target analysis'
$form.Controls.Add($status)

function Mk($text,$y,$action) {
  $b = New-Object System.Windows.Forms.Button
  $b.SetBounds(20,$y,390,38); $b.Text = $text; $b.Add_Click($action)
  $form.Controls.Add($b); return $b
}
Mk 'Open CRISPRme (start + browser)' 60 {
  if (-not (Docker-Up)) { [System.Windows.Forms.MessageBox]::Show('Docker Desktop is not running. Open it, wait, then try again.'); return }
  if (-not (Has-Data))  { [System.Windows.Forms.MessageBox]::Show('No reference data yet. Choose "Download reference data" first.'); return }
  docker compose -f $Compose up -d
  Start-Sleep -Seconds 3; Start-Process 'http://localhost:8080'
} | Out-Null
Mk 'Download reference data (~25 GB, once)' 104 {
  Run-InConsole 'Downloading reference data' "docker run --rm -v `"${DataDir}\crisprme-data:/DATA`" -w /DATA $Image crisprme.py download --what all --path /DATA"
} | Out-Null
Mk 'Download variant index (advanced, 64 GB RAM)' 148 {
  Run-InConsole 'Downloading variant index' "docker run --rm -v `"${DataDir}\crisprme-data:/DATA`" -w /DATA $Image crisprme.py download --what index --index-name NRG_3_hg38-dictless+hg38_1000G_HGDP --path /DATA"
} | Out-Null
Mk 'Update CRISPRme' 192 { Run-InConsole 'Updating CRISPRme' "docker pull $Image" } | Out-Null
Mk 'Stop CRISPRme' 236 { try { docker compose -f $Compose down } catch {}; [System.Windows.Forms.MessageBox]::Show('CRISPRme stopped.') } | Out-Null

if (Has-Data) { $status.Text = 'CRISPRme+  -  reference data: ready' } else { $status.Text = 'CRISPRme+  -  reference data: MISSING (download it first)' }
[void]$form.ShowDialog()
'@ | Set-Content -Encoding UTF8 $AppScript

# ---- 4. Shortcuts (Desktop + Start Menu) ----------------------------------
$ws = New-Object -ComObject WScript.Shell
foreach ($dir in @([Environment]::GetFolderPath('Desktop'),
                   (Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs'))) {
  $lnk = $ws.CreateShortcut((Join-Path $dir 'CRISPRme.lnk'))
  $lnk.TargetPath = 'powershell.exe'
  $lnk.Arguments  = "-NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -File `"$AppScript`""
  $lnk.WorkingDirectory = $DataDir
  $lnk.IconLocation = 'powershell.exe,0'
  $lnk.Save()
}

Write-Host ""
Write-Host "Done! Open 'CRISPRme' from your Desktop or Start Menu." -ForegroundColor Green
Write-Host "First time: click 'Download reference data', then 'Open CRISPRme'."
