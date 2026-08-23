# CRISPRme+ one-line installer (Windows).
#   irm https://pinellolab.github.io/CRISPRme/install.ps1 | iex
#
# Installs a clickable "CRISPRme" app (3 big buttons: Start / Update / Stop). The
# first Start downloads the reference + variant data automatically, then opens the
# web interface — so after this one command the user never needs the terminal
# again. Docker Desktop required.
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

# ---- 3. Write the clickable app (3 big buttons) ---------------------------
New-Item -ItemType Directory -Force -Path $AppDir | Out-Null
$AppScript = Join-Path $AppDir 'CRISPRme.ps1'
@'
Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Drawing
$DataDir = Join-Path $env:USERPROFILE 'CRISPRme'
$Compose = Join-Path $DataDir 'docker-compose.yml'
$Image   = 'pinellolab/crisprme:v2.4.0'
$Mount   = "-v `"$DataDir\crisprme-data:/DATA`" -w /DATA $Image"

function Has-Data { Test-Path (Join-Path $DataDir 'crisprme-data\Genomes') }
function Docker-Running { try { docker info *> $null; return $true } catch { return $false } }
function Ensure-Docker {
  if (Docker-Running) { return $true }
  $dd = Join-Path $env:ProgramFiles 'Docker\Docker\Docker Desktop.exe'
  if (Test-Path $dd) { Start-Process $dd } else { [System.Windows.Forms.MessageBox]::Show('Docker Desktop does not seem to be installed. Install it from https://www.docker.com/products/docker-desktop/ then try again.'); return $false }
  for ($i=0; $i -lt 30; $i++) { Start-Sleep 3; if (Docker-Running) { return $true } }
  [System.Windows.Forms.MessageBox]::Show('Docker Desktop is still starting. Give it a few more seconds, then click Start again.'); return $false
}
function Run-InConsole($title, $cmd) {
  Start-Process powershell -ArgumentList @('-NoExit','-Command',
    "Write-Host '== $title =='; $cmd; Write-Host ''; Write-Host 'DONE - you can close this window.'")
}

$form = New-Object System.Windows.Forms.Form
$form.Text = 'CRISPRme+'
$form.Size = New-Object System.Drawing.Size(440,320)
$form.StartPosition = 'CenterScreen'; $form.FormBorderStyle = 'FixedSingle'; $form.MaximizeBox = $false
$status = New-Object System.Windows.Forms.Label
$status.SetBounds(20,12,400,24); $form.Controls.Add($status)
$big = New-Object System.Drawing.Font('Segoe UI',15,[System.Drawing.FontStyle]::Bold)

function Mk($text,$y,$action) {
  $b = New-Object System.Windows.Forms.Button
  $b.SetBounds(30,$y,380,64); $b.Text = $text; $b.Font = $big; $b.Add_Click($action)
  $form.Controls.Add($b)
}
Mk 'Start' 44 {
  if (-not (Ensure-Docker)) { return }
  if (Has-Data) {
    docker compose -f $Compose up -d; Start-Sleep 3; Start-Process 'http://localhost:8080'
  } else {
    $r = [System.Windows.Forms.MessageBox]::Show("First run: CRISPRme will download the reference genome + 1000G/HGDP variant data (~85 GB), once. A window shows progress, then your browser opens automatically.","CRISPRme+ - first-time setup",'OKCancel')
    if ($r -eq 'OK') {
      Run-InConsole 'First-time setup - downloading ~85 GB, then starting CRISPRme' "docker run --rm $Mount crisprme.py download --what all --path /DATA; docker run --rm $Mount crisprme.py download --what index --index-name NRG_3_hg38-dictless+hg38_1000G_HGDP --path /DATA; docker compose -f `"$Compose`" up -d; Start-Sleep 3; Start-Process 'http://localhost:8080'"
    }
  }
}
Mk 'Update' 120 { if (Ensure-Docker) { Run-InConsole 'Updating CRISPRme' "docker pull $Image" } }
Mk 'Stop' 196 { try { docker compose -f $Compose down } catch {}; [System.Windows.Forms.MessageBox]::Show('CRISPRme stopped.') }

if (Has-Data) { $status.Text = 'Ready. Click Start to open CRISPRme.' } else { $status.Text = 'First run: Start will set everything up (~85 GB, once).' }
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
Write-Host "Done! Open 'CRISPRme' from your Desktop or Start Menu and click Start." -ForegroundColor Green
Write-Host "The first Start downloads the data (once) and opens the web app."
