#!/bin/bash
# CRISPRme+ one-line installer (macOS / Linux).
#   curl -fsSL https://pinellolab.github.io/CRISPRme/install.sh | bash
#
# Installs a clickable "CRISPRme" app that manages the Docker container, downloads
# reference data, and opens the web interface — so after this one command the user
# never needs the terminal again. No build toolchain: the macOS app is generated
# with osacompile (built in). Docker Desktop is required (guided if missing).
set -euo pipefail

IMAGE="pinellolab/crisprme:v2.4.0"
APP_NAME="CRISPRme"
DATA_DIR="${HOME}/CRISPRme"          # compose file + ./crisprme-data live here
BLUE=$'\033[1;34m'; GRN=$'\033[1;32m'; YEL=$'\033[1;33m'; RED=$'\033[1;31m'; NC=$'\033[0m'
say() { printf "%s\n" "$*"; }

say "${BLUE}=== CRISPRme+ installer ===${NC}"

OS="$(uname -s)"
if [ "$OS" != "Darwin" ] && [ "$OS" != "Linux" ]; then
  say "${RED}Unsupported OS ($OS). On Windows use install.ps1.${NC}"; exit 1
fi

# ---- 1. Docker check -------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then
  say "${YEL}Docker is not installed.${NC} CRISPRme runs inside Docker."
  if [ "$OS" = "Darwin" ] && command -v brew >/dev/null 2>&1; then
    say "Installing Docker Desktop via Homebrew..."
    brew install --cask docker || true
  fi
  if ! command -v docker >/dev/null 2>&1; then
    say "Please install Docker Desktop, then re-run this installer:"
    say "  ${BLUE}https://www.docker.com/products/docker-desktop/${NC}"
    [ "$OS" = "Darwin" ] && open "https://www.docker.com/products/docker-desktop/" 2>/dev/null || true
    exit 1
  fi
fi
say "${GRN}Docker found.${NC}"

# ---- 2. Working dir + compose file ----------------------------------------
mkdir -p "${DATA_DIR}/crisprme-data"
cat > "${DATA_DIR}/docker-compose.yml" <<YML
# Managed by the CRISPRme app. Data persists in ./crisprme-data.
services:
  crisprme:
    image: ${IMAGE}
    working_dir: /DATA
    volumes:
      - ./crisprme-data:/DATA
    ports:
      - "8080:8080"
    command: crisprme.py web-interface
    stdin_open: true
    tty: true
    restart: "no"
YML
say "${GRN}Data folder ready:${NC} ${DATA_DIR}"

# ---- 3. Build the clickable app -------------------------------------------
if [ "$OS" = "Darwin" ]; then
  APPS="${HOME}/Applications"; mkdir -p "$APPS"
  APP="${APPS}/${APP_NAME}.app"
  rm -rf "$APP"
  SRC="$(mktemp -t crisprme_app).applescript"
  cat > "$SRC" <<'OSA'
-- CRISPRme+ launcher (clickable app window). Manages Docker + browser.
property dataDir : ""
property dockerBin : "/usr/local/bin/docker"

on dockerPath()
  set candidates to {"/usr/local/bin/docker", "/opt/homebrew/bin/docker", "/Applications/Docker.app/Contents/Resources/bin/docker"}
  repeat with c in candidates
    if (do shell script "test -x " & quoted form of (c as text) & " && echo yes || echo no") is "yes" then return (c as text)
  end repeat
  return "docker"
end dockerPath

on dockerRunning(dk)
  return (do shell script (quoted form of dk) & " info >/dev/null 2>&1 && echo yes || echo no") is "yes"
end dockerRunning

-- start Docker Desktop for the user and wait for the daemon (up to ~90s)
on ensureDocker(dk)
  if dockerRunning(dk) then return true
  try
    do shell script "open -a Docker"
  on error
    display dialog "Docker Desktop does not seem to be installed. Install it from" & return & "https://www.docker.com/products/docker-desktop/" & return & "then open CRISPRme again." buttons {"OK"} default button "OK" with icon caution
    return false
  end try
  repeat 30 times
    delay 3
    if dockerRunning(dk) then return true
  end repeat
  display dialog "Docker Desktop is still starting up. Give it a few more seconds, then click 'Open CRISPRme' again." buttons {"OK"} default button "OK" with icon caution
  return false
end ensureDocker

on hasData()
  return (do shell script "ls " & quoted form of (dataDir & "/crisprme-data/Genomes") & " >/dev/null 2>&1 && echo yes || echo no") is "yes"
end hasData

on run
  set dataDir to (POSIX path of (path to home folder)) & "CRISPRme"
  set dk to dockerPath()
  repeat
    set runningNote to ""
    try
      if (do shell script dk & " compose -f " & quoted form of (dataDir & "/docker-compose.yml") & " ps --status running -q 2>/dev/null | head -1 | wc -l | tr -d ' '") is not "0" then set runningNote to " (running)"
    end try
    set dataNote to "reference data: MISSING — download it first"
    if hasData() then set dataNote to "reference data: ready"
    set choices to {"Open CRISPRme (start + browser)", "Download reference data", "Download variant index (advanced, 64 GB RAM)", "Update CRISPRme", "Stop CRISPRme", "Quit"}
    set pick to (choose from list choices with title "CRISPRme+" with prompt ("CRISPRme+" & runningNote & return & dataNote & return & return & "Choose an action:") default items {"Open CRISPRme (start + browser)"} OK button name "Go" cancel button name "Close")
    if pick is false then return
    set action to item 1 of pick

    if action starts with "Open CRISPRme" then
      if ensureDocker(dk) then
        if not hasData() then
          display dialog "No reference data yet. Choose 'Download reference data' first (one time)." buttons {"OK"} default button "OK" with icon caution
        else
          do shell script dk & " compose -f " & quoted form of (dataDir & "/docker-compose.yml") & " up -d"
          delay 3
          do shell script "open http://localhost:8080"
          display dialog "CRISPRme is starting — your browser is opening http://localhost:8080. First load can take a few seconds." buttons {"OK"} default button "OK" giving up after 6
        end if
      end if

    else if action starts with "Download reference data" then
      runInTerminal(dk & " run --rm -v " & quoted form of (dataDir & "/crisprme-data") & ":/DATA -w /DATA pinellolab/crisprme:v2.4.0 crisprme.py download --what all --path /DATA", "Downloading reference data (~25 GB, one time)")

    else if action starts with "Download variant index" then
      runInTerminal(dk & " run --rm -v " & quoted form of (dataDir & "/crisprme-data") & ":/DATA -w /DATA pinellolab/crisprme:v2.4.0 crisprme.py download --what index --index-name NRG_3_hg38-dictless+hg38_1000G_HGDP --path /DATA", "Downloading variant index (needs 64 GB RAM in Docker Desktop)")

    else if action starts with "Update CRISPRme" then
      runInTerminal(dk & " pull pinellolab/crisprme:v2.4.0", "Updating CRISPRme engine")

    else if action starts with "Stop CRISPRme" then
      try
        do shell script dk & " compose -f " & quoted form of (dataDir & "/docker-compose.yml") & " down"
      end try
      display dialog "CRISPRme stopped." buttons {"OK"} default button "OK" giving up after 3

    else if action is "Quit" then
      return
    end if
  end repeat
end run

on runInTerminal(cmd, title)
  tell application "Terminal"
    activate
    do script "clear; echo '== " & title & " =='; echo 'You can close this window when it finishes.'; echo; " & cmd & "; echo; echo 'DONE — you can close this window.'"
  end tell
end runInTerminal
OSA
  osacompile -o "$APP" "$SRC" 2>/dev/null
  rm -f "$SRC"
  # zero-warning open: strip the quarantine flag (installed by a script the user ran)
  xattr -dr com.apple.quarantine "$APP" 2>/dev/null || true
  say "${GRN}Installed:${NC} ${APP}"
  say ""
  say "${BLUE}Done!${NC} Open ${GRN}${APP_NAME}${NC} from your Applications (Spotlight: type 'CRISPRme')."
  say "First time: choose ${GRN}Download reference data${NC}, then ${GRN}Open CRISPRme${NC}."
  open "$APPS" 2>/dev/null || true
else
  # Linux: a .desktop launcher that runs the same actions via zenity if present
  say "${YEL}Linux clickable-app support is minimal in this MVP.${NC}"
  say "Use: cd ${DATA_DIR} && docker compose up   (then open http://localhost:8080)"
fi
