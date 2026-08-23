#!/bin/bash
# CRISPRme+ one-line installer (macOS / Linux).
#   curl -fsSL https://pinellolab.github.io/CRISPRme/install.sh | bash
#
# Installs a clickable "CRISPRme" app (3 buttons: Start / Update / Stop) that
# manages the Docker container and opens the web interface. On the FIRST Start it
# downloads the reference + variant data automatically, then opens the browser —
# so after this one command the user never needs the terminal again. No build
# toolchain: the macOS app is generated with osacompile (built in). Docker Desktop
# is required (guided if missing).
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

# ---- 1b. Pre-pull the engine image (best-effort, visible progress) --------
# Front-load the ~800 MB image now — with progress here in the terminal — so the
# app's first Start (which probes a bind mount using this image) doesn't appear to
# hang while it pulls. Skipped silently if the daemon isn't up yet; the app pulls
# it on demand in that case.
if docker info >/dev/null 2>&1; then
  say "Fetching the CRISPRme engine image (~800 MB, one time)…"
  docker pull "$IMAGE" || true
fi

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
-- CRISPRme+ launcher. Three actions: Start / Update / Stop.
property img : "pinellolab/crisprme:v2.4.0"
property dataDir : ""

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

-- a bind mount actually works (a freshly-started Docker can report 'info' OK while
-- its file-sharing subsystem is still initializing, which makes the download's
-- makedirs fail with 'No such file or directory')
on mountReady(dk)
  return (do shell script (quoted form of dk) & " run --rm -v " & quoted form of (dataDir & "/crisprme-data") & ":/DATA " & img & " bash -c 'touch /DATA/.probe && rm -f /DATA/.probe' >/dev/null 2>&1 && echo yes || echo no") is "yes"
end mountReady

-- start Docker Desktop for the user and wait for the daemon AND the bind mount
on ensureDocker(dk)
  if not dockerRunning(dk) then
    try
      do shell script "open -a Docker"
    on error
      display dialog "Docker Desktop does not seem to be installed. Install it from" & return & "https://www.docker.com/products/docker-desktop/" & return & "then open CRISPRme again." buttons {"OK"} default button "OK" with icon caution
      return false
    end try
    set ok to false
    repeat 30 times
      delay 3
      if dockerRunning(dk) then
        set ok to true
        exit repeat
      end if
    end repeat
    if not ok then
      display dialog "Docker Desktop is still starting up. Give it a few more seconds, then open CRISPRme again." buttons {"OK"} default button "OK" with icon caution
      return false
    end if
  end if
  -- daemon is up; now wait for the bind-mount / file-sharing subsystem
  repeat 20 times
    if mountReady(dk) then return true
    delay 2
  end repeat
  return true
end ensureDocker

on hasData()
  return (do shell script "ls " & quoted form of (dataDir & "/crisprme-data/Genomes") & " >/dev/null 2>&1 && echo yes || echo no") is "yes"
end hasData

on runInTerminal(cmd, title)
  tell application "Terminal"
    activate
    do script "clear; echo '== " & title & " =='; " & cmd
  end tell
end runInTerminal

on run
  set dataDir to (POSIX path of (path to home folder)) & "CRISPRme"
  set dk to dockerPath()
  set compose to (quoted form of dk) & " compose -f " & quoted form of (dataDir & "/docker-compose.yml") & " "
  set dmount to "-v " & quoted form of (dataDir & "/crisprme-data") & ":/DATA -w /DATA " & img
  set dlRef to (quoted form of dk) & " run --rm " & dmount & " crisprme.py download --what all --path /DATA"
  set dlVar to (quoted form of dk) & " run --rm " & dmount & " crisprme.py download --what index --index-name NRG_3_hg38-dictless+hg38_1000G_HGDP --path /DATA"

  set action to button returned of (display dialog "CRISPRme+  —  CRISPR off-target analysis" & return & return & "Start   —  open the web app (sets everything up the first time)" & return & "Update  —  get the latest CRISPRme" & return & "Stop    —  shut CRISPRme down" buttons {"Stop", "Update", "Start"} default button "Start" with title "CRISPRme+")

  if action is "Start" then
    if not ensureDocker(dk) then return
    if hasData() then
      -- data already downloaded: just start + open the browser
      do shell script compose & "up -d"
      delay 3
      do shell script "open http://localhost:8080"
      display dialog "CRISPRme is starting — your browser is opening http://localhost:8080." buttons {"OK"} default button "OK" giving up after 6
    else
      -- first run: download reference + variant data, then start + open browser
      set r to button returned of (display dialog "First run: CRISPRme will download the reference genome + the 1000G/HGDP variant data (~85 GB) so everything is ready to go. You only do this once. A Terminal window shows the progress, then your browser opens automatically." buttons {"Not now", "Download & Start"} default button "Download & Start" with title "CRISPRme+ — first-time setup")
      if r is "Download & Start" then
        runInTerminal(dlRef & " && " & dlVar & " && " & compose & "up -d && sleep 3 && open http://localhost:8080 && echo && echo 'DONE — CRISPRme is running and your browser is open. You can close this window.'", "First-time setup — downloading ~85 GB, then starting CRISPRme")
      end if
    end if

  else if action is "Update" then
    if not ensureDocker(dk) then return
    runInTerminal((quoted form of dk) & " pull " & img, "Updating CRISPRme")

  else if action is "Stop" then
    try
      do shell script compose & "down"
    end try
    display dialog "CRISPRme stopped." buttons {"OK"} default button "OK" giving up after 3
  end if
end run
OSA
  osacompile -o "$APP" "$SRC" 2>/dev/null
  rm -f "$SRC"
  # zero-warning open: strip the quarantine flag (installed by a script the user ran)
  xattr -dr com.apple.quarantine "$APP" 2>/dev/null || true
  say "${GRN}Installed:${NC} ${APP}"
  say ""
  say "${BLUE}Done!${NC} Open ${GRN}${APP_NAME}${NC} from your Applications (Spotlight: type 'CRISPRme') and click ${GRN}Start${NC}."
  say "The first Start downloads the data (once) and opens the web app; after that Start is instant."
  open "$APPS" 2>/dev/null || true
else
  # Linux: minimal — run the web interface via compose directly
  say "${YEL}Linux clickable-app support is minimal in this MVP.${NC}"
  say "Use: cd ${DATA_DIR} && docker compose up   (then open http://localhost:8080)"
fi
