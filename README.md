# Curated Projects
The repository contains special scripts for affiliated data analysis.

# Deploy `curated_projects` build

## Prepare environment

```
echo Export variables
export TOOL_NAME="curated_projects"

export TOOL_DIR="/opt/${TOOL_NAME}/"
export TOOL_SCRIPT="${TOOL_DIR}${TOOL_NAME}.sh"
export TOOL_SERVICE="/etc/systemd/system/${TOOL_NAME}.service"
export UN="$(whoami)"

echo Create directories
sudo rm \
    --force \
    --recursive \
    --verbose \
    "${TOOL_DIR}"

sudo mkdir \
    --parent \
    --mode 0700 \
    --verbose \
    "${TOOL_DIR}"

sudo chown \
    --recursive \
    --verbose \
    "$(id --user "${UN}")" \
    "${TOOL_DIR}"
```

## Configure and start tool

```
echo Create tool script
curl -fsSL \
    "https://raw.githubusercontent.com/ivasilyev/curated_projects/master/run-ci.sh" \
    -o "${TOOL_SCRIPT}"

sudo chmod a+x "${TOOL_SCRIPT}"
# nano "${TOOL_SCRIPT}"



echo Create system service
cat <<EOF | sudo tee "${TOOL_SERVICE}"
[Unit]
Description=${TOOL_NAME}
Documentation=https://google.com
Wants=network-online.target
After=network-online.target

[Service]
Type=simple
User=${UN}
ExecReload=/usr/bin/env kill -s SIGTERM \$MAINPID
ExecStart=/usr/bin/env bash "${TOOL_SCRIPT}"
SyslogIdentifier=${TOOL_NAME}
Restart=always
RestartSec=5

[Install]
WantedBy=multi-user.target
EOF

# nano "${TOOL_SERVICE}"



echo Activate ${TOOL_NAME} service
sudo systemctl daemon-reload
sudo systemctl enable "${TOOL_NAME}.service"
sudo systemctl restart "${TOOL_NAME}.service"
sudo systemctl status "${TOOL_NAME}.service"
```
