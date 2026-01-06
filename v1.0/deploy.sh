#!/bin/bash

# --- Configuration ---
# Collect all plugin files
PLUGIN_FILES_INX=(*.inx)
PLUGIN_FILES_PY=(*.py)

# Initialize start flag
START_INKSCAPE=false

# Check for arguments
while getopts "n" opt; do
  case $opt in
    n)
      START_INKSCAPE=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# --- Path Detection & Privilege Setup ---
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # System-wide path for Linux
    TARGET_DIR="/usr/share/inkscape/extensions"
    COPY_CMD="sudo cp"
    # Linux launch command
    LAUNCH_CMD="inkscape -g --actions='file-new'"
    echo "Platform: Linux detected."

elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
    # System-wide path for Windows
    TARGET_DIR="$PROGRAMFILES\Inkscape\share\inkscape\extensions"
    COPY_CMD="cp" 
    # Windows launch command (using start to detach process)
    LAUNCH_CMD="start inkscape -g --actions=file-new"
    echo "Platform: Windows detected."
    
    if [ ! -d "$TARGET_DIR" ]; then
        echo "Error: Target directory $TARGET_DIR not found."
        exit 1
    fi
else
    echo "Unknown Operating System. Manual deployment required."
    exit 1
fi

# --- Kill Running Instances ---
# Inkscape must be closed to refresh the extensions cache (extensions.dat)
echo "Closing running Inkscape instances..."
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    pkill -x inkscape 2>/dev/null
else
    taskkill //F //IM inkscape.exe //T 2>/dev/null
fi

# --- Deployment ---
echo "Deploying files to: $TARGET_DIR"

# Copying files with error checking
$COPY_CMD "${PLUGIN_FILES_INX[@]}" "$TARGET_DIR/" && \
$COPY_CMD "${PLUGIN_FILES_PY[@]}" "$TARGET_DIR/"

if [ $? -eq 0 ]; then
    echo "Deployment successful."
    echo "----------------------------------------------------"
    
    # Execute launch only if -n flag was provided
    if [ "$START_INKSCAPE" = true ]; then
        echo "Launching Inkscape with a new document..."
        if [[ "$OSTYPE" == "linux-gnu"* ]]; then
            # Run in background to keep the terminal interactive
            eval $LAUNCH_CMD &
        else
            eval $LAUNCH_CMD
        fi
    else
        echo "Skipping Inkscape launch (no -n argument provided)."
    fi
else
    echo "Deployment failed! Please check if you have Admin/Sudo permissions."
    exit 1
fi