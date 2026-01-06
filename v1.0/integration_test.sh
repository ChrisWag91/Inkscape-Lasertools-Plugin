#!/bin/bash

# --- 1. Configuration & Path Setup ---
# Get the directory where THIS script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set the project root to one level above the script directory
# This allows finding the 'test_Files' folder regardless of absolute location
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TEST_FILES_DIR="$PROJECT_ROOT/test_Files"

# Standard Inkscape command
INKSCAPE_BIN="inkscape"

# List of test SVG files from your specifications
TEST_FILES=(
    "test01_general_function_test_Panther.svg"
    "test02_general_function_test_TestPattern.svg"
    "test03_performance_test_medium.svg"
    #"test04_paths_and_images.svg"
    #"test05_scaling_and_text_Opensource_logo.svg"
    #"test06_line_glitches_Open-source-hardware-logo.svg"
    #"test07_line_glitches_Open-source-hardware-logo_reimport.svg"
    #"test09_line_Glitches.svg"
    #"test11_firstLine glitch.svg"
)

# --- 2. Execute Deployment ---
echo "Step 1: Running deployment..."
if [ -f "$SCRIPT_DIR/deploy.sh" ]; then
    # Execute deploy.sh within its own directory context
    (cd "$SCRIPT_DIR" && bash ./deploy.sh)
else
    echo "Error: deploy.sh not found in $SCRIPT_DIR"
    exit 1
fi

# Check if deployment was successful before continuing
if [ $? -ne 0 ]; then
    echo "Deployment failed. Aborting test cycle."
    exit 1
fi

# --- 3. Run SVG Test Cycle ---
echo "----------------------------------------------------"
echo "Step 2: Starting SVG Test Cycle..."
echo "Project Root: $PROJECT_ROOT"
echo "Test Directory: $TEST_FILES_DIR"
echo "----------------------------------------------------"

for FILE in "${TEST_FILES[@]}"; do
    FULL_PATH="$TEST_FILES_DIR/$FILE"
    
    if [ -f "$FULL_PATH" ]; then
        echo "Launching: $FILE"
        
        if [[ "$OSTYPE" == "linux-gnu"* ]]; then
            # Linux Launch: Standard background process
            
            $INKSCAPE_BIN "$FULL_PATH" -g --actions="org.lasertools.lasertools.noprefs" &
                     
        elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
            # Windows Launch: Use 'start' to detach process and 'cygpath' for Windows-native paths
            WIN_PATH=$(cygpath -w "$FULL_PATH")
            start inkscape "$WIN_PATH" -g --actions=org.lasertools.lasertools.noprefs
        else
            # Fallback for other systems
            $INKSCAPE_BIN "$FULL_PATH" -g --actions="org.lasertools.lasertools.noprefs" &
        fi
        
        # Brief pause to maintain system stability during batch launch
        sleep 60
    else
        echo "Skipping: $FILE (File not found at $FULL_PATH)"
    fi
done

echo "----------------------------------------------------"
echo "All tests initiated. Terminal is now free for logs."