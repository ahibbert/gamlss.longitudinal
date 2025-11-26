#!/bin/bash

# Auto-startup script for R plotting in VS Code
# This script automatically opens the plot viewer

echo "🚀 Starting R development environment..."

# Check if VS Code is running
if command -v code >/dev/null 2>&1; then
    echo "📊 Opening R Plot Viewer..."
    # Open the plot viewer URL in VS Code Simple Browser
    code --command "simpleBrowser.show" "http://127.0.0.1:9001/live?token=vscode-r-plots"
    echo "✅ Plot viewer opened at: http://127.0.0.1:9001/live?token=vscode-r-plots"
else
    echo "⚠️  VS Code not found in PATH"
    echo "📋 Manual setup: Press Ctrl+Alt+P or open http://127.0.0.1:9001/live?token=vscode-r-plots"
fi

echo "💡 Remember to start R and run: source('start_httpgd.R')"
