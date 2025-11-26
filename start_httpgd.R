# Quick httpgd startup for VS Code with persistent token
# Run this once at the start of your R session

library(httpgd)

# Start httpgd with a fixed token and port for consistency
cat("Starting httpgd server with persistent settings...\n")

# Use fixed settings for consistent access
hgd(host = "127.0.0.1",
    port = 9001,
    token = "vscode-r-plots",
    silent = FALSE)

# Wait for server to start
Sys.sleep(1)

# Fixed URL that will work every session
fixed_url <- "http://127.0.0.1:9001/live?token=vscode-r-plots"

cat("🔗 httpgd server is running at FIXED URL:\n")
cat("�", fixed_url, "\n\n")

cat("✅ PERSISTENT SETUP COMPLETE!\n")
cat("💡 This URL will work every time you start R:\n")
cat("   ", fixed_url, "\n\n")

cat("📋 To open plot viewer (ONE TIME SETUP):\n")
cat("   1. Press Ctrl+Shift+P in VS Code\n")
cat("   2. Type 'Simple Browser' and select 'Simple Browser: Show'\n")
cat("   3. Paste this URL: ", fixed_url, "\n")
cat("   4. BOOKMARK this tab - it will work forever!\n\n")

# Set httpgd as default device
options(device = function() hgd(host = "127.0.0.1", port = 9001, token = "vscode-r-plots"))

cat("💾 Pro tip: Bookmark this URL in VS Code and you'll never need to set it up again!\n\n")

# Test plot to verify
cat("🧪 Creating test plot...\n")
plot(1:10, main = "Persistent httpgd Test - Same URL Every Time!")
cat("📊 Test plot created! Check your viewer at the fixed URL above.\n")
