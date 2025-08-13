@echo off
cd /d "%~dp0"

python --version >nul 2>&1
if errorlevel 1 (
    echo Python is not installed or not in PATH. Please install Python >= 3.10
    pause
    exit /b
)

for /f "tokens=2 delims= " %%v in ('python --version') do set PYVER=%%v
for /f "tokens=1,2 delims=." %%a in ("%PYVER%") do (
    set PYMAJOR=%%a
    set PYMINOR=%%b
)
if %PYMAJOR% LSS 3 (
    echo Python version too low. Please install >= 3.10
    pause
    exit /b
)
if %PYMAJOR%==3 if %PYMINOR% LSS 10 (
    echo Python version too low. Please install >= 3.10
    pause
    exit /b
)

echo Python detected: %PYVER%

set PACKAGES=flask aptamerutils

for %%p in (%PACKAGES%) do (
    python -c "import %%p" 2>nul
    if errorlevel 1 (
        echo Installing Python package: %%p
        python -m pip install --user %%p
        if errorlevel 1 (
            echo Failed to install %%p. Please install it manually
            pause
            exit /b
        )
    ) else (
        echo Python package already installed: %%p
    )
)

start "" python app.py

:wait_loop
powershell -command "try { $c = New-Object Net.Sockets.TcpClient('127.0.0.1', 5000); $c.Close(); exit 0 } catch { exit 1 }"
if errorlevel 1 (
    timeout /t 1 >nul
    goto wait_loop
)

start "" http://127.0.0.1:5000