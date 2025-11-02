# Python venv setup for the cancer progression script

This workspace contains a simple tumor progression simulation script (`import json.py`).

Quick steps (PowerShell, run from the workspace root `c:\Users\reryh\Music\New folder`):

1) Create a virtual environment

    python -m venv venv

2) (Optional) If PowerShell blocks activation due to execution policy, run:

    Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force

3) Activate the venv (PowerShell)

    .\venv\Scripts\Activate.ps1

   Or for cmd.exe:

    .\venv\Scripts\activate.bat

4) Upgrade pip and install dependencies

    python -m pip install --upgrade pip
    pip install -r requirements.txt

5) Run the script

    python ".\import json.py"

   Note: The script filename contains a space and starts with the words `import json`, which can be confusing. I recommend renaming the file to something like `cancer_sim.py` and running:

    python cancer_sim.py

6) Verify output

The script should write `cancer_progression.json` in the same folder. Check that file to confirm the simulation ran.

Why rename? The current filename `import json.py` is valid but awkward: it contains a space and looks like an import statement. Renaming avoids quoting and confusion and is safer for imports or reuse.

If you'd like, I can:
- Create the venv and install packages for you (I can run the commands here if you allow me to run terminal commands),
- Or rename the file to a safer name and update the README and run a quick test.

Tell me which of those you'd like me to do next.