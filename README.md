# Backspace App - Cancer Simulation Visualization

A full-stack application for visualizing cancer progression simulations with a Python backend and Next.js frontend.

## Quick Start for New Contributors

This guide will help you get the project running locally.

### Prerequisites

- **Node.js 18+** and npm (for the Next.js frontend)
- **Python 3.8+** (for the simulation scripts)
- **Git** (to clone the repository)

### Setup Instructions

#### 1. Clone the Repository

```bash
git clone https://github.com/kevinjaypatel/backspace-app.git
cd backspace-app
```

#### 2. Set Up Python Environment

Navigate to the Python directory and set up the virtual environment:

```bash
cd python
python3 -m venv venv
```

Activate the virtual environment:

**On macOS/Linux:**

```bash
source venv/bin/activate
```

**On Windows:**

```bash
venv\Scripts\activate
```

Install Python dependencies:

```bash
pip install -r requirements.txt
```

#### 3. Run Python Simulations (Optional)

To generate simulation data and images:

```bash
# Run the latest simulation
python python-latest.py

# Or run other simulations
python cancer_sim2.py
python cancer_sim.py
```

This will create output files in the `../outputs/` directory.

#### 4. Set Up Next.js Frontend

Go back to the root directory:

```bash
cd ..
```

Install Node.js dependencies:

```bash
npm install
```

#### 5. Run the Development Server

Start the Next.js development server:

```bash
npm run dev
```

The application will be available at **http://localhost:3000**

### Project Structure

```
backspace-app/
├── app/                    # Next.js frontend (App Router)
│   ├── api/               # API routes
│   │   ├── simulation/    # Serves JSON data
│   │   └── simulation-image/  # Serves PNG images
│   ├── page.tsx           # Main page component
│   └── layout.tsx         # Root layout
├── python/                # Python simulation scripts
│   ├── python-latest.py   # Latest simulation script
│   ├── cancer_sim.py      # Cancer progression simulator
│   ├── cancer_sim2.py     # Stage I cancer model
│   └── venv/              # Python virtual environment
├── outputs/               # Simulation outputs
│   ├── python_latest_simulation.png
│   ├── python_latest_simulation.json
│   └── stage1_surgery_repop_outputs.json
└── package.json           # Node.js dependencies
```

### Available Scripts

**Node.js:**

- `npm run dev` - Start development server
- `npm run build` - Build for production
- `npm start` - Start production server
- `npm run lint` - Run ESLint

**Python:**

- `python python-latest.py` - Run latest stochastic simulation
- `python cancer_sim2.py` - Run Stage I cancer model with multiple therapies
- `python cancer_sim.py` - Run cancer progression simulator

### API Endpoints

- `GET /api/simulation` - Returns Stage I surgery repopulation data (JSON)
- `GET /api/simulation-latest` - Returns latest stochastic simulation data (JSON)
- `GET /api/simulation-image` - Returns latest simulation graph (PNG)

### Troubleshooting

**Python Virtual Environment Issues:**

- If activation fails, make sure you're using the correct path separator for your OS
- On Windows PowerShell, you may need: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`

**Node.js Issues:**

- Make sure Node.js version is 18 or higher: `node --version`
- Clear `node_modules` and reinstall: `rm -rf node_modules && npm install`

**Missing Output Files:**

- Run the Python scripts first to generate the JSON and PNG files
- Check that files are in the `outputs/` directory

### Contributing

1. Create a new branch for your feature
2. Make your changes
3. Test locally
4. Commit and push to your branch
5. Open a pull request

For more details, see the individual README files:

- `README_APP.md` - Frontend documentation
- `python/README.md` - Python scripts documentation
