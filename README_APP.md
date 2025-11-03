# Backspace App - Cancer Simulation Visualization

A Next.js 14 application that visualizes cancer progression simulation data with interactive charts.

## Features

- **API Endpoint**: Fetches simulation JSON data from `outputs/stage1_surgery_repop_outputs.json`
- **Interactive Charts**: Visualize tumor cell growth and clonogen counts over time using Recharts
- **Multiple Therapy Types**: View different treatment regimens:
  - Pulse Chemo (q3w ×6)
  - Metronomic Chemo (daily ×6wk)
  - Weekly Dose-Dense Chemo (q1w ×8)
  - Radiotherapy (50 Gy / 25 fx)
- **Log Scale Visualization**: Charts use logarithmic scale for better visualization of cell counts

## Getting Started

### Prerequisites

- Node.js 18+
- npm or yarn

### Installation

1. Install dependencies:

```bash
npm install
```

2. Make sure you have the simulation data:
   - Run the Python simulation scripts to generate `outputs/stage1_surgery_repop_outputs.json`
   - Or ensure the file exists in the `outputs/` directory

### Running the Development Server

```bash
npm run dev
```

Open [http://localhost:3000](http://localhost:3000) in your browser to see the visualization.

### Building for Production

```bash
npm run build
npm start
```

## Project Structure

```
├── app/
│   ├── api/
│   │   └── simulation/
│   │       └── route.ts          # API endpoint to serve JSON data
│   ├── layout.tsx                 # Root layout
│   ├── page.tsx                   # Main page with chart visualization
│   └── globals.css                # Global styles
├── outputs/
│   └── stage1_surgery_repop_outputs.json  # Simulation data
└── package.json                   # Dependencies
```

## API Endpoint

### GET `/api/simulation`

Returns the simulation data from `outputs/stage1_surgery_repop_outputs.json`

**Response:**

```json
{
  "days": [0, 1, 2, ...],
  "tumor": {
    "chemo_pulse": [...],
    "chemo_metronomic": [...],
    "chemo_weekly": [...],
    "radiotherapy": [...]
  },
  "clonogens_fraction": 0.01,
  "surgery_fraction_t0": 0.95,
  "config": {...}
}
```

## Technologies Used

- **Next.js 14** - React framework with App Router
- **TypeScript** - Type safety
- **Tailwind CSS** - Styling
- **Recharts** - Chart visualization library
