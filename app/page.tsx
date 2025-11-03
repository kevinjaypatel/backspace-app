'use client';

import { useEffect, useState } from 'react';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceLine,
  Area,
  ComposedChart,
  ReferenceArea,
} from 'recharts';

interface SimulationData {
  days: number[];
  tumor: {
    chemo_pulse: number[];
    chemo_metronomic: number[];
    chemo_weekly: number[];
    radiotherapy: number[];
  };
  clonogens_fraction: number;
  surgery_fraction_t0: number;
  config: any;
}

interface LatestSimulationData {
  days: number[];
  tumor: {
    median: number[];
    p25: number[];
    p75: number[];
  };
  clonogens: {
    median: number[];
    p25: number[];
    p75: number[];
  };
  threshold: number;
  treatment_days: number[];
  parameters: {
    N0_tumor: number;
    clonogenic_fraction: number;
    r_mean: number;
    r_std: number;
    K_mean: number;
    K_std: number;
    kill_mean: number;
    kill_std: number;
    n_runs: number;
  };
}

export default function Home() {
  const [data, setData] = useState<SimulationData | null>(null);
  const [latestData, setLatestData] = useState<LatestSimulationData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [selectedTherapy, setSelectedTherapy] = useState<string>('chemo_pulse');
  const [viewMode, setViewMode] = useState<'original' | 'latest'>('latest');

  useEffect(() => {
    // Fetch both datasets
    Promise.all([
      fetch('/api/simulation').then(res => res.ok ? res.json() : null).catch(() => null),
      fetch('/api/simulation-latest').then(res => res.ok ? res.json() : null).catch(() => null)
    ]).then(([originalData, latest]) => {
      if (originalData) setData(originalData);
      if (latest) setLatestData(latest);
      setLoading(false);
    }).catch((err) => {
      setError(err.message);
      setLoading(false);
    });
  }, []);

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="text-xl">Loading simulation data...</div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="text-red-500 text-xl">Error: {error}</div>
      </div>
    );
  }

  const therapyLabels: Record<string, string> = {
    chemo_pulse: 'Pulse Chemo (q3w ×6)',
    chemo_metronomic: 'Metronomic Chemo (daily ×6wk)',
    chemo_weekly: 'Weekly Dose-Dense Chemo (q1w ×8)',
    radiotherapy: 'Radiotherapy (50 Gy / 25 fx)',
  };

  // Prepare chart data for original view
  const chartData = data?.days.map((day, index) => {
    const clonogensFraction = data?.clonogens_fraction || 0.01;
    return {
      day,
      tumorCells: data?.tumor[selectedTherapy as keyof typeof data.tumor]?.[index] || 0,
      clonogens: (data?.tumor[selectedTherapy as keyof typeof data.tumor]?.[index] || 0) * clonogensFraction,
    };
  }) || [];

  // Prepare chart data for latest simulation (with IQR bands)
  // Use actual P25 and P75 values - Recharts handles log scale stacking internally
  const latestChartData = latestData?.days.map((day, index) => ({
    day,
    tumorMedian: latestData.tumor.median[index],
    tumorP25: latestData.tumor.p25[index],
    tumorP75: latestData.tumor.p75[index],
    clonoMedian: latestData.clonogens.median[index],
    clonoP25: latestData.clonogens.p25[index],
    clonoP75: latestData.clonogens.p75[index],
  })) || [];

  return (
    <div className="container mx-auto p-8">
      <h1 className="text-4xl font-bold mb-8 text-center">
        Cancer Simulation Visualization
      </h1>

      {/* View Mode Selector */}
      <div className="mb-6 flex gap-4 items-center">
        <label className="text-sm font-medium">View:</label>
        <button
          onClick={() => setViewMode('latest')}
          className={`px-4 py-2 rounded-md ${
            viewMode === 'latest'
              ? 'bg-blue-500 text-white'
              : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
          }`}
        >
          Latest Simulation (Stochastic)
        </button>
        {data && (
          <button
            onClick={() => setViewMode('original')}
            className={`px-4 py-2 rounded-md ${
              viewMode === 'original'
                ? 'bg-blue-500 text-white'
                : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
            }`}
          >
            Original Therapies
          </button>
        )}
      </div>

      {viewMode === 'latest' && latestData ? (
        <>
          <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
            <h2 className="text-2xl font-semibold mb-4">
              Tumor and Clonogenic Cell Dynamics with Stochastic Treatments
            </h2>
            
            {/* Display the PNG image */}
            <div className="w-full flex justify-center">
              <img
                src="/api/simulation-image"
                alt="Cancer Simulation Graph"
                className="max-w-full h-auto"
                onError={(e) => {
                  console.error('Error loading image:', e);
                  (e.target as HTMLImageElement).src = '';
                }}
              />
            </div>
          </div>

          <div className="bg-gray-50 rounded-lg p-6">
            <h3 className="text-xl font-semibold mb-4">Simulation Information</h3>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <span className="font-medium">Monte Carlo Runs:</span>{' '}
                {latestData.parameters.n_runs}
              </div>
              <div>
                <span className="font-medium">Initial Tumor Cells:</span>{' '}
                {latestData.parameters.N0_tumor.toExponential(2)}
              </div>
              <div>
                <span className="font-medium">Clonogenic Fraction:</span>{' '}
                {(latestData.parameters.clonogenic_fraction * 100).toFixed(3)}%
              </div>
              <div>
                <span className="font-medium">Control Threshold:</span>{' '}
                {latestData.threshold} clonogen
              </div>
              <div>
                <span className="font-medium">Growth Rate (mean):</span>{' '}
                {latestData.parameters.r_mean.toFixed(3)}/day
              </div>
              <div>
                <span className="font-medium">Treatment Days:</span>{' '}
                {latestData.treatment_days.join(', ')}
              </div>
            </div>
          </div>
        </>
      ) : viewMode === 'original' && data ? (
        <>
          <div className="mb-6">
            <label className="block text-sm font-medium mb-2">
              Select Therapy Type:
            </label>
            <select
              value={selectedTherapy}
              onChange={(e) => setSelectedTherapy(e.target.value)}
              className="block w-full max-w-xs px-4 py-2 border border-gray-300 rounded-md shadow-sm focus:ring-blue-500 focus:border-blue-500"
            >
              {Object.keys(data.tumor).map((therapy) => (
                <option key={therapy} value={therapy}>
                  {therapyLabels[therapy] || therapy}
                </option>
              ))}
            </select>
          </div>

          <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
            <h2 className="text-2xl font-semibold mb-4">
              {therapyLabels[selectedTherapy] || selectedTherapy}
            </h2>
            
            <ResponsiveContainer width="100%" height={500}>
              <LineChart data={chartData}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis 
                  dataKey="day" 
                  label={{ value: 'Days since diagnosis', position: 'insideBottom', offset: -5 }}
                />
                <YAxis 
                  scale="log" 
                  domain={['auto', 'auto']}
                  label={{ value: 'Cell count (log scale)', angle: -90, position: 'insideLeft' }}
                />
                <Tooltip 
                  formatter={(value: number, name: string) => [
                    typeof value === 'number' ? value.toExponential(2) : value,
                    name === 'tumorCells' ? 'Tumor Cells' : 'Clonogens',
                  ]}
                />
                <Legend />
                <Line
                  type="monotone"
                  dataKey="tumorCells"
                  stroke="#8884d8"
                  strokeWidth={2}
                  name="Tumor Cells"
                  dot={false}
                />
                <Line
                  type="monotone"
                  dataKey="clonogens"
                  stroke="#82ca9d"
                  strokeWidth={2}
                  strokeDasharray="5 5"
                  name={`Clonogens (${(data.clonogens_fraction * 100).toFixed(1)}%)`}
                  dot={false}
                />
              </LineChart>
            </ResponsiveContainer>
          </div>

          <div className="bg-gray-50 rounded-lg p-6">
            <h3 className="text-xl font-semibold mb-4">Simulation Information</h3>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <span className="font-medium">Clonogens Fraction:</span>{' '}
                {(data.clonogens_fraction * 100).toFixed(1)}%
              </div>
              <div>
                <span className="font-medium">Surgery Fraction (t=0):</span>{' '}
                {(data.surgery_fraction_t0 * 100).toFixed(1)}%
              </div>
              <div>
                <span className="font-medium">Total Days Simulated:</span>{' '}
                {data.days.length}
              </div>
              <div>
                <span className="font-medium">Duration:</span>{' '}
                {data.days[data.days.length - 1]} days
              </div>
            </div>
          </div>
        </>
      ) : (
        <div className="text-center text-gray-500">
          No data available
        </div>
      )}
    </div>
  );
}

