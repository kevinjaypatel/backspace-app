import { NextResponse } from 'next/server';
import fs from 'fs';
import path from 'path';

export async function GET() {
  try {
    // Path to the JSON output file
    const filePath = path.join(process.cwd(), 'outputs', 'stage1_surgery_repop_outputs.json');
    
    // Check if file exists
    if (!fs.existsSync(filePath)) {
      return NextResponse.json(
        { error: 'Simulation data file not found' },
        { status: 404 }
      );
    }

    // Read and parse the JSON file
    const fileContents = fs.readFileSync(filePath, 'utf-8');
    const data = JSON.parse(fileContents);

    return NextResponse.json(data);
  } catch (error) {
    console.error('Error reading simulation data:', error);
    return NextResponse.json(
      { error: 'Failed to read simulation data' },
      { status: 500 }
    );
  }
}

