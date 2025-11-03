import { NextResponse } from 'next/server';
import fs from 'fs';
import path from 'path';

export async function GET() {
  try {
    // Path to the PNG image file
    const imagePath = path.join(process.cwd(), 'outputs', 'python_latest_simulation.png');
    
    // Check if file exists
    if (!fs.existsSync(imagePath)) {
      return NextResponse.json(
        { error: 'Simulation image file not found' },
        { status: 404 }
      );
    }

    // Read the image file
    const imageBuffer = fs.readFileSync(imagePath);
    
    // Return the image with proper headers
    return new NextResponse(imageBuffer, {
      headers: {
        'Content-Type': 'image/png',
        'Cache-Control': 'public, max-age=3600',
      },
    });
  } catch (error) {
    console.error('Error reading simulation image:', error);
    return NextResponse.json(
      { error: 'Failed to read simulation image' },
      { status: 500 }
    );
  }
}

