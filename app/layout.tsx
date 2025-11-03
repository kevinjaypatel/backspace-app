import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Backspace - Cancer Simulation Visualization",
  description: "Visualize cancer progression simulation data",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}

