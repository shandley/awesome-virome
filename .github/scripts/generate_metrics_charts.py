#!/usr/bin/env python3
"""
Script to generate metrics visualization charts.
This creates SVG charts to visualize the metrics history.
"""

import os
import sys
import json
import datetime
from pathlib import Path

def ensure_output_dir():
    """Ensure the charts output directory exists."""
    charts_dir = Path("metrics_history/charts")
    charts_dir.mkdir(parents=True, exist_ok=True)
    return charts_dir

def load_summary():
    """Load the metrics summary file."""
    summary_file = "metrics_history/metrics_summary.json"
    if not os.path.exists(summary_file):
        print(f"Error: Summary file {summary_file} not found")
        return None
    
    with open(summary_file, 'r') as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as e:
            print(f"Error parsing summary file: {e}")
            return None

def generate_performance_chart(summary, charts_dir):
    """Generate a performance trends chart."""
    if not summary or "performance_trends" not in summary:
        return
    
    data_load_times = summary["performance_trends"].get("data_load_time", [])
    readme_parse_times = summary["performance_trends"].get("readme_parse_time", [])
    
    # Create simple placeholder chart if there's no data yet
    if not data_load_times and not readme_parse_times:
        # Generate a simple empty chart with message
        width, height = 800, 400
        svg = [
            f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
            f'  <rect width="{width}" height="{height}" fill="white"/>',
            f'  <text x="{width/2}" y="20" text-anchor="middle" font-family="Arial" font-size="16">Performance Trends</text>',
            f'  <text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="Arial" font-size="14">No performance data available yet</text>',
            f'  <text x="{width/2}" y="{height/2 + 30}" text-anchor="middle" font-family="Arial" font-size="12">Data will appear after workflow runs</text>',
            '</svg>'
        ]
        
        chart_path = os.path.join(charts_dir, "performance_trend.svg")
        with open(chart_path, 'w') as f:
            f.write('\n'.join(svg))
        
        print(f"Empty performance chart generated: {chart_path}")
        return
    
    # Limit to last 20 entries for readability
    data_load_times = data_load_times[-20:]
    readme_parse_times = readme_parse_times[-20:]
    
    # Create a simple SVG line chart
    width, height = 800, 400
    padding = 50
    chart_width = width - 2 * padding
    chart_height = height - 2 * padding
    
    # Determine y-axis scale based on max value
    data_max = max((point["value"] for point in data_load_times), default=0.001)
    readme_max = max((point["value"] for point in readme_parse_times), default=0.001)
    y_max = max(data_max, readme_max, 0.001) * 1.1  # Add 10% margin, ensure non-zero
    
    # Generate SVG content
    svg = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
        f'  <rect width="{width}" height="{height}" fill="white"/>',
        f'  <text x="{width/2}" y="20" text-anchor="middle" font-family="Arial" font-size="16">Performance Trends</text>',
        f'  <text x="{width/2}" y="{height-10}" text-anchor="middle" font-family="Arial" font-size="12">Last {len(data_load_times)} Runs</text>',
        f'  <text x="10" y="{height/2}" text-anchor="middle" font-family="Arial" font-size="12" transform="rotate(-90 10,{height/2})">Time (seconds)</text>'
    ]
    
    # Draw axes
    svg.extend([
        f'  <line x1="{padding}" y1="{height-padding}" x2="{width-padding}" y2="{height-padding}" stroke="black" stroke-width="1"/>',
        f'  <line x1="{padding}" y1="{padding}" x2="{padding}" y2="{height-padding}" stroke="black" stroke-width="1"/>'
    ])
    
    # Draw data_load_time line
    data_points = []
    x_step = chart_width / (len(data_load_times) - 1) if len(data_load_times) > 1 else 0
    
    for i, point in enumerate(data_load_times):
        x = padding + i * x_step
        y = height - padding - (point["value"] / y_max * chart_height)
        data_points.append(f"{x},{y}")
    
    if data_points:
        svg.append(f'  <polyline points="{" ".join(data_points)}" fill="none" stroke="blue" stroke-width="2"/>')
    
    # Draw readme_parse_time line
    readme_points = []
    for i, point in enumerate(readme_parse_times):
        x = padding + i * x_step
        y = height - padding - (point["value"] / y_max * chart_height)
        readme_points.append(f"{x},{y}")
    
    if readme_points:
        svg.append(f'  <polyline points="{" ".join(readme_points)}" fill="none" stroke="red" stroke-width="2"/>')
    
    # Add legend
    svg.extend([
        f'  <rect x="{width-200}" y="40" width="12" height="12" fill="blue"/>',
        f'  <text x="{width-180}" y="50" font-family="Arial" font-size="12">Data Loading Time</text>',
        f'  <rect x="{width-200}" y="60" width="12" height="12" fill="red"/>',
        f'  <text x="{width-180}" y="70" font-family="Arial" font-size="12">README Parsing Time</text>'
    ])
    
    # Close SVG
    svg.append('</svg>')
    
    # Write to file
    chart_path = os.path.join(charts_dir, "performance_trend.svg")
    with open(chart_path, 'w') as f:
        f.write('\n'.join(svg))
    
    print(f"Performance chart generated: {chart_path}")

def generate_validation_chart(summary, charts_dir):
    """Generate a validation success rate chart."""
    if not summary or "validation_stats" not in summary:
        return
    
    total_runs = summary["validation_stats"].get("total_runs", 0)
    successful_runs = summary["validation_stats"].get("successful_runs", 0)
    
    if total_runs == 0:
        # Generate a simple empty chart with message
        width, height = 400, 400
        svg = [
            f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
            f'  <rect width="{width}" height="{height}" fill="white"/>',
            f'  <text x="{width/2}" y="20" text-anchor="middle" font-family="Arial" font-size="16">Validation Success Rate</text>',
            f'  <text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="Arial" font-size="14">No validation data available yet</text>',
            f'  <text x="{width/2}" y="{height/2 + 30}" text-anchor="middle" font-family="Arial" font-size="12">Data will appear after workflow runs</text>',
            '</svg>'
        ]
        
        chart_path = os.path.join(charts_dir, "validation_success.svg")
        with open(chart_path, 'w') as f:
            f.write('\n'.join(svg))
        
        print(f"Empty validation chart generated: {chart_path}")
        return
    
    success_rate = (successful_runs / total_runs) * 100
    failure_rate = 100 - success_rate
    
    # Create a simple SVG pie chart
    width, height = 400, 400
    center_x, center_y = width / 2, height / 2
    radius = min(center_x, center_y) - 50
    
    # Calculate slice angles for pie chart
    success_angle = success_rate * 3.6  # 360 degrees = 100%
    
    # Generate SVG content
    svg = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
        f'  <rect width="{width}" height="{height}" fill="white"/>',
        f'  <text x="{width/2}" y="20" text-anchor="middle" font-family="Arial" font-size="16">Validation Success Rate</text>'
    ]
    
    # Draw success slice (up to 100%)
    if success_rate > 0:
        if success_rate < 100:
            end_x = center_x + radius * math.cos(math.radians(success_angle - 90))
            end_y = center_y + radius * math.sin(math.radians(success_angle - 90))
            large_arc_flag = 1 if success_angle > 180 else 0
            
            svg.append(f'  <path d="M {center_x} {center_y} L {center_x} {center_y - radius} ' +
                       f'A {radius} {radius} 0 {large_arc_flag} 1 {end_x} {end_y} Z" ' +
                       f'fill="green" stroke="white" stroke-width="1"/>')
        else:
            # Full circle if 100%
            svg.append(f'  <circle cx="{center_x}" cy="{center_y}" r="{radius}" fill="green"/>')
    
    # Draw failure slice
    if failure_rate > 0:
        start_x = center_x + radius * math.cos(math.radians(success_angle - 90))
        start_y = center_y + radius * math.sin(math.radians(success_angle - 90))
        large_arc_flag = 1 if failure_rate > 50 else 0
        
        svg.append(f'  <path d="M {center_x} {center_y} L {start_x} {start_y} ' +
                   f'A {radius} {radius} 0 {large_arc_flag} 1 {center_x} {center_y - radius} Z" ' +
                   f'fill="red" stroke="white" stroke-width="1"/>')
    
    # Add labels
    svg.extend([
        f'  <text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="Arial" font-size="24">{success_rate:.1f}%</text>',
        f'  <text x="{width/2}" y="{height/2 + 30}" text-anchor="middle" font-family="Arial" font-size="12">{successful_runs} of {total_runs} runs</text>',
        f'  <rect x="50" y="{height-70}" width="12" height="12" fill="green"/>',
        f'  <text x="70" y="{height-60}" font-family="Arial" font-size="12">Success ({success_rate:.1f}%)</text>',
        f'  <rect x="50" y="{height-50}" width="12" height="12" fill="red"/>',
        f'  <text x="70" y="{height-40}" font-family="Arial" font-size="12">Failure ({failure_rate:.1f}%)</text>'
    ])
    
    # Close SVG
    svg.append('</svg>')
    
    # Write to file
    chart_path = os.path.join(charts_dir, "validation_success.svg")
    with open(chart_path, 'w') as f:
        f.write('\n'.join(svg))
    
    print(f"Validation chart generated: {chart_path}")

def generate_link_health_chart(summary, charts_dir):
    """Generate a link health trend chart."""
    if not summary or "link_health" not in summary or "history" not in summary["link_health"]:
        return
    
    link_history = summary["link_health"]["history"]
    if not link_history:
        # Generate a simple empty chart with message
        width, height = 800, 400
        svg = [
            f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
            f'  <rect width="{width}" height="{height}" fill="white"/>',
            f'  <text x="{width/2}" y="20" text-anchor="middle" font-family="Arial" font-size="16">Link Health Trends</text>',
            f'  <text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="Arial" font-size="14">No link health data available yet</text>',
            f'  <text x="{width/2}" y="{height/2 + 30}" text-anchor="middle" font-family="Arial" font-size="12">Data will appear after workflow runs</text>',
            '</svg>'
        ]
        
        chart_path = os.path.join(charts_dir, "link_health.svg")
        with open(chart_path, 'w') as f:
            f.write('\n'.join(svg))
        
        print(f"Empty link health chart generated: {chart_path}")
        return
    
    # Limit to last 20 entries for readability
    link_history = link_history[-20:]
    
    # Create a simple SVG stacked bar chart
    width, height = 800, 400
    padding = 50
    chart_width = width - 2 * padding
    chart_height = height - 2 * padding
    bar_width = chart_width / len(link_history) - 5
    
    # Generate SVG content
    svg = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
        f'  <rect width="{width}" height="{height}" fill="white"/>',
        f'  <text x="{width/2}" y="20" text-anchor="middle" font-family="Arial" font-size="16">Link Health Trends</text>',
        f'  <text x="{width/2}" y="{height-10}" text-anchor="middle" font-family="Arial" font-size="12">Last {len(link_history)} Runs</text>',
        f'  <text x="10" y="{height/2}" text-anchor="middle" font-family="Arial" font-size="12" transform="rotate(-90 10,{height/2})">Link Count</text>'
    ]
    
    # Draw axes
    svg.extend([
        f'  <line x1="{padding}" y1="{height-padding}" x2="{width-padding}" y2="{height-padding}" stroke="black" stroke-width="1"/>',
        f'  <line x1="{padding}" y1="{padding}" x2="{padding}" y2="{height-padding}" stroke="black" stroke-width="1"/>'
    ])
    
    # Find max total links for scaling
    max_total = max((entry["total_links"] for entry in link_history), default=1)
    
    # Draw bars
    for i, entry in enumerate(link_history):
        x = padding + i * (bar_width + 5)
        
        # Calculate heights proportionally
        total = entry["total_links"]
        if total == 0:
            continue
            
        scale_factor = chart_height / max_total
        
        success_height = entry["successful_links"] * scale_factor
        broken_height = entry["broken_links"] * scale_factor
        excluded_height = entry["excluded_links"] * scale_factor
        
        # Draw stacked bar segments
        current_y = height - padding
        
        # Success (bottom)
        if success_height > 0:
            svg.append(f'  <rect x="{x}" y="{current_y - success_height}" width="{bar_width}" height="{success_height}" fill="green"/>')
            current_y -= success_height
        
        # Excluded (middle)
        if excluded_height > 0:
            svg.append(f'  <rect x="{x}" y="{current_y - excluded_height}" width="{bar_width}" height="{excluded_height}" fill="gray"/>')
            current_y -= excluded_height
        
        # Broken (top)
        if broken_height > 0:
            svg.append(f'  <rect x="{x}" y="{current_y - broken_height}" width="{bar_width}" height="{broken_height}" fill="red"/>')
    
    # Add legend
    svg.extend([
        f'  <rect x="{width-150}" y="40" width="12" height="12" fill="green"/>',
        f'  <text x="{width-130}" y="50" font-family="Arial" font-size="12">Successful</text>',
        f'  <rect x="{width-150}" y="60" width="12" height="12" fill="gray"/>',
        f'  <text x="{width-130}" y="70" font-family="Arial" font-size="12">Excluded</text>',
        f'  <rect x="{width-150}" y="80" width="12" height="12" fill="red"/>',
        f'  <text x="{width-130}" y="90" font-family="Arial" font-size="12">Broken</text>'
    ])
    
    # Close SVG
    svg.append('</svg>')
    
    # Write to file
    chart_path = os.path.join(charts_dir, "link_health.svg")
    with open(chart_path, 'w') as f:
        f.write('\n'.join(svg))
    
    print(f"Link health chart generated: {chart_path}")

def generate_index_html(charts_dir):
    """Generate an index.html file to view all charts."""
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Awesome-Virome Repository Metrics</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 0; padding: 20px; }
            h1 { color: #333; }
            .chart-container { margin-bottom: 30px; }
            .chart { max-width: 100%; border: 1px solid #ddd; }
            .metrics-info { margin-top: 10px; color: #666; }
        </style>
    </head>
    <body>
        <h1>Awesome-Virome Repository Metrics</h1>
        <p>This page displays historical metrics for the Awesome-Virome repository.</p>
        
        <div class="chart-container">
            <h2>Performance Trends</h2>
            <img src="performance_trend.svg" alt="Performance Trends" class="chart">
            <div class="metrics-info">
                Shows the time taken for key operations over time. Lower is better.
            </div>
        </div>
        
        <div class="chart-container">
            <h2>Validation Success Rate</h2>
            <img src="validation_success.svg" alt="Validation Success Rate" class="chart">
            <div class="metrics-info">
                Shows the percentage of validation checks that passed.
            </div>
        </div>
        
        <div class="chart-container">
            <h2>Link Health</h2>
            <img src="link_health.svg" alt="Link Health" class="chart">
            <div class="metrics-info">
                Shows the status of links in the repository over time.
            </div>
        </div>
        
        <footer>
            <p>Last updated: <span id="update-time"></span></p>
            <script>
                document.getElementById('update-time').textContent = new Date().toISOString();
            </script>
        </footer>
    </body>
    </html>
    """
    
    index_path = os.path.join(charts_dir, "index.html")
    with open(index_path, 'w') as f:
        f.write(html_content)
    
    print(f"Index HTML generated: {index_path}")

def main():
    # Ensure output directory exists
    charts_dir = ensure_output_dir()
    
    # Load summary data
    summary = load_summary()
    if not summary:
        print("Error: Could not load metrics summary")
        return 1
    
    # Import math module here for pie chart calculations
    global math
    import math
    
    # Generate charts
    generate_performance_chart(summary, charts_dir)
    generate_validation_chart(summary, charts_dir)
    generate_link_health_chart(summary, charts_dir)
    
    # Generate index.html
    generate_index_html(charts_dir)
    
    print("Metrics visualization charts generated successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())