#!/usr/bin/env python3
"""
Citation Report Generator for Awesome-Virome Repository

This script generates comprehensive citation reports for the tools in the
Awesome-Virome collection, including:
- Most-cited tools analysis
- Citation trends over time
- Citation network visualization
- Field impact analysis

Usage:
    python citation_report.py [--output OUTPUT] [--format FORMAT]
"""

import os
import json
import argparse
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Set
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Directory constants
METADATA_DIR = os.path.join("metadata", "academic_impact")
REPORTS_DIR = os.path.join("reports", "citations")
ACADEMIC_IMPACT_FILE = os.path.join(METADATA_DIR, "academic_impact.json")
SUMMARY_FILE = os.path.join(METADATA_DIR, "summary.json")


class CitationReportGenerator:
    """Generates citation reports for Awesome-Virome tools."""
    
    def __init__(self, output_dir: str = REPORTS_DIR):
        """Initialize the citation report generator."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Load academic impact data
        self.academic_impact = {}
        if os.path.exists(ACADEMIC_IMPACT_FILE):
            try:
                with open(ACADEMIC_IMPACT_FILE, 'r') as f:
                    self.academic_impact = json.load(f)
            except json.JSONDecodeError:
                logger.error(f"Could not parse academic impact data from {ACADEMIC_IMPACT_FILE}")
        
        # Load summary data
        self.summary = {}
        if os.path.exists(SUMMARY_FILE):
            try:
                with open(SUMMARY_FILE, 'r') as f:
                    self.summary = json.load(f)
            except json.JSONDecodeError:
                logger.error(f"Could not parse summary data from {SUMMARY_FILE}")
    
    def generate_most_cited_report(self) -> Dict[str, Any]:
        """Generate a report of the most-cited tools."""
        logger.info("Generating most-cited tools report")
        
        # Extract citation data for all tools
        citation_data = []
        for key, data in self.academic_impact.items():
            name = data.get('name', '')
            url = data.get('url', '')
            doi = data.get('doi', '')
            total_citations = data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0)
            influential_citations = data.get('citation_metrics', {}).get('metrics', {}).get('influential_citations', 0)
            
            if total_citations > 0:
                citation_data.append({
                    'name': name,
                    'url': url,
                    'doi': doi,
                    'total_citations': total_citations,
                    'influential_citations': influential_citations
                })
        
        # Sort by total citations
        most_cited = sorted(citation_data, key=lambda x: x['total_citations'], reverse=True)
        
        # Generate report
        report = {
            'title': 'Most-Cited Tools in Awesome-Virome Repository',
            'generated': datetime.now().isoformat(),
            'total_tools_with_citations': len(citation_data),
            'most_cited_tools': most_cited[:20],  # Top 20 most cited
            'total_citations_all_tools': sum(tool['total_citations'] for tool in citation_data),
            'average_citations_per_tool': sum(tool['total_citations'] for tool in citation_data) / len(citation_data) if citation_data else 0
        }
        
        # Save report
        report_file = os.path.join(self.output_dir, "most_cited_report.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Generate visualization
        if most_cited:
            self._plot_most_cited(most_cited[:15])  # Top 15 for visualization
        
        return report
    
    def _plot_most_cited(self, most_cited: List[Dict[str, Any]]):
        """Create a bar chart of the most-cited tools."""
        plt.figure(figsize=(12, 8))
        names = [tool['name'] for tool in most_cited]
        citations = [tool['total_citations'] for tool in most_cited]
        influential = [tool['influential_citations'] for tool in most_cited]
        
        plt.bar(names, citations, color='steelblue', label='Total Citations')
        plt.bar(names, influential, color='firebrick', label='Influential Citations')
        
        plt.title('Most-Cited Bioinformatics Tools in Awesome-Virome Repository', fontsize=16)
        plt.xlabel('Tool Name', fontsize=14)
        plt.ylabel('Citation Count', fontsize=14)
        plt.xticks(rotation=45, ha='right')
        plt.legend()
        plt.tight_layout()
        
        # Save plot
        plt.savefig(os.path.join(self.output_dir, "most_cited_tools.png"), dpi=300)
        plt.savefig(os.path.join(self.output_dir, "most_cited_tools.svg"), format='svg')
        plt.close()
    
    def generate_citation_trends(self) -> Dict[str, Any]:
        """Generate a report of citation trends over time."""
        logger.info("Generating citation trends report")
        
        # Collect citation data by year for all tools
        yearly_data = defaultdict(int)
        tool_yearly_data = {}
        
        for key, data in self.academic_impact.items():
            name = data.get('name', '')
            citations_by_year = data.get('citation_metrics', {}).get('metrics', {}).get('citations_by_year', {})
            
            if citations_by_year:
                tool_yearly_data[name] = citations_by_year
                for year, count in citations_by_year.items():
                    yearly_data[year] += count
        
        # Sort years
        sorted_years = sorted(yearly_data.keys())
        yearly_totals = [yearly_data[year] for year in sorted_years]
        
        # Find top tools to include in the trend visualization
        total_citations_by_tool = {}
        for name, yearly in tool_yearly_data.items():
            total_citations_by_tool[name] = sum(yearly.values())
        
        top_tools = sorted(total_citations_by_tool.items(), key=lambda x: x[1], reverse=True)[:10]
        top_tool_names = [name for name, _ in top_tools]
        
        # Generate report
        report = {
            'title': 'Citation Trends Over Time for Awesome-Virome Tools',
            'generated': datetime.now().isoformat(),
            'years': sorted_years,
            'yearly_totals': yearly_totals,
            'top_tools': top_tools,
            'tool_yearly_data': {name: tool_yearly_data[name] for name in top_tool_names if name in tool_yearly_data}
        }
        
        # Save report
        report_file = os.path.join(self.output_dir, "citation_trends_report.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Generate visualizations
        self._plot_citation_trends(sorted_years, yearly_totals, tool_yearly_data, top_tool_names)
        
        return report
    
    def _plot_citation_trends(self, years, totals, tool_data, top_tools):
        """Create visualizations of citation trends over time."""
        # Plot overall citation trend
        plt.figure(figsize=(12, 6))
        plt.plot(years, totals, 'o-', linewidth=2, markersize=8, color='steelblue')
        
        plt.title('Citation Trends Over Time for All Awesome-Virome Tools', fontsize=16)
        plt.xlabel('Year', fontsize=14)
        plt.ylabel('Number of Citations', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save overall trend plot
        plt.savefig(os.path.join(self.output_dir, "citation_trend_overall.png"), dpi=300)
        plt.savefig(os.path.join(self.output_dir, "citation_trend_overall.svg"), format='svg')
        plt.close()
        
        # Plot individual tool trends
        plt.figure(figsize=(14, 8))
        
        # Color map for different tools
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                 '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        # Plot each top tool
        for idx, tool_name in enumerate(top_tools):
            if tool_name in tool_data:
                tool_years = sorted(tool_data[tool_name].keys())
                tool_citations = [tool_data[tool_name][year] for year in tool_years]
                
                plt.plot(tool_years, tool_citations, 'o-', linewidth=2, 
                         label=tool_name, color=colors[idx % len(colors)])
        
        plt.title('Citation Trends for Top 10 Most-Cited Tools', fontsize=16)
        plt.xlabel('Year', fontsize=14)
        plt.ylabel('Number of Citations', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        
        # Save tool trends plot
        plt.savefig(os.path.join(self.output_dir, "citation_trend_top_tools.png"), dpi=300)
        plt.savefig(os.path.join(self.output_dir, "citation_trend_top_tools.svg"), format='svg')
        plt.close()
    
    def generate_citation_network(self) -> Dict[str, Any]:
        """Generate a citation network visualization and report."""
        logger.info("Generating citation network report")
        
        # Create a graph
        G = nx.DiGraph()
        
        # Collect tools with DOIs and related papers
        tools_with_relations = []
        
        for key, data in self.academic_impact.items():
            name = data.get('name', '')
            doi = data.get('doi', '')
            related_papers = data.get('related_papers', [])
            
            if doi and related_papers:
                # Add tool node
                G.add_node(name, type='tool', citations=data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0))
                
                # Add related paper nodes and edges
                for paper in related_papers:
                    paper_title = paper.get('title', 'Unknown')
                    paper_id = paper.get('paperId', '')
                    
                    if paper_id:
                        # Add paper node
                        G.add_node(paper_id, type='paper', title=paper_title, 
                                   citations=paper.get('citationCount', 0),
                                   year=paper.get('year', ''))
                        
                        # Add edge from tool to paper
                        G.add_edge(name, paper_id, type='related')
                
                tools_with_relations.append({
                    'name': name,
                    'doi': doi,
                    'related_count': len(related_papers)
                })
        
        # Find connections between tools via related papers
        connections = []
        for tool1 in G.nodes():
            if G.nodes[tool1].get('type') == 'tool':
                # Get related papers for this tool
                tool1_papers = [n for n in G.successors(tool1)]
                
                for tool2 in G.nodes():
                    if G.nodes[tool2].get('type') == 'tool' and tool1 != tool2:
                        # Get related papers for tool2
                        tool2_papers = [n for n in G.successors(tool2)]
                        
                        # Find common papers
                        common_papers = set(tool1_papers).intersection(set(tool2_papers))
                        
                        if common_papers:
                            connections.append({
                                'tool1': tool1,
                                'tool2': tool2,
                                'common_papers': list(common_papers),
                                'connection_strength': len(common_papers)
                            })
                            
                            # Add edge between tools
                            G.add_edge(tool1, tool2, weight=len(common_papers), 
                                     papers=list(common_papers), type='tool_tool')
        
        # Generate report
        report = {
            'title': 'Citation Network Analysis for Awesome-Virome Tools',
            'generated': datetime.now().isoformat(),
            'tools_in_network': len([n for n, attr in G.nodes(data=True) if attr.get('type') == 'tool']),
            'papers_in_network': len([n for n, attr in G.nodes(data=True) if attr.get('type') == 'paper']),
            'tool_connections': connections,
            'most_connected_tools': sorted(
                [n for n, attr in G.nodes(data=True) if attr.get('type') == 'tool'],
                key=lambda x: G.degree(x),
                reverse=True
            )[:10]
        }
        
        # Save report
        report_file = os.path.join(self.output_dir, "citation_network_report.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Generate network visualization if we have enough connections
        if len(connections) > 1:
            self._plot_citation_network(G)
        
        return report
    
    def _plot_citation_network(self, G):
        """Create a visualization of the citation network."""
        # Create a subgraph of only tools and their connections
        tool_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'tool']
        tool_graph = G.subgraph(tool_nodes).copy()
        
        # Remove isolated nodes
        for node in list(tool_graph.nodes()):
            if tool_graph.degree(node) == 0:
                tool_graph.remove_node(node)
        
        if len(tool_graph.nodes()) > 1:
            plt.figure(figsize=(14, 12))
            
            # Node positions using spring layout
            pos = nx.spring_layout(tool_graph, seed=42)
            
            # Node sizes based on citation count
            node_sizes = [G.nodes[n].get('citations', 10) * 2 + 100 for n in tool_graph.nodes()]
            
            # Edge weights
            edge_weights = [G[u][v].get('weight', 1) for u, v in tool_graph.edges()]
            
            # Draw the network
            nx.draw_networkx_nodes(tool_graph, pos, node_size=node_sizes, 
                                  node_color='skyblue', alpha=0.8)
            nx.draw_networkx_edges(tool_graph, pos, width=edge_weights, 
                                  edge_color='gray', alpha=0.6, arrows=True)
            nx.draw_networkx_labels(tool_graph, pos, font_size=10)
            
            plt.title('Citation Connections Between Awesome-Virome Tools', fontsize=16)
            plt.axis('off')
            plt.tight_layout()
            
            # Save network visualization
            plt.savefig(os.path.join(self.output_dir, "citation_network.png"), dpi=300)
            plt.savefig(os.path.join(self.output_dir, "citation_network.svg"), format='svg')
            plt.close()
    
    def generate_all_reports(self) -> Dict[str, Any]:
        """Generate all citation reports."""
        logger.info("Generating all citation reports")
        
        # Check if we have academic impact data
        if not self.academic_impact:
            logger.error("No academic impact data found. Run academic_impact.py first.")
            return {
                'error': 'No academic impact data found',
                'status': 'failed'
            }
        
        # Generate individual reports
        most_cited_report = self.generate_most_cited_report()
        trends_report = self.generate_citation_trends()
        network_report = self.generate_citation_network()
        
        # Generate an overall summary
        summary = {
            'title': 'Academic Impact Analysis for Awesome-Virome Repository',
            'generated': datetime.now().isoformat(),
            'total_tools_analyzed': len(self.academic_impact),
            'total_tools_with_doi': sum(1 for data in self.academic_impact.values() if data.get('doi')),
            'total_tools_with_citations': sum(1 for data in self.academic_impact.values() 
                                             if data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0) > 0),
            'total_citations': sum(data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0) 
                                  for data in self.academic_impact.values()),
            'citation_reports': {
                'most_cited': most_cited_report.get('title'),
                'trends': trends_report.get('title'),
                'network': network_report.get('title')
            },
            'status': 'success'
        }
        
        # Save overall summary
        summary_file = os.path.join(self.output_dir, "citation_reports_summary.json")
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Generate an HTML report
        self._generate_html_report(summary, most_cited_report, trends_report, network_report)
        
        return summary
    
    def _generate_html_report(self, summary, most_cited, trends, network):
        """Generate an HTML report combining all citation data."""
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Awesome-Virome Citation Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3 {{
            color: #0066cc;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .summary-box {{
            background: #f5f5f5;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .citation-report {{
            margin-bottom: 40px;
        }}
        .citation-chart {{
            text-align: center;
            margin: 20px 0;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px 15px;
            border-bottom: 1px solid #ddd;
            text-align: left;
        }}
        th {{
            background-color: #0066cc;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .footer {{
            margin-top: 40px;
            text-align: center;
            font-size: 0.9em;
            color: #666;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Awesome-Virome Academic Impact Report</h1>
        <p>Generated on {datetime.now().strftime('%B %d, %Y')}</p>
        
        <div class="summary-box">
            <h2>Summary</h2>
            <p>This report provides an analysis of academic impact for tools in the Awesome-Virome repository.</p>
            <ul>
                <li><strong>Total Tools Analyzed:</strong> {summary['total_tools_analyzed']}</li>
                <li><strong>Tools with DOIs:</strong> {summary['total_tools_with_doi']}</li>
                <li><strong>Tools with Citations:</strong> {summary['total_tools_with_citations']}</li>
                <li><strong>Total Citations:</strong> {summary['total_citations']}</li>
            </ul>
        </div>
        
        <div class="citation-report">
            <h2>Most-Cited Tools</h2>
            <p>The following tools have received the most citations in academic literature.</p>
            
            <div class="citation-chart">
                <img src="most_cited_tools.png" alt="Most-Cited Tools" style="max-width: 100%;">
            </div>
            
            <h3>Top 10 Most-Cited Tools</h3>
            <table>
                <tr>
                    <th>Rank</th>
                    <th>Tool Name</th>
                    <th>Citations</th>
                    <th>Influential Citations</th>
                    <th>DOI</th>
                </tr>
"""

        # Add most-cited tools table rows
        for i, tool in enumerate(most_cited.get('most_cited_tools', [])[:10], 1):
            html_content += f"""
                <tr>
                    <td>{i}</td>
                    <td>{tool.get('name', '')}</td>
                    <td>{tool.get('total_citations', 0)}</td>
                    <td>{tool.get('influential_citations', 0)}</td>
                    <td>{tool.get('doi', '-')}</td>
                </tr>"""

        html_content += """
            </table>
        </div>
        
        <div class="citation-report">
            <h2>Citation Trends Over Time</h2>
            <p>This section shows how citations to Awesome-Virome tools have evolved over time.</p>
            
            <div class="citation-chart">
                <img src="citation_trend_overall.png" alt="Overall Citation Trend" style="max-width: 100%;">
            </div>
            
            <div class="citation-chart">
                <img src="citation_trend_top_tools.png" alt="Citation Trends for Top Tools" style="max-width: 100%;">
            </div>
        </div>
        
        <div class="citation-report">
            <h2>Citation Network</h2>
            <p>This visualization shows connections between tools based on shared citations and related papers.</p>
            
            <div class="citation-chart">
                <img src="citation_network.png" alt="Citation Network" style="max-width: 100%;">
            </div>
            
            <h3>Most Connected Tools</h3>
            <p>These tools share the most citations and related papers with other tools in the repository.</p>
            <ul>
"""

        # Add most connected tools
        for tool in network.get('most_connected_tools', []):
            html_content += f"""                <li>{tool}</li>\n"""

        html_content += """
            </ul>
        </div>
        
        <div class="footer">
            <p>This report was automatically generated by the Awesome-Virome citation analysis system.</p>
        </div>
    </div>
</body>
</html>
"""

        # Save HTML report
        with open(os.path.join(self.output_dir, "citation_report.html"), 'w') as f:
            f.write(html_content)


def main():
    """Main function to run the citation report generation."""
    parser = argparse.ArgumentParser(description='Generate citation reports for Awesome-Virome tools')
    parser.add_argument('--output', help='Output directory for reports', default=REPORTS_DIR)
    parser.add_argument('--format', help='Report format (json, html, or all)', choices=['json', 'html', 'all'], default='all')
    args = parser.parse_args()
    
    # Update output directory if specified
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize report generator
    generator = CitationReportGenerator(output_dir=output_dir)
    
    # Generate all reports
    summary = generator.generate_all_reports()
    
    if summary.get('status') == 'success':
        logger.info(f"Citation reports generated successfully in {output_dir}")
        logger.info(f"Total tools analyzed: {summary['total_tools_analyzed']}")
        logger.info(f"Tools with DOIs: {summary['total_tools_with_doi']}")
        logger.info(f"Tools with citations: {summary['total_tools_with_citations']}")
        logger.info(f"Total citations: {summary['total_citations']}")
    else:
        logger.error(f"Failed to generate citation reports: {summary.get('error', 'Unknown error')}")


if __name__ == "__main__":
    main()