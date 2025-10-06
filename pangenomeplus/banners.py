"""Professional banners for PanGenomePlus CLI."""

from rich.console import Console
from rich.table import Table

from pangenomeplus import __version__


def print_startup_banner(update_available: str = None) -> None:
    """Display clean startup banner with attribution and citation.

    Args:
        update_available: If set, version string of newer release (e.g., "0.2.0")
    """
    console = Console()

    # Create table with fixed width for perfect right edge alignment
    table = Table(show_header=False, box=None, padding=(0, 2), width=70, border_style="cyan")
    table.add_column(width=66)  # Content width (70 - 4 for borders/padding)

    # Add banner content
    table.add_row(f"[bold cyan]PanGenomePlus v{__version__}[/bold cyan]")
    table.add_row("[dim]Comprehensive Pangenome Analysis Pipeline[/dim]")
    table.add_row("")
    table.add_row("[bold]Author:[/bold] James McInerney")
    table.add_row("[bold]GitHub:[/bold] https://github.com/mol-evol/pangenomeplus")
    table.add_row("")
    table.add_row("[bold]Citation:[/bold] McInerney JO (2025). PanGenomePlus:")
    table.add_row("Comprehensive pangenome analysis for all genomic features.")
    table.add_row("GitHub: https://github.com/mol-evol/pangenomeplus")

    if update_available:
        table.add_row("")
        table.add_row(f"[bold yellow]⚠️  Update available:[/bold yellow] v{update_available} (you have v{__version__})")
        table.add_row("[dim]Update: pip install --upgrade pangenomeplus[/dim]")
        table.add_row("[dim]    or: conda update pangenomeplus[/dim]")

    # Print with box border
    from rich.box import ROUNDED
    bordered_table = Table(show_header=False, box=ROUNDED, border_style="cyan", padding=0)
    bordered_table.add_column(width=70)
    bordered_table.add_row(table)

    console.print(bordered_table)
