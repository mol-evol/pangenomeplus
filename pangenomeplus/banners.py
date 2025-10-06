"""Professional banners for PanGenomePlus CLI."""

from rich.console import Console
from rich.panel import Panel

from pangenomeplus import __version__


def print_startup_banner(update_available: str = None) -> None:
    """Display clean startup banner with attribution and citation.

    Args:
        update_available: If set, version string of newer release (e.g., "0.2.0")
    """
    banner_text = f"""[bold cyan]PanGenomePlus v{__version__}[/bold cyan]
[dim]Comprehensive Pangenome Analysis Pipeline[/dim]

[bold]Author:[/bold] James McInerney
[bold]GitHub:[/bold] https://github.com/mol-evol/pangenomeplus

[bold]Citation:[/bold] McInerney JO (2025). PanGenomePlus:
Comprehensive pangenome analysis for all genomic features.
GitHub: https://github.com/mol-evol/pangenomeplus"""

    if update_available:
        banner_text += f"""

[bold yellow]⚠️  Update available:[/bold yellow] v{update_available} (you have v{__version__})
[dim]Update: pip install --upgrade pangenomeplus
    or: conda update pangenomeplus[/dim]"""

    console = Console()
    console.print(Panel(banner_text, width=70, border_style="cyan"))
