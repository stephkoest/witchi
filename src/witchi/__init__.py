import os
import subprocess

__version__ = "0.2.0rc1"


def _git_sha():
    """Return short git SHA when running from a checkout, else None."""
    try:
        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=pkg_dir,
            capture_output=True,
            text=True,
            timeout=2,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return None


def version_banner():
    """One-line provenance banner: 'witchi <version> [<sha>]'."""
    sha = _git_sha()
    return f"witchi {__version__}" + (f" ({sha})" if sha else "")
