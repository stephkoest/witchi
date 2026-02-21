from setuptools import setup, find_packages
from pathlib import Path
import re


def read_version():
    init_path = Path("src/witchi/__init__.py")
    version_pattern = re.compile(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]')
    with init_path.open("r") as f:
        for line in f:
            match = version_pattern.match(line)
            if match:
                return match.group(1)
    raise RuntimeError("Version string not found.")


setup(
    name="witchi",
    version=read_version(),
    description="A compositional bias pruning tool for multiple sequence alignments.",
    author="Stephan Koestlbacher",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            "witchi=witchi.witchi:main",
        ],
    },
    install_requires=["numpy", "biopython", "joblib"],
)
